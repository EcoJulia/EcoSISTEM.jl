using EcoSISTEM
using MPI
using Diversity
using HCubature
using Unitful
using EcoSISTEM.Units
using Missings

import Diversity: _calcabundance

function EcoSISTEM.MPIEcosystem(abundances::MPIGL,
    spplist::SL, abenv::Part, ordinariness::Union{Matrix{Float64},
    Missing}, relationship::TR, lookup::Vector{EcoSISTEM.Lookup},
    sppcounts::Vector, firstsp::Int64, sccounts::Vector,
    firstsc::Int64, cache::EcoSISTEM.Cache) where {MPIGL <: MPIGridLandscape,
    Part <: EcoSISTEM.AbstractAbiotic, SL <: SpeciesList,
    TR <: EcoSISTEM.AbstractTraitRelationship}
    tematch(spplist, abenv) || error("Traits do not match habitats")
    trmatch(spplist, relationship) || error("Traits do not match trait functions")
    #_mcmatch(abundances.matrix, spplist, abenv) ||
    #  error("Dimension mismatch")
    return MPIEcosystem{MPIGL, Part, SL, TR}(abundances, spplist, abenv, ordinariness, relationship, lookup, sppcounts, firstsp, sccounts, firstsc,
    cache)
end

"""
    MPIEcosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
        rel::AbstractTraitRelationship)

Function to create an `MPIEcosystem` given a species list, an abiotic environment and trait relationship.
"""
function EcoSISTEM.MPIEcosystem(popfun::F, spplist::SpeciesList{SpeciesTypes{T, Req, MO, TY}, P, PA},
  abenv::GridAbioticEnv, rel::EcoSISTEM.AbstractTraitRelationship
  ) where {F<:Function, T, Req, MO, TY, P, PA}
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    totalsize = MPI.Comm_size(comm)
    numspp = length(spplist.species.names)
    numsc = countsubcommunities(abenv.habitat)

    count = div(numspp, totalsize)
    sppcounts = Int32.(fill(count, totalsize))
    sppcounts[1:(numspp - sum(sppcounts))] .+= 1
    sppindices = vcat([0], cumsum(sppcounts))
    firstsp = sppindices[rank + 1] + 1

    sccount = div(numsc, totalsize)
    sccounts = Int32.(fill(sccount, totalsize))
    sccounts[1:(numsc - sum(sccounts))] .+= 1
    scindices = vcat([0], cumsum(sccounts))
    firstsc = scindices[rank + 1] + 1

    # Create matrix landscape of zero abundances
    ml = emptyMPIgridlandscape(sppcounts, sccounts)

    # Populate this matrix with species abundances
    popfun(ml, spplist, abenv, rel)

    rankspp = firstsp : sppindices[rank + 2]
    lookup_tab = collect(map(k -> EcoSISTEM.genlookups(abenv.habitat, k), @view getkernels(spplist.species.movement)[rankspp]))
    nm = zeros(Int64, (sppcounts[rank + 1], numsc))
    totalE = zeros(Float64, (numsc, EcoSISTEM.numrequirements(Req)))
    MPIEcosystem{typeof(ml), typeof(abenv),
    typeof(spplist), typeof(rel)}(ml, spplist, abenv, missing, rel,
    lookup_tab, sppcounts, firstsp, sccounts, firstsc,
    EcoSISTEM.Cache(nm, totalE, false))
end

function EcoSISTEM.MPIEcosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
   rel::EcoSISTEM.AbstractTraitRelationship)
   return MPIEcosystem(populate!, spplist, abenv, rel)
end

"""
    gather_abundance(eco::MPIEcosystem)

Gather full abundances matrix on root node.
"""
function EcoSISTEM.gather_abundance(eco::MPIEcosystem)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    true_abuns = zeros(Int64, counttypes(eco), countsubcommunities(eco))
    if rank == 0
        output_vbuf = VBuffer(true_abuns, Int32.(eco.sppcounts .* sum(eco.sccounts)))
    else
        output_vbuf = VBuffer(nothing)
    end
    MPI.Gatherv!(vcat(eco.abundances.reshaped_cols...)[1:end], output_vbuf, 0, comm)
    return true_abuns
end


function EcoSISTEM.gettransitions(eco::MPIEcosystem)
    return nothing
end

import Diversity.API: _getabundance
function Diversity.API._getabundance(eco::MPIEcosystem, raw::Bool)
    if raw
        return eco.abundances.rows_matrix
    else
        return _calcabundance(Diversity.API._gettypes(eco), eco.abundances.rows_matrix / sum(eco.abundances.rows_matrix))[1]
    end
end

import Diversity.API: _getmetaabundance
function Diversity.API._getmetaabundance(eco::MPIEcosystem)
    comm = MPI.COMM_WORLD
    ab = sum(_getabundance(eco), dims = 2)
    return MPI.Allgatherv(MPI.VBuffer(ab, eco.sppcounts), comm)
end

import Diversity.API: _getweight
function Diversity.API._getweight(eco::MPIEcosystem)
    comm = MPI.COMM_WORLD
    w = sum(_getabundance(eco, false), dims = 1)
    return MPI.Allreduce(w, +, comm)[1, :]
end

import Diversity.API: _getordinariness!
function Diversity.API._getordinariness!(eco::MPIEcosystem)
    if ismissing(eco.ordinariness)
        eco.ordinariness = _calcordinariness(eco)
    end
    return eco.ordinariness
end

import Diversity.API: _calcordinariness
function Diversity.API._calcordinariness(eco::MPIEcosystem)
    relab = getabundance(eco, false)
    sp_rng = eco.abundances.rows_tuple.first:eco.abundances.rows_tuple.last
    return Diversity.API._calcsimilarity(eco.spplist.species.types, one(eltype(relab)))[sp_rng, sp_rng] * relab
end

function EcoSISTEM.gather_diversity(eco::MPIEcosystem, divmeasure::F, q) where {F<:Function}
    comm = MPI.COMM_WORLD
    totalsize = MPI.Comm_size(comm)
    div = divmeasure(eco, q)
    totalabun = MPI.Gather(sum(eco.abundances.rows_matrix), 0, comm)
    mpidivs = MPI.Gather(div[!, :diversity], 0, comm)
    if rank == 0
        mpidivs = vcat(reshape(mpidivs, countsubcommunities(eco), totalsize), div[!, :q])
        div[!, :diversity] .= mapslices(x -> Diversity.powermean(x[:, 1:end-1], 1 - x[:, end], totalabun .* 1.0), mpidivs, dims = 2)[:, 1]
        return div
    else
        return div
    end
end
