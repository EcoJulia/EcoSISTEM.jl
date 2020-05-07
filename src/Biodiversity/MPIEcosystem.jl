using MPI
using Diversity
if VERSION > v"1.0.0"
    using HCubature
else
    using Cubature
end
using Unitful
using Simulation.Units
using Missings
using Compat

import Diversity: _calcabundance
"""
    MPIEcosystem{MPIGL <: MPIGridLandscape, Part <: AbstractAbiotic,
    SL <: SpeciesList, TR <: AbstractTraitRelationship} <: AbstractEcosystem{Part, SL, TR}

MPIEcosystem houses information on species and their interaction with their environment. It houses all information of a normal `Ecosystem` (see documentation for more details), with additional fields to describe which species are calculated on which machine. This includes: `sppcounts` - a vector of number of species per node, `firstsp` - the identity of the first species held by that particular node.
"""
mutable struct MPIEcosystem{MPIGL <: MPIGridLandscape, Part <: AbstractAbiotic, SL <: SpeciesList, TR <: AbstractTraitRelationship} <: AbstractEcosystem{Part, SL, TR}
  abundances::MPIGL
  spplist::SL
  abenv::Part
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::TR
  lookup::Vector{Lookup}
  sppcounts::Vector{Int32}
  firstsp::Int64
  sccounts::Vector{Int32}
  firstsc::Int64
  cache::Cache

  function MPIEcosystem{MPIGL, Part, SL, TR}(abundances::MPIGL,
    spplist::SL, abenv::Part, ordinariness::Union{Matrix{Float64},
    Missing}, relationship::TR, lookup::Vector{Lookup},
    sppcounts::Vector, firstsp::Int64, sccounts::Vector,
    firstsc::Int64, cache::Cache) where {MPIGL <: MPIGridLandscape,
    Part <: AbstractAbiotic, SL <: SpeciesList,
    TR <: AbstractTraitRelationship}
    tematch(spplist, abenv) || error("Traits do not match habitats")
    trmatch(spplist, relationship) || error("Traits do not match trait functions")
    #_mcmatch(abundances.matrix, spplist, abenv) ||
    #  error("Dimension mismatch")
    new{MPIGL, Part, SL, TR}(abundances, spplist, abenv, ordinariness, relationship, lookup, sppcounts, firstsp, sccounts, firstsc,
    cache)
  end
end

"""
    MPIEcosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
        rel::AbstractTraitRelationship)

Function to create an `MPIEcosystem` given a species list, an abiotic environment and trait relationship.
"""
function MPIEcosystem(popfun::F, spplist::SpeciesList{T, Req},
  abenv::GridAbioticEnv, rel::AbstractTraitRelationship
  ) where {F<:Function, T, Req}
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    totalsize = MPI.Comm_size(comm)
    numspp = length(spplist.names)
    numsc = countsubcommunities(abenv.habitat)

    count = div(numspp + totalsize - 1, totalsize)
    sppcounts = Int32.(fill(count, totalsize))
    sppcounts[end] = numspp - sum(sppcounts) + count
    sppindices = vcat([0], cumsum(sppcounts))
    firstsp = sppindices[rank + 1] + 1

    sccount = div(numsc + totalsize - 1, totalsize)
    sccounts = Int32.(fill(sccount, totalsize))
    sccounts[end] = numsc - sum(sccounts) + sccount
    scindices = vcat([0], cumsum(sccounts))
    firstsc = scindices[rank + 1] + 1

    # Create matrix landscape of zero abundances
    ml = emptyMPIgridlandscape(sppcounts, sccounts)

    # Populate this matrix with species abundances
    popfun(ml, spplist, abenv, rel)

    rankspp = firstsp : sppindices[rank + 2]
    lookup_tab = collect(map(k -> genlookups(abenv.habitat, k), @view getkernels(spplist.movement)[rankspp]))
    nm = zeros(Int64, (sppcounts[rank + 1], numsc))
    totalE = zeros(Float64, (numsc, numrequirements(Req)))
    MPIEcosystem{typeof(ml), typeof(abenv),
    typeof(spplist), typeof(rel)}(ml, spplist, abenv, missing, rel,
    lookup_tab, sppcounts, firstsp, sccounts, firstsc,
    Cache(nm, totalE, false))
end

function MPIEcosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship)
   return MPIEcosystem(populate!, spplist, abenv, rel)
end

# function gather_abundance!(eco::MPIEcosystem)
#     comm = MPI.COMM_WORLD
#     rank = MPI.Comm_rank(comm)
#     abun = MPI.Gatherv(eco.abundances.rows_matrix, Int32.(eco.sppcounts .* sum(eco.sccounts)), 0, comm)
#     if rank == 0
#         eco.fullabun = reshape(abun, counttypes(eco), countsubcommunities(eco))
#     else
#         eco.fullabun = missing
#     end
# end

import Diversity.API: _getabundance
function _getabundance(eco::MPIEcosystem, raw::Bool)
    if raw
        return eco.abundances.rows_matrix
    else
        return _calcabundance(_gettypes(eco), eco.abundances.rows_matrix / sum(eco.abundances.rows_matrix))[1]
    end
end

import Diversity.API: _getmetaabundance
function _getmetaabundance(eco::MPIEcosystem)
    comm = MPI.COMM_WORLD
    ab = sum(_getabundance(eco), dims = 2)
    return MPI.Allgatherv(ab, eco.sppcounts, comm)
end

import Diversity.API: _getweight
function _getweight(eco::MPIEcosystem)
    comm = MPI.COMM_WORLD
    w = sum(_getabundance(eco, false), dims = 1)
    return MPI.Allreduce(w, +, comm)[1, :]
end

import Diversity.API: _getordinariness!
function _getordinariness!(eco::MPIEcosystem)
    if ismissing(eco.ordinariness)
        eco.ordinariness = _calcordinariness(eco)
    end
    return eco.ordinariness
end

import Diversity.API: _calcordinariness
function _calcordinariness(eco::MPIEcosystem)
    relab = getabundance(eco, false)
    sp_rng = eco.abundances.rows_tuple.first:eco.abundances.rows_tuple.last
    return _calcsimilarity(eco.spplist.types, one(eltype(relab)))[sp_rng, sp_rng] * relab
end

function gather_diversity(eco::MPIEcosystem, divmeasure::F, q) where {F<:Function}
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
