using MPI
using Diversity
if VERSION > v"1.0.0"
    using HCubature
else
    using Cubature
end
using Unitful
using MyUnitful
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
function MPIEcosystem(popfun::Function, spplist::SpeciesList{T, Req},
  abenv::GridAbioticEnv, rel::AbstractTraitRelationship) where {T, Req}
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
