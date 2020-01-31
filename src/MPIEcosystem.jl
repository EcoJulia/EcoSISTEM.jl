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
    MPIEcosystem{Part <: AbstractAbiotic} <:
       AbstractEcosystem{Part, SL, TR}

MPIEcosystem houses information on species and their interaction with their
environment. For species, it holds abundances and locations, `abundances`,
as well as properties such as trait information, `spplist`, and movement types,
`lookup`. For environments, it provides information on environmental conditions
and available resources,`abenv`. Finally, there is a slot for the relationship
between the environment and the characteristics of the species, `relationship`.
"""
mutable struct MPIEcosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
    TR <: AbstractTraitRelationship} <: AbstractEcosystem{Part, SL, TR}
  abundances::GridLandscape
  spplist::SL
  abenv::Part
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::TR
  lookup::Vector{Lookup}
  counts::Vector{Int32}
  firstspecies::Int64
  rank::Int64
  comm::MPI.Comm
  cache::Cache

  function MPIEcosystem{Part, SL, TR}(abundances::GridLandscape,
    spplist::SL, abenv::Part, ordinariness::Union{Matrix{Float64}, Missing},
    relationship::TR, lookup::Vector{Lookup}, counts::Vector, firstspecies::Int64, rank::Int64,
    comm::MPI.Comm, cache::Cache) where {Part <:
     AbstractAbiotic,
    SL <: SpeciesList, TR <: AbstractTraitRelationship}
    tematch(spplist, abenv) || error("Traits do not match habitats")
    trmatch(spplist, relationship) || error("Traits do not match trait functions")
    #_mcmatch(abundances.matrix, spplist, abenv) ||
    #  error("Dimension mismatch")
    new{Part, SL, TR}(abundances, spplist, abenv, ordinariness, relationship, lookup, counts, firstspecies, rank, comm, cache)
  end
end

"""
    MPIEcosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
        rel::AbstractTraitRelationship)

Function to create an `MPIEcosystem` given a species list, an abiotic environment
and trait relationship.
"""
function MPIEcosystem(popfun::Function, spplist::SpeciesList{T, Req}, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship, comm::MPI.Comm = MPI.COMM_WORLD) where {T, Req}
    rank = MPI.Comm_rank(comm)
    totalsize = MPI.Comm_size(comm)
    numspecies = length(spplist.names)
    count = div(numspecies + totalsize - 1, totalsize)
    counts = Int32.(fill(count, totalsize))
    counts[end] = numspecies - sum(counts) + count
    # Create matrix landscape of zero abundances
    ml = emptygridlandscape(abenv, spplist)
    # Populate this matrix with species abundances
    if rank == 0
        popfun(ml, spplist, abenv, rel)
    end
    MPI.Bcast!(ml.matrix, 0, comm)
    # Create lookup table of all moves and their probabilities
    indices = vcat([0], cumsum(counts))
    firstspecies = indices[rank + 1] + 1
    rankspp = firstspecies : indices[rank + 2]
    lookup_tab = collect(map(k -> genlookups(abenv.habitat, k), @view getkernels(spplist.movement)[rankspp]))
    nm = zeros(Int64, size(ml.matrix))
    totalE = zeros(Float64, (size(ml.matrix, 2), numrequirements(Req)))
    MPIEcosystem{typeof(abenv), typeof(spplist), typeof(rel)}(ml, spplist, abenv,
    missing, rel, lookup_tab, counts, firstspecies, rank, comm, Cache(nm, totalE, false))
end

function MPIEcosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship)
   return MPIEcosystem(populate!, spplist, abenv, rel)
end
