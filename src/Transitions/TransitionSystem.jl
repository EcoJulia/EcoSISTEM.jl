using Diversity
using DataFrames
using Unitful
using EcoSISTEM.Units
using Missings
using RecipesBase


import Diversity: _calcabundance

"""
    AbstractCache

Abstract type for `Ecosystem` caches.
"""
abstract type AbstractCache end

"""
    Cache

Cache houses an integer array of moves made by all species in a timestep for the
update! function, `netmigration`, a matrix of current resource used in the
Ecosystem, `totalE`, and a Bool to say if these caches are `valid`.
"""
mutable struct Cache <: AbstractCache
  netmigration::Array{Int64, 2}
  totalE::Matrix{Float64}
  valid::Bool
end

"""
    EpiCache

EpiCache houses a float array of moves made by all force of infection categories in a timestep for the
update! function, `virusmigration`, an integer of how many initially infected individuals
are introduced to the system, an array of locations that are active in the system, `ordered_active`,
and a Bool to say if these caches are `valid`.
"""
mutable struct EpiCache <: AbstractCache
  virusmigration::Array{Float64, 2}
  initial_infected::Int64
  ordered_active::Vector{Int64}
  valid::Bool
end

"""
    PlantCache

Cache houses an integer array of seed production made by all species in the Ecosystem,
`seedbank`, moves made by all dispersing seeds in a timestep, `netmigration`,
a matrix of current resource used in the Ecosystem, `totalE`,
and a Bool to say if these caches are `valid`.
"""
mutable struct PlantCache <: AbstractCache
  netmigration::Array{Int64, 2}
  seedbank::Array{Int64, 2}
  totalE::Matrix{Float64}
  valid::Bool
end

"""
    create_cache(sppl::SpeciesList, ml::GridLandscape)

Function to create a `Cache` for a SpeciesList and GridLandscape.
"""
function create_cache(sppl::SpeciesList, ml::GridLandscape)
  nm = zeros(Int64, size(ml.matrix))
  totalE = zeros(Float64, (size(ml.matrix, 2), numrequirements(typeof(sppl.species.requirement))))
  return Cache(nm, totalE, false)
end

"""
    create_cache(sppl::SpeciesList{SpeciesTypes{TR, R, MO, T}},
        ml::GridLandscape) where {TR, R, MO <: BirthOnlyMovement, T}

Function to create a `PlantCache` for a SpeciesList with movement only via
seed production (`BirthOnlyMovement`) and GridLandscape.
"""
function create_cache(sppl::SpeciesList{SpeciesTypes{TR, R, MO, T}},
    ml::GridLandscape) where {TR, R, MO <: BirthOnlyMovement, T}
  nm = zeros(Int64, size(ml.matrix))
  sb = zeros(Int64, size(ml.matrix))
  totalE = zeros(Float64, (size(ml.matrix, 2), numrequirements(typeof(sppl.species.requirement))))
  return PlantCache(nm, sb, totalE, false)
end

"""
    create_cache(sppl::SpeciesList, ml::EpiLandscape)

Function to create an `EpiCache` for a SpeciesList and EpiLandscape.
"""
function create_cache(sppl::SpeciesList, ml::EpiLandscape)
  vm = zeros(Float64, size(ml.matrix))
  return EpiCache(vm, false)
end


"""
    tematch(sppl::SpeciesList, abenv::AbstractAbiotic)

Function to check that the types of a trait list and habitat list are
the same for a species list (`sppl`) and abiotic environment (`abenv`).
"""
function tematch(sppl::SpeciesList, abenv::AbstractAbiotic)
    (eltype(sppl.species.traits) == eltype(abenv.habitat)) &&
    (iscontinuous(sppl.species.traits) == iscontinuous(abenv.habitat))
end

"""
    trmatch(sppl::SpeciesList, traitrel::AbstractTraitRelationship)

Function to check that the types of a trait list and trait relationship list are
the same for a species list (`sppl`) and trait relationship (`traitrel`).
"""
function trmatch(sppl::SpeciesList, traitrel::AbstractTraitRelationship)
    eltype(sppl.species.traits) == eltype(traitrel) &&
    (iscontinuous(sppl.species.traits) == iscontinuous(traitrel))
end

"""
    AbstractEcosystem{L <: AbstractLandscape, Part <: AbstractPartition, SL <: AbstractSpeciesList,
    TR <: AbstractTraitRelationship, LU <: AbstractLookup, C <: AbstractCache} <:
    AbstractMetacommunity{Float64, Matrix{Int64}, Matrix{Float64}, SL, Part}

Abstract type for `Ecosystem`s which is a subtype of `AbstractMetacommunity`.
"""
abstract type
    AbstractEcosystem{L <: AbstractLandscape, Part <: AbstractPartition, SL <: AbstractSpeciesList,
        TR <: AbstractTraitRelationship, LU <: AbstractLookup, C <: AbstractCache} <: AbstractMetacommunity{Float64, Matrix{Int64},
                                        Matrix{Float64}, SL, Part}
end

"""
    Ecosystem{L <: AbstractLandscape, Part <: AbstractPartition, SL <: AbstractSpeciesList,
    TR <: AbstractTraitRelationship, LU <: AbstractLookup, C <: AbstractCache} <:
    AbstractEcosystem{L, Part, SL, TR, LU, C}

`Ecosystem` type which is a subtype of `AbstractEcosystem`. It houses information
on the abundances of species in the landscape, `abundances`, information about
those species, `spplist`, information about the abiotic environment, `abenv`,
the `ordinariness` (for Diversity.jl calculations), the `relationship` between
the species and their environment, a pre-calculated `lookup` of all possible moves that could
be made by each species, a `cache`, and a list of `transitions`.
"""
mutable struct Ecosystem{L <: AbstractLandscape, Part <: AbstractPartition, SL <: AbstractSpeciesList,
    TR <: AbstractTraitRelationship, LU <: AbstractLookup, C <: AbstractCache} <: AbstractEcosystem{L, Part, SL, TR, LU, C}
  abundances::L
  spplist::SL
  abenv::Part
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::TR
  lookup::LU
  cache::C
  transitions::Union{Nothing, TransitionList}
end

"""
   Ecosystem(popfun::F, spplist::SpeciesList{T, Req}, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship; transitions::Union{Nothing, TransitionList} = nothing)
   where {F<:Function, T, Req}

Function to create an `Ecosystem` with species from a `SpeciesList`, environment from
 a `GridAbioticEnvironment`, an `AbstractTraitRelationship` between the two, and
 a list of `transitions`. The `Ecosystem` can be populated with a `popfun.`
"""
function Ecosystem(popfun::F, spplist::SpeciesList{T, Req}, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship; transitions::Union{Nothing, TransitionList} = nothing) where {F<:Function, T, Req}

    # Create matrix landscape of zero abundances
  ml = emptygridlandscape(abenv, spplist)
  # Populate this matrix with species abundances
  popfun(ml, spplist, abenv, rel)
  # Create lookup table of all moves and their probabilities
  lookup = SpeciesLookup(collect(map(k -> genlookups(abenv.habitat, k), getkernels(spplist.species.movement))))
  cache = create_cache(spplist, ml)
  return Ecosystem{typeof(ml), typeof(abenv), typeof(spplist), typeof(rel), typeof(lookup), typeof(cache)}(ml, spplist, abenv,
  missing, rel, lookup, cache, transitions)
end

"""
   Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv, rel::AbstractTraitRelationship;
   transitions::Union{Nothing, TransitionList} = nothing)

Function to create an `Ecosystem` with species from a `SpeciesList`, environment from
 a `GridAbioticEnvironment`, an `AbstractTraitRelationship` between the two, and
 a list of `transitions`. The `Ecosystem` is automatically and randomly populated
 from the start abundances in `spplist`.
"""
function Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship; transitions::Union{Nothing, TransitionList} = nothing)
   return Ecosystem(populate!, spplist, abenv, rel, transitions = transitions)
end

"""
   Ecosystem(abundances::EpiLandscape{U, VecRNGType}, epilist::EL, epienv::EE,
       ordinariness::Union{Matrix{Float64}, Missing}, relationship::ER, lookup::EpiLookup,
       vm::Array{Float64, 2}, initial_infected::Int64, valid::Bool, transitions::Union{Nothing, TransitionList}
       ) where {U <: Integer, VecRNGType <: AbstractVector{<:Random.AbstractRNG},
       EE <: AbstractEpiEnv, EL <: SpeciesList, ER <: AbstractTraitRelationship}

Function to create an `Ecosystem` with epi categories from a `SpeciesList`, environment from
 a `AbstractEpiEnvironment`, an `AbstractTraitRelationship` between the two,
 an `EpiLookup` of pre-determined moves, information for the `EpiCache`
 (`initial_infected`, `valid` and `vm`) and a list of `transitions`.
 An `EpiCache` is automatically generated.
"""
function Ecosystem(abundances::EpiLandscape{U, VecRNGType}, epilist::EL, epienv::EE,
    ordinariness::Union{Matrix{Float64}, Missing}, relationship::ER, lookup::EpiLookup,
    vm::Array{Float64, 2}, initial_infected::Int64, valid::Bool, transitions::Union{Nothing, TransitionList}
    ) where {U <: Integer, VecRNGType <: AbstractVector{<:Random.AbstractRNG},
    EE <: AbstractEpiEnv, EL <: SpeciesList, ER <: AbstractTraitRelationship}
  total_pop = sum(abundances.matrix, dims = 1)[1, :]
  sorted_grid_ids = sortperm(total_pop, rev = true)
  sorted_grid_ids = sorted_grid_ids[total_pop[sorted_grid_ids] .> 0]
  cache = EpiCache(vm, initial_infected, sorted_grid_ids, valid)
  return Ecosystem(abundances, epilist, epienv, ordinariness, relationship, lookup, cache, transitions)
end

"""
   Ecosystem(popfun::F, epilist::SpeciesList, epienv::GridEpiEnv,
         rel::AbstractTraitRelationship, intnum::U; initial_infected = 0,
         rngtype::Type{R} = Random.MersenneTwister,
         transitions = nothing) where {F<:Function, U <: Integer, R <: Random.AbstractRNG}

Function to create an `Ecosystem` with epi categories from a `SpeciesList`, environment from
 a `GridEpiEnv`, an `AbstractTraitRelationship` between the two,
 an indication of what integer type the abundances should be stored in `intnum`,
 an optional number of `initial_infected` to be seeded, and optional type for rand calls,
`rngtype`, and an optional list of `transitions`. The `Ecosystem` can be populated with a `popfun.`
"""
function Ecosystem(popfun::F, epilist::SpeciesList, epienv::GridEpiEnv,
      rel::AbstractTraitRelationship, intnum::U; initial_infected = 0,
      rngtype::Type{R} = Random.MersenneTwister,
      transitions = nothing) where {F<:Function, U <: Integer, R <: Random.AbstractRNG}

  # Create matrix landscape of zero abundances
  ml = emptyepilandscape(epienv, epilist, intnum, rngtype)

  # Populate this matrix with species abundances
  popfun(ml, epilist, epienv, rel)
  initial_pop = sum(ml.matrix, dims = 1)
  # Create lookup table of all moves and their probabilities
  local_lookup = genlookups(epienv, epilist.species.movement.localmoves)
  region_lookup = genlookups(epienv, epilist.species.movement.regionmoves, initial_pop)
  lookup = EpiLookup(local_lookup, region_lookup)
  vm = zeros(Float64, size(ml.matrix))
  return Ecosystem(ml, epilist, epienv, missing, rel, lookup, vm, initial_infected, false, transitions)
end

"""
   Ecosystem(epilist::SpeciesList, epienv::GridEpiEnv, rel::AbstractTraitRelationship,
           intnum::U = Int64(1); initial_infected = 0, rngtype::Type{R} = Random.MersenneTwister,
           transitions = nothing
           ) where {U <: Integer, R <: Random.AbstractRNG}

Function to create an `Ecosystem` with epi categories from a `SpeciesList`, environment from
 a `GridEpiEnv`, an `AbstractTraitRelationship` between the two,
 an indication of what integer type the abundances should be stored in `intnum`,
 an optional number of `initial_infected` to be seeded, and optional type for rand calls,
`rngtype`, and an optional list of `transitions`. The `Ecosystem` is automatically and randomly populated
from the start abundances in `spplist`.
"""
function Ecosystem(epilist::SpeciesList, epienv::GridEpiEnv, rel::AbstractTraitRelationship,
        intnum::U = Int64(1); initial_infected = 0, rngtype::Type{R} = Random.MersenneTwister,
        transitions = nothing
        ) where {U <: Integer, R <: Random.AbstractRNG}
    return Ecosystem(populate!, epilist, epienv, rel, intnum, initial_infected = initial_infected, rngtype = rngtype, transitions = transitions)
end

"""
    isapprox(epi_1::AbstractEpiSystem, epi_2::AbstractEpiSystem; kwargs...)

Compare two `Ecosystem`s for approximate equality. Specifically, compares the
`Landscape`s of the two systems.

## Keyword arguments
- Anything to pass to `Base.isapprox`.

!!! note
    You may want to pass in `atol` or `rtol` to loosen the equality tolerance.
"""
function Base.isapprox(eco_1::AbstractEcosystem, eco_2::AbstractEcosystem; kwargs...)
    return isapprox(eco_1.abundances, eco_2.abundances; kwargs...)
end

save(path::String, system::Ecosystem) = JLSO.save(path, :ecosystem => system)
load(path::String, obj_type::Type{Ecosystem}) = JLSO.load(path)[:ecosystem]


import Diversity.API: _getabundance
function _getabundance(eco::AbstractEcosystem, input::Bool)
    if input
        return eco.abundances.matrix
    else
        return _calcabundance(_gettypes(eco), eco.abundances.matrix / sum(eco.abundances.matrix))[1]
    end
end

import Diversity.API: _getmetaabundance
function _getmetaabundance(eco::AbstractEcosystem)
  return sumoversubcommunities(eco, _getabundance(eco))
end

import Diversity.API: _getpartition
function _getpartition(eco::AbstractEcosystem)
  return eco.abenv
end

import Diversity.API: _gettypes
function _gettypes(eco::AbstractEcosystem)
    return eco.spplist
end

import Diversity.API: _getordinariness!
function _getordinariness!(eco::AbstractEcosystem)
    if ismissing(eco.ordinariness)
        relab = getabundance(eco, false)
        eco.ordinariness = _calcordinariness(eco.spplist, relab)
    end
    return eco.ordinariness
end

import Diversity.API._getscale
function _getscale(eco::AbstractEcosystem)
    return _calcabundance(_gettypes(eco), getabundance(eco, false))[2]
end


"""
    CachedEcosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
        TR <: AbstractTraitRelationship} <: AbstractEcosystem{Part, SL, TR}

CachedEcosystem houses the same information as Ecosystem (see ?Ecosystem), but
holds the time period abundances as a CachedGridLandscape, so that they may
be present or missing.
"""
mutable struct CachedEcosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
    TR <: AbstractTraitRelationship} <: AbstractEcosystem{CachedGridLandscape, Part, SL, TR, SpeciesLookup, Cache}
  abundances::CachedGridLandscape
  spplist::SL
  abenv::Part
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::TR
  lookup::SpeciesLookup
  cache::Cache
end

"""
    CachedEcosystem(eco::Ecosystem, outputfile::String, rng::StepRangeLen)

Function to create a CachedEcosystem given an existing ecosystem, `eco`,
output folder to which the simulations are saved, `outputfile`, and a range of
times over which to simulate, `rng`.
"""
function CachedEcosystem(eco::Ecosystem, outputfile::String, rng::StepRangeLen)
    if size(eco.abenv.habitat, 3) > 1
        size(eco.abenv.habitat, 3) == length(rng) || error("Time range does not match habitat")
    end
    abundances = CachedGridLandscape(outputfile, rng)
    abundances.matrix[1] = eco.abundances
  CachedEcosystem{typeof(eco.abenv), typeof(eco.spplist), typeof(eco.relationship)}(abundances,
  eco.spplist, eco.abenv, eco.ordinariness, eco.relationship, eco.lookup, eco.cache)
end
