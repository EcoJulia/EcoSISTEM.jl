using Diversity
using DataFrames
using Unitful
using EcoSISTEM.Units
using Missings
using Compat
using RecipesBase


import Diversity: _calcabundance


abstract type AbstractCache end
"""
    Cache

Cache houses an integer array of moves made by all species in a timestep for the
update! function, `netmigration`.
"""
mutable struct Cache <: AbstractCache
  netmigration::Array{Int64, 2}
  totalE::Matrix{Float64}
  valid::Bool
end
mutable struct EpiCache <: AbstractCache
  virusmigration::Array{Float64, 2}
  initial_infected::Int64
  ordered_active::Vector{Int64}
  valid::Bool
end
mutable struct PlantCache <: AbstractCache
  netmigration::Array{Int64, 2}
  seedbank::Array{Int64, 2}
  totalE::Matrix{Float64}
  valid::Bool
end

function create_cache(sppl::SpeciesList{SpeciesTypes{TR, R, MO, T}},
    ml::GridLandscape) where {TR, R, MO <: AlwaysMovement, T}
  nm = zeros(Int64, size(ml.matrix))
  totalE = zeros(Float64, (size(ml.matrix, 2), numrequirements(typeof(sppl.species.requirement))))
  return Cache(nm, totalE, false)
end

function create_cache(sppl::SpeciesList{SpeciesTypes{TR, R, MO, T}},
    ml::GridLandscape) where {TR, R, MO <: BirthOnlyMovement, T}
  nm = zeros(Int64, size(ml.matrix))
  sb = zeros(Int64, size(ml.matrix))
  totalE = zeros(Float64, (size(ml.matrix, 2), numrequirements(typeof(sppl.species.requirement))))
  return PlantCache(nm, sb, totalE, false)
end

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


abstract type
    AbstractEcosystem{L <: AbstractLandscape, Part <: AbstractPartition, SL <: AbstractSpeciesList,
        TR <: AbstractTraitRelationship, LU <: AbstractLookup, C <: AbstractCache} <: AbstractMetacommunity{Float64, Matrix{Int64},
                                        Matrix{Float64}, SL, Part}
end


mutable struct Ecosystem{L <: AbstractLandscape, Part <: AbstractPartition, SL <: AbstractSpeciesList,
    TR <: AbstractTraitRelationship, LU <: AbstractLookup, C <: AbstractCache} <: AbstractEcosystem{L, Part, SL, TR, LU, C}
  abundances::L
  spplist::SL
  abenv::Part
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::TR
  lookup::LU
  cache::C
  transitions::Union{Missing, TransitionList}
end

function Ecosystem(popfun::F, spplist::SpeciesList{T, Req}, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship; transitions::Union{Missing, TransitionList} = missing) where {F<:Function, T, Req}

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

function Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship; transitions::Union{Missing, TransitionList} = missing)
   return Ecosystem(populate!, spplist, abenv, rel, transitions = transitions)
end

function Ecosystem(abundances::EpiLandscape{U, VecRNGType}, epilist::EL, epienv::EE,
    ordinariness::Union{Matrix{Float64}, Missing}, relationship::ER, lookup::EpiLookup,
    vm::Array{Float64, 2}, initial_infected::Int64, valid::Bool, transitions::Union{Missing, TransitionList}
    ) where {U <: Integer, VecRNGType <: AbstractVector{<:Random.AbstractRNG},
    EE <: AbstractEpiEnv, EL <: SpeciesList, ER <: AbstractTraitRelationship}
  total_pop = sum(abundances.matrix, dims = 1)[1, :]
  sorted_grid_ids = sortperm(total_pop, rev = true)
  sorted_grid_ids = sorted_grid_ids[total_pop[sorted_grid_ids] .> 0]
  cache = EpiCache(vm, initial_infected, sorted_grid_ids, valid)
  return Ecosystem(abundances, epilist, epienv, ordinariness, relationship, lookup, cache, transitions)
end

function Ecosystem(popfun::F, epilist::SpeciesList, epienv::GridEpiEnv,
      rel::AbstractTraitRelationship, intnum::U; initial_infected = 0,
      rngtype::Type{R} = Random.MersenneTwister,
      transitions = missing) where {F<:Function, U <: Integer, R <: Random.AbstractRNG}

  # Create matrix landscape of zero abundances
  ml = emptyepilandscape(epienv, epilist, intnum, rngtype)

  # Populate this matrix with species abundances
  popfun(ml, epilist, epienv, rel)
  initial_pop = sum(ml.matrix, dims = 1)
  # Create lookup table of all moves and their probabilities
  home_lookup = genlookups(epienv, epilist.species.movement.home)
  work_lookup = genlookups(epienv, epilist.species.movement.work, initial_pop)
  lookup = EpiLookup(home_lookup, work_lookup)
  vm = zeros(Float64, size(ml.matrix))
  return Ecosystem(ml, epilist, epienv, missing, rel, lookup, vm, initial_infected, false, transitions)
end

function Ecosystem(epilist::SpeciesList, epienv::GridEpiEnv, rel::AbstractTraitRelationship,
        intnum::U = Int64(1); initial_infected = 0, rngtype::Type{R} = Random.MersenneTwister,
        transitions = missing
        ) where {U <: Integer, R <: Random.AbstractRNG}
    return Ecosystem(populate!, epilist, epienv, rel, intnum, initial_infected = initial_infected, rngtype = rngtype, transitions = transitions)
end

function Ecosystem(epilist::SpeciesList, epienv::GridEpiEnv, rel::AbstractTraitRelationship,
        initial_population::A, intnum::U = Int64(1); initial_infected = 0,
        rngtype::Type{R} = Random.MersenneTwister,
        transitions = missing) where {U <: Integer, A <: AbstractArray, R <: Random.AbstractRNG}
    if size(initial_population) != size(epienv.active)
        msg = "size(initial_population)==$(size(initial_population)) != " *
            "size(epienv.active)==$(size(epienv.active))"
        throw(DimensionMismatch(msg))
    end
    epienv.active .&= .!_inactive.(initial_population)

    # Create matrix landscape of zero abundances
    ml = emptyepilandscape(epienv, epilist, intnum, rngtype)

    # Create lookup table of all moves and their probabilities
    home_lookup = genlookups(epienv, epilist.species.movement.home)
    work_lookup = genlookups(epienv, epilist.species.movement.work, initial_population[1:end])
    lookup = EpiLookup(home_lookup, work_lookup)

    vm = zeros(Float64, size(ml.matrix))

    epi = Ecosystem(ml, epilist, epienv, missing, rel, lookup, vm, initial_infected, false, transitions)

    # Add in the initial susceptible population
    # TODO Need to fix code so it doesn't rely on name of susceptible class
    idx = findfirst(occursin.("Susceptible", epilist.species.names))
    if idx == nothing
        msg = "epilist has no Susceptible category. epilist.names = $(epilist.species.names)"
        throw(ArgumentError(msg))
    end
    # Modify active cells based on new population
    initial_population = convert_population(initial_population, intnum)
    epi.abundances.grid[idx, :, :] .+= initial_population

    total_pop = sum(human(epi.abundances), dims = 1)[1, :]
    sorted_grid_ids = sortperm(total_pop, rev = true)
    sorted_grid_ids = sorted_grid_ids[total_pop[sorted_grid_ids] .> 0]
    epi.cache.ordered_active = sorted_grid_ids
    return epi
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
