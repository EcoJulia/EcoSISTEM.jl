using JLSO

"""
    AbstractEpiSystem

Abstract supertype for all disease system types.
"""
abstract type AbstractEpiSystem{Part <: AbstractEpiEnv, EL <: EpiList, TR <: AbstractTraitRelationship} <: AbstractMetacommunity{Float64, Matrix{Int64},
                                    Matrix{Float64}, EL, Part} end


mutable struct EpiCache
  virusmigration::Array{Float64, 2}
  valid::Bool
end



mutable struct EpiLookup
  x::Vector{Int64}
  y::Vector{Int64}
  p::Vector{Float64}
  pnew::Vector{Float64}
  moves::Vector{Float64}
end
EpiLookup(df::DataFrame) = EpiLookup(df[!, :X], df[!, :Y], df[!, :Prob],
zeros(Float64, nrow(df)), zeros(Int64, nrow(df)))

"""
    EpiSystem{EE <: AbstractEpiEnv, EL <: EpiList, ER <: AbstractRelationship} <: AbstractEpiSystem{EE, EL, ER}

EpiSystem houses information on different disease classes, `epilist`, the environment, `epienv`, and their relationship to one another, `relationship`.

See `help?>plot_epidynamics` and `help?>plot_epiheatmaps` for relevant plotting functions.
"""
mutable struct EpiSystem{U <: Integer, EE <: AbstractEpiEnv, EL <: EpiList, ER <: AbstractTraitRelationship} <: AbstractEpiSystem{EE, EL, ER}
  abundances::EpiLandscape{U}
  epilist::EL
  epienv::EE
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::ER
  lookup::Vector{EpiLookup}
  cache::EpiCache

  function EpiSystem{U, EE, EL, ER}(abundances::EpiLandscape{U},
    epilist::EL, epienv::EE, ordinariness::Union{Matrix{Float64}, Missing}, relationship::ER, lookup::Vector{EpiLookup}, cache::EpiCache) where {U <: Integer, EE <:
     AbstractEpiEnv,
    EL <: EpiList, ER <: AbstractTraitRelationship}
    new{U, EE, EL, ER}(abundances, epilist, epienv, ordinariness, relationship, lookup, cache)
  end
end

function EpiSystem(popfun::F, epilist::EpiList, epienv::GridEpiEnv,
    rel::AbstractTraitRelationship, intnum::U) where {F<:Function, U <: Integer}

  # Create matrix landscape of zero abundances
  ml = emptyepilandscape(epienv, epilist, intnum)
  # Populate this matrix with species abundances
  popfun(ml, epilist, epienv, rel)
  # Create lookup table of all moves and their probabilities
  lookup_tab = collect(map(k -> genlookups(epienv, k), getkernels(epilist.human.movement)))
  vm = zeros(Float64, size(ml.matrix))
  EpiSystem{U, typeof(epienv), typeof(epilist), typeof(rel)}(ml, epilist, epienv, missing, rel, lookup_tab, EpiCache(vm, false))
end

function EpiSystem(epilist::EpiList, epienv::GridEpiEnv, rel::AbstractTraitRelationship, intnum::U = Int64(1)) where U <: Integer
    return EpiSystem(populate!, epilist, epienv, rel, intnum)
end

function EpiSystem(epilist::EpiList, epienv::GridEpiEnv, rel::AbstractTraitRelationship, initial_population, intnum::U = Int64(1)) where U <: Integer
    if size(initial_population) != size(epienv.active)
        msg = "size(initial_population)==$(size(initial_population)) != " *
            "size(epienv.active)==$(size(epienv.active))"
        throw(DimensionMismatch(msg))
    end
    epi = EpiSystem(epilist, epienv, rel, intnum)
    # Add in the initial susceptible population
    idx = findfirst(epilist.human.names .== "Susceptible")
    if idx == nothing
        msg = "epilist has no Susceptible category. epilist.names = $(epilist.human.names)"
        throw(ArgumentError(msg))
    end
    # Modify active cells based on new population
    epi.epienv.active .&= .!_inactive.(initial_population)
    initial_population = convert_population(initial_population, intnum)
    epi.abundances.grid[idx, :, :] .+= initial_population
    return epi
end

"""
    isapprox(epi_1::AbstractEpiSystem, epi_2::AbstractEpiSystem; kwargs...)

Compare two `EpiSystem`s for approximate equality. Specifically, compares the
`EpiLandscape`s of the two systems.

## Keyword arguments
- Anything to pass to `Base.isapprox`.

!!! note
    You may want to pass in `atol` or `rtol` to loosen the equality tolerance.
"""
function Base.isapprox(epi_1::AbstractEpiSystem, epi_2::AbstractEpiSystem; kwargs...)
    return isapprox(epi_1.abundances, epi_2.abundances; kwargs...)
end

save(path::String, system::EpiSystem) = JLSO.save(path, :episystem => system)
load(path::String, obj_type::Type{EpiSystem}) = JLSO.load(path)[:episystem]


"""
    genlookups(hab::AbstractHabitat, mov::GaussianMovement)

Function to generate lookup tables, which hold information on the probability
of moving to neighbouring squares.
"""
function genlookups(epienv::AbstractEpiEnv, mov::GaussianKernel)
    hab = epienv.habitat
    sd = (2 * mov.dist) / sqrt(pi)
    relsize =  _getgridsize(hab) ./ sd
    m = maximum(_getdimension(hab))
    p = mov.thresh
    return EpiLookup(_lookup(relsize, m, p, _gaussian_disperse))
end
function genlookups(epienv::AbstractEpiEnv, mov::LongTailKernel)
    hab = epienv.habitat
    sd = (2 * mov.dist) / sqrt(pi)
    relsize =  _getgridsize(hab) ./ sd
    m = maximum(_getdimension(hab))
    p = mov.thresh
    b = mov.shape
    return EpiLookup(_lookup(relsize, m, p, b, _2Dt_disperse))
end

function getsize(epi::AbstractEpiSystem)
  return _getsize(epi.epienv.habitat)
end

function getgridsize(epi::AbstractEpiSystem)
  return _getgridsize(epi.epienv.habitat)
end

function getdimension(epi::AbstractEpiSystem)
    return _getdimension(epi.epienv.habitat)
end


function gettraitrel(epi::AbstractEpiSystem)
  return epi.relationship
end

function gethabitat(epi::AbstractEpiSystem)
  return epi.epienv.habitat
end

import Diversity.API: _getabundance
function _getabundance(epi::AbstractEpiSystem, input::Bool)
    if input
        return epi.abundances.matrix
    else
        return _calcabundance(_gettypes(epi), epi.abundances.matrix / sum(epi.abundances.matrix))[1]
    end
end

import Diversity.API: _getmetaabundance
function _getmetaabundance(epi::AbstractEpiSystem)
  return sumoversubcommunities(epi, _getabundance(epi))
end


import Diversity.API: _getpartition
function _getpartition(epi::AbstractEpiSystem)
  return epi.epienv
end
import Diversity.API: _gettypes
function _gettypes(epi::AbstractEpiSystem)
    return epi.epilist
end
import Diversity.API: _getordinariness!
function _getordinariness!(epi::AbstractEpiSystem)
    if ismissing(epi.ordinariness)
        relab = getabundance(epi, false)
        epi.ordinariness = _calcordinariness(epi.epilist, relab)
    end
    return epi.ordinariness
end

import Diversity.API._getscale
function _getscale(epi::AbstractEpiSystem)
    return _calcabundance(_gettypes(epi), getabundance(epi, false))[2]
end

function invalidatecaches!(epi::AbstractEpiSystem)
    epi.ordinariness = missing
    epi.cache.virusmigration .= 0
    epi.cache.valid = false
end

function getdispersaldist(epi::AbstractEpiSystem, sp::Int64)
  dist = epi.epilist.human.movement.kernels[sp].dist
  return dist
end
function getdispersaldist(epi::AbstractEpiSystem, sp::String)
  num = Compat.findall(epi.epilist.human.names.==sp)[1]
  getdispersaldist(epi, num)
end

function getdispersalvar(epi::AbstractEpiSystem, sp::Int64)
    var = (epi.epilist.human.movement.kernels[sp].dist)^2 * pi / 4
    return var
end
function getdispersalvar(epi::AbstractEpiSystem, sp::String)
    num = Compat.findall(epi.epilist.human.names.==sp)[1]
    getdispersalvar(epi, num)
end

function getlookup(epi::AbstractEpiSystem, sp::Int64)
    return epi.lookup[sp]
end
function getlookup(epi::AbstractEpiSystem, sp::String)
    num = Compat.findall(epi.epilist.human.names.==sp)[1]
    getlookup(epi, num)
end
