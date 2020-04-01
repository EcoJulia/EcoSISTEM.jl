"""
    AbstractEpiSystem

Abstract supertype for all disease system types.
"""
abstract type AbstractEpiSystem{Part <: AbstractEpiEnv, EL <: EpiList, TR <: AbstractTraitRelationship} <: AbstractMetacommunity{Float64, Matrix{Int64},
                                    Matrix{Float64}, EL, Part} end


mutable struct EpiCache
  netmigration::Array{Int64, 2}
  valid::Bool
end
"""
    EpiSystem{EE <: AbstractEpiEnv, EL <: EpiList, ER <: AbstractRelationship} <: AbstractEpiSystem{EE, EL, ER}

EpiSystem houses information on different disease classes, `epilist`, the environment, `epienv`, and their relationship to one another, `relationship`.
"""
mutable struct EpiSystem{EE <: AbstractEpiEnv, EL <: EpiList, ER <: AbstractTraitRelationship} <: AbstractEpiSystem{EE, EL, ER}
  abundances::EpiLandscape
  epilist::EL
  epienv::EE
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::ER
  lookup::Vector{Lookup}
  cache::EpiCache

  function EpiSystem{EE, EL, ER}(abundances::EpiLandscape,
    epilist::EL, epienv::EE, ordinariness::Union{Matrix{Float64}, Missing}, relationship::ER, lookup::Vector{Lookup}, cache::EpiCache) where {EE <:
     AbstractEpiEnv,
    EL <: EpiList, ER <: AbstractTraitRelationship}
    new{EE, EL, ER}(abundances, epilist, epienv, ordinariness, relationship, lookup, cache)
  end
end

function EpiSystem(popfun::Function, epilist::EpiList, epienv::GridEpiEnv, rel::AbstractTraitRelationship)

  # Create matrix landscape of zero abundances
  ml = emptyepilandscape(epienv, epilist)
  # Populate this matrix with species abundances
  popfun(ml, epilist, epienv, rel)
  # Create lookup table of all moves and their probabilities
  lookup_tab = collect(map(k -> genlookups(epienv.habitat, k), getkernels(epilist.movement)))
  nm = zeros(Int64, size(ml.matrix))
  EpiSystem{typeof(epienv), typeof(epilist), typeof(rel)}(ml, epilist, epienv, missing, rel, lookup_tab, EpiCache(nm, false))
end

function EpiSystem(epilist::EpiList, epienv::GridEpiEnv,
   rel::AbstractTraitRelationship)
   return EpiSystem(populate!, epilist, epienv, rel)
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

function invalidatecaches!(epi::EpiSystem)
    epi.ordinariness = missing
    epi.cache.netmigration .= 0
    epi.cache.valid = false
end

function getdispersaldist(epi::AbstractEpiSystem, sp::Int64)
  dist = epi.epilist.movement.kernels[sp].dist
  return dist
end
function getdispersaldist(epi::AbstractEpiSystem, sp::String)
  num = Compat.findall(epi.epilist.names.==sp)[1]
  getdispersaldist(epi, num)
end

function getdispersalvar(epi::AbstractEpiSystem, sp::Int64)
    var = (epi.epilist.movement.kernels[sp].dist)^2 * pi / 4
    return var
end
function getdispersalvar(epi::AbstractEpiSystem, sp::String)
    num = Compat.findall(epi.epilist.names.==sp)[1]
    getdispersalvar(epi, num)
end

function getlookup(epi::AbstractEpiSystem, sp::Int64)
    return epi.lookup[sp]
end
function getlookup(epi::AbstractEpiSystem, sp::String)
    num = Compat.findall(epi.epilist.names.==sp)[1]
    getlookup(epi, num)
end
