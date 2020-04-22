
"""
    EpiList{TR <: AbstractTraits,
                 MO <: AbstractMovement,
                 T <: AbstractTypes,
                 P <: AbstractParams} <: AbstractTypes
Epi list houses all disease class specific information including trait information and movement types.
"""
mutable struct EpiList{TR <: AbstractTraits,
                 MO <: AbstractMovement,
                 T <: AbstractTypes,
                 P <: AbstractParams} <: AbstractTypes
  names::Vector{String}
  traits::TR
  abun::Vector{Int64}
  types::T
  movement::MO
  params::P

  function EpiList{TR, MO, T, P}(names:: Vector{String}, traits::TR, abun::Vector{Int64}, types::T, movement::MO, params::P) where {
                       TR <: AbstractTraits,
                       MO <: AbstractMovement,
                       T <: AbstractTypes,
                       P <: AbstractParams}
      new{TR, MO, T, P}(names, traits, abun, types,
       movement, params)
  end
  function EpiList{TR, MO, T, P}(traits::TR, abun::Vector{Int64}, types::T, movement::MO, params::P) where {
                       TR <: AbstractTraits,
                       MO <: AbstractMovement,
                       T <: AbstractTypes,
                       P <: AbstractParams}
      names = map(x -> "$x", 1:length(abun))
      new{TR, MO, T, P}(names, traits, abun, types,
       movement, params)
  end
end

import Diversity.API: _gettypenames
function _gettypenames(el::EpiList, input::Bool)
    return _gettypenames(el.types, input)
end
import Diversity.API: _counttypes
function _counttypes(el::EpiList, input::Bool)
    return _counttypes(el.types, input)
end
import Diversity.API: _calcsimilarity
function _calcsimilarity(el::EpiList, a::AbstractArray)
    return _calcsimilarity(el.types, a)
end
import Diversity.API: floattypes
function floattypes(::EpiList)
    return Set([Float64])
end
import Diversity.API: _calcordinariness
function _calcordinariness(el::EpiList, a::AbstractArray)
    _calcordinariness(el.types, a, one(eltype(a)))
end
import Diversity.API: _calcabundance
function _calcabundance(el::EpiList, a::AbstractArray)
  return _calcabundance(el.types, a)
end
import Diversity.API._getdiversityname
function _getdiversityname(el::EpiList)
    return _getdiversityname(el.types)
end


function SIR(traits::TR, abun::Vector{Int64},
    movement::MO, params::P) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    names = ["Virus", "Susceptible", "Infected", "Recovered"]
    types = UniqueTypes(length(names))
    length(abun) == length(names) || throw(DimensionMismatch("Abundance vector doesn't match number of disease classes"))
  EpiList{typeof(traits), typeof(movement), typeof(types), typeof(params)}(names, traits, abun, types, movement, params)
end

function SEI2HRD(traits::TR, abun::Vector{Int64},
    movement::MO, params::P) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    names = ["Virus", "Susceptible", "Exposed", "AsymptomaticInfected", "SymptomaticInfected", "Hospitalised", "Recovered", "Dead"]
    types = UniqueTypes(length(names))
    length(abun) == length(names) || throw(DimensionMismatch("Abundance vector doesn't match number of disease classes"))
    length(params.birth) == length(params.death) || throw(DimensionMismatch("Birth and death rates do not have the same number of classes"))
    length(params.birth) == length(names) || throw(DimensionMismatch("Birth and death rates do not have the correct number of classes"))
  EpiList{typeof(traits), typeof(movement), typeof(types), typeof(params)}(names, traits, abun, types, movement, params)
end
