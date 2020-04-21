
abstract type ModelClass end

mutable struct SIR <: ModelClass
    dict::Dict
    function SIR()
        names = ["Virus", "Susceptible", "Infected", "Recovered"]
        dict = Dict(zip(names, 1:4))
        new(dict)
    end
end

mutable struct SEIR <: ModelClass
    dict::Dict
    function SEIR()
        names = ["Virus", "Exposed", "Susceptible", "Infected", "Recovered"]
        dict = Dict(zip(names, 1:5))
        new(dict)
    end
end

mutable struct SEI2HRD <: ModelClass
    dict::Dict
    function SEI2HRD()
        names = ["Virus", "Susceptible", "Exposed", "AsymptomaticInfected", "SymptomaticInfected", "Hospitalised", "Recovered", "Dead"]
        dict = Dict(zip(names, 1:8))
        new(dict)
    end
end

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
                 P <: AbstractParams,
                 MC <: ModelClass} <: AbstractTypes
  names::Vector{String}
  traits::TR
  abun::Vector{Int64}
  types::T
  movement::MO
  params::P
  model::MC

  function EpiList{TR, MO, T, P, MC}(names:: Vector{String}, traits::TR, abun::Vector{Int64}, types::T, movement::MO, params::P, model::MC) where {
                       TR <: AbstractTraits,
                       MO <: AbstractMovement,
                       T <: AbstractTypes,
                       P <: AbstractParams,
                       MC <: ModelClass}
      new{TR, MO, T, P, MC}(names, traits, abun, types,
       movement, params, model)
  end
  function EpiList{TR, MO, T, P, MC}(traits::TR, abun::Vector{Int64}, types::T, movement::MO, params::P, model::MC) where {
                       TR <: AbstractTraits,
                       MO <: AbstractMovement,
                       T <: AbstractTypes,
                       P <: AbstractParams,
                       MC <: ModelClass}
      names = map(x -> "$x", 1:length(abun))
      new{TR, MO, T, P, MC}(names, traits, abun, types,
       movement, params, model)
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

abstract type DiseaseClass end

mutable struct Virus <: DiseaseClass
end
mutable struct Susceptible <: DiseaseClass
end
mutable struct Infected <: DiseaseClass
end
mutable struct Recovered <: DiseaseClass
end

classDict = Dict("Virus" => Virus(), "Susceptible" => Susceptible(), "Infected" => Infected(), "Recovered" => Recovered())
function getclass(class::String)
    return classDict[class]
end
