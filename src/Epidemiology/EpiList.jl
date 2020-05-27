"""
    VirusTypes{TR <: AbstractTraits,
                 T <: AbstractTypes} <: AbstractTypes
VirusTypes holds information on the virus classes, such as the name of each class, their trait match to the environment, initial abundances and types.
"""
mutable struct VirusTypes{TR <: AbstractTraits,
                 T <: AbstractTypes} <: AbstractTypes
  names::Vector{String}
  traits::TR
  abun::Vector{Int64}
  types::T

  function VirusTypes{TR, T}(names:: Vector{String}, traits::TR, abun::Vector{Int64}, types::T) where {
                       TR <: AbstractTraits,
                       T <: AbstractTypes}
      new{TR, T}(names, traits, abun, types)
  end
  function VirusTypes{TR, T}(traits::TR, abun::Vector{Int64}, types::T) where {TR <: AbstractTraits, T <: AbstractTypes}
      names = map(x -> "$x", 1:length(abun))
      new{TR, T}(names, traits, abun, types)
  end
end

"""
    HumanTypes{MO <: AbstractMovement,
                 T <: AbstractTypes} <: AbstractTypes
HumanTypes holds information on the human disease classes, such as the name of each class, their initial abundances and types, as well as how they disperse virus across the landscape.
"""
mutable struct HumanTypes{MO <: AbstractMovement,
                 T <: AbstractTypes} <: AbstractTypes
  names::Vector{String}
  abun::Vector{Int64}
  types::T
  movement::MO

  function HumanTypes{MO, T}(names:: Vector{String}, abun::Vector{Int64}, types::T, movement::MO) where {
                       MO <: AbstractMovement,
                       T <: AbstractTypes}
      new{MO, T}(names, abun, types, movement)
  end
  function HumanTypes{MO, T}(abun::Vector{Int64}, types::T, movement::MO) where {
                       MO <: AbstractMovement,
                       T <: AbstractTypes}
      names = map(x -> "$x", 1:length(abun))
      new{MO, T}(names, abun, types, movement)
  end
end

"""
    EpiList{P <: AbstractParams} <: AbstractTypes
Epi list houses all disease and virus class specific information, as well as parameters for model runs.
"""
mutable struct EpiList{P <: AbstractParams} <: AbstractTypes
    virus::VirusTypes
    human::HumanTypes
    params::P

    function EpiList{P}(virus::VirusTypes, human::HumanTypes, params::P) where {P <: AbstractParams}
      new{P}(virus, human, params)
    end
end

import Diversity.API: _gettypenames
function _gettypenames(el::EpiList, input::Bool)
    return _gettypenames(el.human.types, input)
end

import Diversity.API: _counttypes
function _counttypes(el::EpiList, input::Bool)
    return _counttypes(el.human.types, input) + _counttypes(el.virus.types, input)
end
function _counttypes(hm::HumanTypes, input::Bool)
    return _counttypes(hm.types, input)
end
function _counttypes(vr::VirusTypes, input::Bool)
    return _counttypes(vr.types, input)
end

import Diversity.API: _calcsimilarity
function _calcsimilarity(el::EpiList, a::AbstractArray)
    return _calcsimilarity(el.human.types, a)
end
import Diversity.API: floattypes
function floattypes(::EpiList)
    return Set([Float64])
end
import Diversity.API: _calcordinariness
function _calcordinariness(el::EpiList, a::AbstractArray)
    _calcordinariness(el.human.types, a, one(eltype(a)))
end
import Diversity.API: _calcabundance
function _calcabundance(el::EpiList, a::AbstractArray)
  return _calcabundance(el.human.types, a)
end
import Diversity.API._getdiversityname
function _getdiversityname(el::EpiList)
    return _getdiversityname(el.human.types)
end

"""
    SIS(traits::TR, virus_abun::Vector{Int64}, human_abun::Vector{Int64}, movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits, MO <: AbstractMovement, P <: AbstractParams}

Function to create an `EpiList` for an SIS model - creating the correct number of classes and checking dimensions.
"""
function SIS(traits::TR, virus_abun::Vector{Int64}, human_abun::Vector{Int64}, movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    names = ["Susceptible", "Infected", "Dead"]
    new_names = [ifelse(i == 1, "$j", "$j$i") for i in 1:age_categories, j in names][1:end]
    ht = UniqueTypes(length(new_names))
    human = HumanTypes{typeof(movement), typeof(ht)}(new_names, human_abun, ht, movement)

    virus_names = [ifelse(i == 1, "Virus", "Virus$i") for i in 1:length(virus_abun)]
    vt = UniqueTypes(length(virus_names))
    virus = VirusTypes{typeof(traits), typeof(vt)}(virus_names, traits, virus_abun, vt)

    length(movement.kernels) == length(new_names) || throw(DimensionMismatch("Movement vector doesn't match number of disease classes"))
    length(traits.mean) == length(virus_names) || throw(DimensionMismatch("Trait vector doesn't match number of virus classes"))
    size(params.transition, 1) == length(new_names) || throw(DimensionMismatch("Transition matrix doesn't match number of disease classes"))
  EpiList{typeof(params)}(virus, human, params)
end

"""
    SIR(traits::TR, virus_abun::Vector{Int64}, human_abun::Vector{Int64}, movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits, MO <: AbstractMovement, P <: AbstractParams}

Function to create an `EpiList` for an SIR model - creating the correct number of classes and checking dimensions.
"""
function SIR(traits::TR, virus_abun::Vector{Int64}, human_abun::Vector{Int64}, movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    names = ["Susceptible", "Infected", "Recovered", "Dead"]
    new_names = [ifelse(i == 1, "$j", "$j$i") for i in 1:age_categories, j in names][1:end]
    ht = UniqueTypes(length(new_names))
    human = HumanTypes{typeof(movement), typeof(ht)}(new_names, human_abun, ht, movement)

    virus_names = [ifelse(i == 1, "Virus", "Virus$i") for i in 1:length(virus_abun)]
    vt = UniqueTypes(length(virus_names))
    virus = VirusTypes{typeof(traits), typeof(vt)}(virus_names, traits, virus_abun, vt)

    length(movement.kernels) == length(new_names) || throw(DimensionMismatch("Movement vector doesn't match number of disease classes"))
    length(traits.mean) == length(virus_names) || throw(DimensionMismatch("Trait vector doesn't match number of virus classes"))
    size(params.transition, 1) == length(new_names) || throw(DimensionMismatch("Transition matrix doesn't match number of disease classes"))
  EpiList{typeof(params)}(virus, human, params)
end

"""
    SEIR(traits::TR, virus_abun::Vector{Int64}, human_abun::Vector{Int64}, movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits, MO <: AbstractMovement, P <: AbstractParams}

Function to create an `EpiList` for an SEIR model - creating the correct number of classes and checking dimensions.
"""
function SEIR(traits::TR, virus_abun::Vector{Int64}, human_abun::Vector{Int64}, movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    names = ["Susceptible", "Exposed", "Infected", "Recovered", "Dead"]
    new_names = [ifelse(i == 1, "$j", "$j$i") for i in 1:age_categories, j in names][1:end]
    ht = UniqueTypes(length(new_names))
    human = HumanTypes{typeof(movement), typeof(ht)}(new_names, human_abun, ht, movement)

    virus_names = [ifelse(i == 1, "Virus", "Virus$i") for i in 1:length(virus_abun)]
    vt = UniqueTypes(length(virus_names))
    virus = VirusTypes{typeof(traits), typeof(vt)}(virus_names, traits, virus_abun, vt)

    length(movement.kernels) == length(new_names) || throw(DimensionMismatch("Movement vector doesn't match number of disease classes"))
    length(traits.mean) == length(virus_names) || throw(DimensionMismatch("Trait vector doesn't match number of virus classes"))
    size(params.transition, 1) == length(new_names) || throw(DimensionMismatch("Transition matrix doesn't match number of disease classes"))
  EpiList{typeof(params)}(virus, human, params)
end

"""
    SEIRS(traits::TR, virus_abun::Vector{Int64}, human_abun::Vector{Int64}, movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits, MO <: AbstractMovement, P <: AbstractParams}

Function to create an `EpiList` for an SEIRS model - creating the correct number of classes and checking dimensions.
"""
function SEIRS(traits::TR, virus_abun::Vector{Int64}, human_abun::Vector{Int64}, movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    return SEIR(traits, virus_abun, human_abun, movement, params, age_categories)
end

"""
    SEI2HRD(traits::TR, virus_abun::Vector{Int64}, human_abun::Vector{Int64}, movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits, MO <: AbstractMovement, P <: AbstractParams}

Function to create an `EpiList` for an SEI2HRD model - creating the correct number of classes and checking dimensions.
"""
function SEI2HRD(traits::TR, virus_abun::Vector{Int64}, human_abun::Vector{Int64}, movement::MO, params::P, age_categories::Int64 = 1) where {TR <: AbstractTraits,
        MO <: AbstractMovement, P <: AbstractParams}

    names = ["Susceptible", "Exposed", "AsymptomaticInfected", "SymptomaticInfected", "Hospitalised", "Recovered", "Dead"]
    new_names = [ifelse(i == 1, "$j", "$j$i") for i in 1:age_categories, j in names][1:end]
    ht = UniqueTypes(length(new_names))
    human = HumanTypes{typeof(movement), typeof(ht)}(new_names, human_abun, ht, movement)

    virus_names = [ifelse(i == 1, "Virus", "Virus$i") for i in 1:length(virus_abun)]
    vt = UniqueTypes(length(virus_names))
    virus = VirusTypes{typeof(traits), typeof(vt)}(virus_names, traits, virus_abun, vt)

    length(movement.kernels) == length(new_names) || throw(DimensionMismatch("Movement vector doesn't match number of disease classes"))
    length(traits.mean) == length(virus_names) || throw(DimensionMismatch("Trait vector doesn't match number of virus classes"))
    size(params.transition, 1) == length(new_names) || throw(DimensionMismatch("Transition matrix doesn't match number of disease classes"))
    length(params.births) == length(new_names) || throw(DimensionMismatch("Birth vector doesn't match number of disease classes"))
  EpiList{typeof(params)}(virus, human, params)
end
