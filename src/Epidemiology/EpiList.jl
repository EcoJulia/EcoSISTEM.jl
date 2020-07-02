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
  susceptible::Vector{Int64}
  infectious::Vector{Int64}

  function HumanTypes{MO, T}(names:: Vector{String}, abun::Vector{Int64}, types::T, movement::MO, susceptible::Vector{Int64}, infectious::Vector{Int64}) where {
                       MO <: AbstractMovement,
                       T <: AbstractTypes}
      new{MO, T}(names, abun, types, movement, susceptible, infectious)
  end
  function HumanTypes{MO, T}(abun::Vector{Int64}, types::T, movement::MO, susceptible::Vector{Int64}, infectious::Vector{Int64}) where {
                       MO <: AbstractMovement,
                       T <: AbstractTypes}
      names = map(x -> "$x", 1:length(abun))
      new{MO, T}(names, abun, types, movement, susceptible, infectious)
  end
end

"""
    EpiList{P <: AbstractParams} <: AbstractTypes
Epi list houses all disease and virus class specific information, as well as parameters for model runs.
"""
mutable struct EpiList{P <: AbstractParams, V <: VirusTypes, H <: HumanTypes} <: AbstractTypes
    virus::V
    human::H
    params::P

    function EpiList{P, V, H}(virus::V, human::H, param::P) where {P <: AbstractParams, V <: VirusTypes, H <: HumanTypes}
      new{P, V, H}(virus, human, param)
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
    EpiList(traits::TR, virus_abun::NamedTuple, human_abun::NamedTuple, disease_classes::NamedTuple, movement::MO, param::P, age_categories::Int64 = 1) where {TR <: AbstractTraits, MO <: AbstractMovement, P <: AbstractParams}

Function to create an `EpiList` for any type of epidemiological model - creating the correct number of classes and checking dimensions.
"""
function EpiList(traits::TR, virus_abun::NamedTuple, human_abun::NamedTuple,
                 disease_classes::NamedTuple, movement::MO, param::P,
                 age_categories::Int64 = 1) where {TR <: AbstractTraits, MO <: AbstractMovement, P <: AbstractParams}
    # Test for susceptibility/infectiousness categories
    haskey(disease_classes, :infectious) ||
        error("Missing 'infectious' key - vector of infectious categories")
    haskey(disease_classes, :susceptible) ||
        error("Missing 'susceptible' key - vector of susceptible categories")

    # Extract infectious/susceptible categories
    susceptible = disease_classes.susceptible
    infectious = disease_classes.infectious

    # Find their index locations in the names list
    names = collect(string.(keys(human_abun)))
    abuns = vcat(collect(human_abun)...)
#    rm_idx = indexin([susceptible; infectious], abuns)
#    deleteat!(abuns, rm_idx)

    new_names = [ifelse(i == 1, "$j", "$j$i") for i in 1:age_categories,
                                                   j in names][1:end]
    sus = findall(occursin.(susceptible, new_names))
    inf = vcat(map(i -> findall(occursin.(infectious[i], new_names)),
                   eachindex(infectious))...)
    ht = UniqueTypes(length(new_names))
    human = HumanTypes{typeof(movement), typeof(ht)}(new_names, Int64.(abuns),
                                                     ht, movement, sus, inf)

    virus_names = collect(string.(keys(virus_abun)))

    vt = UniqueTypes(length(virus_names))
    virus = VirusTypes{typeof(traits), typeof(vt)}(virus_names, traits, vcat(collect(virus_abun)...), vt)

    length(sus) == length(susceptible) * age_categories ||
        throw(DimensionMismatch("Number of susceptible categories is incorrect"))
    length(inf) == length(infectious) * age_categories ||
        throw(DimensionMismatch("Number of infectious categories is incorrect"))
    length(traits.mean) == length(virus_names) ||
        throw(DimensionMismatch("Trait vector doesn't match number of virus classes"))
    size(param.transition, 1) == length(new_names) ||
        throw(DimensionMismatch("Transition matrix doesn't match number of disease classes"))
    return EpiList{typeof(param), typeof(virus), typeof(human)}(virus, human, param)
end
