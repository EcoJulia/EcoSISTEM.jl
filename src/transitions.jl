
abstract type AbstractTransition end

abstract type AbstractStateTransition <: AbstractTransition end

abstract type AbstractPlaceTransition <: AbstractTransition end

function getprob(rule::S) where S <: AbstractStateTransition
    return rule.prob
end

function getspecies(rule::S) where S <: AbstractStateTransition
    return rule.species
end

function getlocation(rule::S) where S <: AbstractStateTransition
    return rule.location
end

function getspecies(rule::P) where P <: AbstractPlaceTransition
    return rule.species
end

function getlocation(rule::P) where P <: AbstractPlaceTransition
    return rule.location
end

mutable struct BirthProcess{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct DeathProcess{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct BirthDeathProcess{U <: Unitful.Units} <: AbstractStateTransition
    birth::BirthProcess{U}
    death::DeathProcess{U}
end

mutable struct AllDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

mutable struct TransitionList{T1 <: AbstractStateTransition, T2 <: AbstractPlaceTransition}
    state::Array{T1, 1}
    place::Array{T2, 1}
end

function create_transitions(spplist::SpeciesList, abenv::GridAbioticEnv)

    state_list = [BirthDeathProcess(BirthProcess(spp, loc, spplist.params.birth[spp]), DeathProcess(spp, loc, spplist.params.death[spp])) for spp in eachindex(spplist.names) for loc in eachindex(abenv.habitat.matrix)]

    place_list = [AllDisperse(spp, loc) for spp in eachindex(spplist.names) for loc in eachindex(abenv.habitat.matrix)]

    return TransitionList(state_list, place_list)
end
