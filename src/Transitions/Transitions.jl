
abstract type AbstractTransition end

abstract type AbstractStateTransition <: AbstractTransition end

abstract type AbstractPlaceTransition <: AbstractTransition end

abstract type AbstractSetUp <: AbstractTransition end

abstract type AbstractWindDown <: AbstractTransition end

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

function getdestination(rule::S) where S <: AbstractStateTransition
    return rule.destination
end

mutable struct TransitionList{T1 <: Union{Missing, AbstractSetUp}, T2 <: AbstractStateTransition,
    T3 <: AbstractPlaceTransition, T4 <: Union{Missing, AbstractWindDown}}
    setup::T1
    state::Array{T2, 1}
    place::Array{T3, 1}
    winddown::T4
end
