
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

mutable struct TransitionList{T1 <: AbstractSetUp, T2 <: AbstractStateTransition,
    T3 <: AbstractPlaceTransition, T4 <: AbstractWindDown}
    setup::Array{T1, 1}
    state::Array{T2, 1}
    place::Array{T3, 1}
    winddown::Array{T4, 1}
end

function create_transition_list()
    before = Vector{AbstractSetUp}(undef,0)
    state_list = Vector{AbstractStateTransition}(undef,0)
    place_list = Vector{AbstractPlaceTransition}(undef,0)
    after = Vector{AbstractWindDown}(undef,0)
    return TransitionList(before, state_list, place_list, after)
end

function addtransition!(tl::TransitionList, rule::R) where R <: AbstractStateTransition
    push!(tl.state, rule)
end

function addtransition!(tl::TransitionList, rule::R) where R <: AbstractPlaceTransition
    push!(tl.place, rule)
end

function addtransition!(tl::TransitionList, rule::R) where R <: AbstractSetUp
    push!(tl.setup, rule)
end

function addtransition!(tl::TransitionList, rule::R) where R <: AbstractWindDown
    push!(tl.winddown, rule)
end
