"""
    AbstractTransition

Abstract type for all transitions.
"""
abstract type AbstractTransition end

"""
    AbstractStateTransition <: AbstractTransition

Abstract type for state transitions. State transitions occur between types within
a single grid square, e.g. growth between age categories,
or transitions between disease categories.
"""
abstract type AbstractStateTransition <: AbstractTransition end

"""
    AbstractPlaceTransition <: AbstractTransition

Abstract type for place transitions. Place transitions occur across multiple grid
squares for a single type, e.g. dispersal of seeds across the landscape.
"""
abstract type AbstractPlaceTransition <: AbstractTransition end

"""
    AbstractSetUp <: AbstractTransition

Abstract type for transitions that set up the `Ecosystem` at
     the start of the timestep.
"""
abstract type AbstractSetUp <: AbstractTransition end

"""
    AbstractWindDown <: AbstractTransition

Abstract type for transitions that wind down the `Ecosystem` at
     the end of the timestep.
"""
abstract type AbstractWindDown <: AbstractTransition end

"""
    getprob(rule::S) where S <: AbstractStateTransition

Generic function to get probability of a state transition happening.
"""
function getprob(rule::S) where S <: AbstractStateTransition
    return rule.prob
end

"""
    getspecies(rule::S) where S <: AbstractStateTransition

Generic function to get species to which a state transition happening.
"""
function getspecies(rule::S) where S <: AbstractStateTransition
    return rule.species
end

"""
    getlocation(rule::S) where S <: AbstractStateTransition

Generic function to get location where a state transition happening.
"""
function getlocation(rule::S) where S <: AbstractStateTransition
    return rule.location
end

"""
    getspecies(rule::P) where P <: AbstractPlaceTransition

Generic function to get species to which a place transition happening.
"""
function getspecies(rule::P) where P <: AbstractPlaceTransition
    return rule.species
end

"""
    getlocation(rule::P) where P <: AbstractPlaceTransition

Generic function to get location where a place transition happening.
"""
function getlocation(rule::P) where P <: AbstractPlaceTransition
    return rule.location
end

"""
    getdestination(rule::S) where S <: AbstractStateTransition

Generic function to get a destination to which a state transition happening.
"""
function getdestination(rule::S) where S <: AbstractStateTransition
    return rule.destination
end

"""
    TransitionList{T1 <: AbstractSetUp, T2 <: AbstractStateTransition,
    T3 <: AbstractPlaceTransition, T4 <: AbstractWindDown}

`TransitionList` type which houses information on `setup`, `state`,
`place` and `winddown` transitions.
"""
mutable struct TransitionList{T1 <: AbstractSetUp, T2 <: AbstractStateTransition,
    T3 <: AbstractPlaceTransition, T4 <: AbstractWindDown}
    setup::Vector{T1}
    state::Vector{T2}
    place::Vector{T3}
    winddown::Vector{T4}
end

"""
    create_transition_list()

Create an empty `TransitionList`.
"""
function create_transition_list()
    before = Vector{AbstractSetUp}(undef,0)
    state_list = Vector{AbstractStateTransition}(undef,0)
    place_list = Vector{AbstractPlaceTransition}(undef,0)
    after = Vector{AbstractWindDown}(undef,0)
    return TransitionList(before, state_list, place_list, after)
end

"""
    addtransition!(tl::TransitionList, rule::R) where R <: AbstractStateTransition

Add a state transition to the `TransitionList`.
"""
function addtransition!(tl::TransitionList, rule::R) where R <: AbstractStateTransition
    push!(tl.state, rule)
end

"""
    addtransition!(tl::TransitionList, rule::R) where R <: AbstractPlaceTransition

Add a place transition to the `TransitionList`.
"""
function addtransition!(tl::TransitionList, rule::R) where R <: AbstractPlaceTransition
    push!(tl.place, rule)
end

"""
    addtransition!(tl::TransitionList, rule::R) where R <: AbstractSetUp

Add a set up transition to the `TransitionList`.
"""
function addtransition!(tl::TransitionList, rule::R) where R <: AbstractSetUp
    push!(tl.setup, rule)
end

"""
    addtransition!(tl::TransitionList, rule::R) where R <: AbstractWindDown

Add a wind down transition to the `TransitionList`.
"""
function addtransition!(tl::TransitionList, rule::R) where R <: AbstractWindDown
    push!(tl.winddown, rule)
end
