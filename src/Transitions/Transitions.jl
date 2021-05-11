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
    getprob(rule::AbstractStateTransition)

Generic function to get probability of a state transition happening.
"""
function getprob(rule::AbstractStateTransition)
    return rule.prob
end

"""
    getspecies(rule::AbstractStateTransition)

Generic function to get species to which a state transition happening.
"""
function getspecies(rule::AbstractStateTransition)
    return rule.species
end

"""
    getlocation(rule::AbstractStateTransition)

Generic function to get location where a state transition happening.
"""
function getlocation(rule::AbstractStateTransition)
    return rule.location
end

"""
    getspecies(rule::AbstractPlaceTransition)

Generic function to get species to which a place transition happening.
"""
function getspecies(rule::AbstractPlaceTransition)
    return rule.species
end

"""
    getlocation(rule::AbstractPlaceTransition)

Generic function to get location where a place transition happening.
"""
function getlocation(rule::AbstractPlaceTransition)
    return rule.location
end

"""
    getdestination(rule::AbstractStateTransition)

Generic function to get a destination to which a state transition happening.
"""
function getdestination(rule::AbstractStateTransition)
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
    addtransition!(tl::TransitionList, rule::AbstractStateTransition)

Add a state transition to the `TransitionList`.
"""
function addtransition!(tl::TransitionList, rule::AbstractStateTransition)
    push!(tl.state, rule)
end

"""
    addtransition!(tl::TransitionList, rule::AbstractPlaceTransition)

Add a place transition to the `TransitionList`.
"""
function addtransition!(tl::TransitionList, rule::AbstractPlaceTransition)
    push!(tl.place, rule)
end

"""
    addtransition!(tl::TransitionList, rule::AbstractSetUp)

Add a set up transition to the `TransitionList`.
"""
function addtransition!(tl::TransitionList, rule::AbstractSetUp)
    push!(tl.setup, rule)
end

"""
    addtransition!(tl::TransitionList, rule::AbstractWindDown)

Add a wind down transition to the `TransitionList`.
"""
function addtransition!(tl::TransitionList, rule::AbstractWindDown)
    push!(tl.winddown, rule)
end
