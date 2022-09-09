using InteractiveUtils
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

function TransitionList(specialise=false)
    if specialise
        before = Vector{Union{rsubtypes(AbstractSetUp)...}}(undef,0)
        state_list = Vector{Union{rsubtypes(AbstractStateTransition)...}}(undef,0)
        place_list = Vector{Union{rsubtypes(AbstractPlaceTransition)...}}(undef,0)
        after = Vector{Union{rsubtypes(AbstractWindDown)...}}(undef,0)
    else
        before = Vector{AbstractSetUp}(undef,0)
        state_list = Vector{AbstractStateTransition}(undef,0)
        place_list = Vector{AbstractPlaceTransition}(undef,0)
        after = Vector{AbstractWindDown}(undef,0)
    end
    return TransitionList{eltype(before), eltype(state_list), eltype(place_list), eltype(after)}(before, state_list, place_list, after)
end

function specialise_transition_list(tl::TransitionList)
    setup = Vector{Union{unique(typeof.(tl.setup))...}}(tl.setup)
    state = Vector{Union{unique(typeof.(tl.state))...}}(tl.state)
    place = Vector{Union{unique(typeof.(tl.place))...}}(tl.place)
    winddown = Vector{Union{unique(typeof.(tl.winddown))...}}(tl.winddown)
    return TransitionList{eltype(setup), eltype(state), eltype(place), eltype(winddown)}(setup, state, place, winddown)
end

function rsubtypes(t)
    return rsubtypes([t], [])
end

function rsubtypes(types::AbstractVector, out)
    if isempty(types) 
        return out 
    else
        sub = subtypes.(types)
        return rsubtypes(filter(isabstracttype, vcat(sub...)), filter(!isabstracttype, vcat(out, types, sub...))) 
    end
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
