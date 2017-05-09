"""
    AbstractTraits{T}

Abstract supertype for all trait types, parameterised by traits of any type `T`.
"""
abstract AbstractTraits{T}
"""
    BasicTrait{T} <: AbstractTraits{T}

Basic trait type that holds information on a single trait for each species, of any
type `T`.
"""
type BasicTrait{T} <: AbstractTraits{T}
  trait::Vector{T}
end
