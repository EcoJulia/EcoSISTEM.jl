"""
    AbstractTraits{T}

Abstract supertype for all trait types, parameterised by traits of any type `T`.
"""
abstract type AbstractTraits{T} end
"""
    BasicTrait{T} <: AbstractTraits{T}

Basic trait type that holds information on a single trait for each species, of any
type `T`.
"""
mutable struct BasicTrait{T} <: AbstractTraits{T}
  trait::Array{T, 1}
end

mutable struct TempTrait{T <: Unitful.Temperature{Float64}} <: AbstractTraits{T}
  mean::Array{T, 1}
  var::Array{T, 1}
end
