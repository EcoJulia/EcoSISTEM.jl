import Base.eltype
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
mutable struct DiscreteTrait{D} <: AbstractTraits{D}
  val::Array{D, 1}
end


iscontinuous(trait::DiscreteTrait) = false
function eltype{D}(trait::DiscreteTrait{D})
    return D
end
mutable struct ContinuousTrait{C <: Number} <: AbstractTraits{C}
  mean::Array{C, 1}
  var::Array{C, 1}
end
iscontinuous(trait::ContinuousTrait) = true
function eltype{C}(trait::ContinuousTrait{C})
    return C
end

mutable struct TraitCollection2{T1, T2} <: AbstractTraits{Tuple{T1, T2}}
    t1::T1
    t2::T2
end
iscontinuous(trait::TraitCollection2) = [iscontinuous(trait.t1), iscontinuous(trait.t2)]
function eltype(trait::TraitCollection2)
    return [eltype(trait.t1), eltype(trait.t2)]
end
mutable struct TraitCollection3{T1, T2, T3} <: AbstractTraits{Tuple{T1, T2, T3}}
    t1::T1
    t2::T2
    t3::T3
end
iscontinuous(trait::TraitCollection3) = [iscontinuous(trait.t1),
iscontinuous(trait.t2), iscontinuous(trait.t3)]
function eltype(trait::TraitCollection3)
    return [eltype(trait.t1), eltype(trait.t2), eltype(trait.t3)]
end
