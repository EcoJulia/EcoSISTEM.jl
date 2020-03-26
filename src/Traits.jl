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
function eltype(trait::DiscreteTrait{D}) where D
    return D
end


"""
    DiscreteEvolve(numTraits::Int64, tree::BinaryTree)

Function to evolve a discrete switching trait along a BinaryTree, `tree`. Takes in a number of traits, `numTraits` to be switched between and rate to switch between traits, `switch_rate` with default value of 0.5.
"""
function DiscreteEvolve(numTraits::Int64, tree::BinaryTree, switch_rate= 0.5)
    # Create traits and assign to tips
    trts = DataFrame(trait1 = collect(1:numTraits))
    assign_traits!(tree, 0.5, trts)
    # Get traits from tree
    return DiscreteTrait(Array(get_traits(tree, true)[:, 1]))
end

"""
    ContinuousEvolve(val::Union{Float64, Unitful.Quantity{Float64}}, var::Union{Float64, Unitful.Quantity{Float64}}, tree::BinaryTree)

Function to evolve a continuous trait along a BinaryTree, `tree` via Brownian motion. Takes in a starting value, `val` and a variance, `var`.
"""
function ContinuousEvolve(val::Union{Float64, Unitful.Quantity{Float64}},
    var::Union{Float64, Unitful.Quantity{Float64}}, tree::BinaryTree)
    # Create traits and assign to tips
    numspecies = length(getleafnames(tree))
    trts = DataFrame(start = ustrip(val), σ² = ustrip(var))
    assign_traits!(tree, trts)
    # Get traits from tree
    newtrts = get_traits(tree, true)
    newtrts[:start] = newtrts[:start] .* unit(val)
    return GaussTrait(newtrts[:start], newtrts[:σ²])
end

"""
    ContinuousTrait{C <: Number} <: AbstractTraits{T}

Abstract trait type that holds information on a single continuous trait for each species, of any Number type `C`.
"""
abstract type ContinuousTrait{C <: Number} <: AbstractTraits{C}
end

"""
    GaussTrait{C <: Number} <: ContinuousTrait{C}

Trait type that holds Gaussian mean and variance trait information for each species, of any number type `C`.
"""
mutable struct GaussTrait{C <: Number} <: ContinuousTrait{C}
  mean::Array{C, 1}
  var::Array{C, 1}
end
iscontinuous(trait::GaussTrait{C}) where C = true
function eltype(trait::GaussTrait{C}) where C
    return C
end
function GaussTrait(mean::Array{C, 1}, var::Array{C, 1}) where C  <: Unitful.Temperature
    meanK = uconvert.(K, mean)
    varK = ustrip.(var) .* K
    return GaussTrait{typeof(1.0K)}(meanK, varK)
end

"""
    TempBin{C <: Int} <: ContinuousTrait{C}

Trait type that holds binned temperature preference information created through ClimatePref. Holds an array of counts per temperature band (°C).
"""
mutable struct TempBin{C <: Int}<: ContinuousTrait{C}
  dist::Array{C, 2}
end
iscontinuous(trait::TempBin{C}) where C = true
function eltype(trait::TempBin{C}) where C
    return typeof(1.0°C)
end

"""
    RainBin{C <: Int} <: ContinuousTrait{C}

Trait type that holds binned rainfall preference information created through ClimatePref. Holds an array of counts per rainfall band (mm).
"""
mutable struct RainBin{C <: Int} <: ContinuousTrait{C}
  dist::Array{C, 2}
end
iscontinuous(trait::RainBin{C}) where C = true
function eltype(trait::RainBin{C}) where C
    return typeof(1.0mm)
end

"""
    TraitCollection2{T1, T2} <: AbstractTraits{Tuple{T1, T2}}

Trait collection that holds two trait types, `TR1` and `TR2`.
"""
mutable struct TraitCollection2{T1, T2} <: AbstractTraits{Tuple{T1, T2}}
    t1::T1
    t2::T2
end

iscontinuous(trait::TraitCollection2{T1, T2}) where {T1, T2} =
    [iscontinuous(trait.t1), iscontinuous(trait.t2)]
function eltype(trait::TraitCollection2)
    return [eltype(trait.t1), eltype(trait.t2)]
end

"""
    TraitCollection3{T1, T2, T3} <: AbstractTraits{Tuple{T1, T2, T3}}

Trait collection that holds three trait types, `TR1`, `TR2` and `TR3`.
"""
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
