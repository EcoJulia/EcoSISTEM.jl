# SPDX-License-Identifier: LGPL-3.0-or-later

import Base.eltype

"""
    AbstractTraits{T}

Abstract supertype for all trait types, parameterised by traits of any type `T`.
"""
abstract type AbstractTraits{T} end

eltype(::AbstractTraits{D}) where {D} = D

"""
    BasicTrait{T} <: AbstractTraits{T}

Basic trait type that holds information on a single trait for each species, of
any type `T`.
"""
struct DiscreteTrait{D} <: AbstractTraits{D}
    val::Vector{D}
end

iscontinuous(trait::DiscreteTrait) = false

struct LCtrait{D <: Number} <: AbstractTraits{D}
    vals::Vector{Vector{D}}
end

iscontinuous(trait::LCtrait) = false

function LCtrait(vals::Vector{Vector{<:AbstractFloat}})
    return LCtrait{typeof(1.0)}(vals)
end

"""
    DiscreteEvolve(numTraits::Int64, tree::BinaryTree)

Evolve a discrete switching trait along a BinaryTree, `tree`. Takes in a number
of traits, `numTraits` to be switched between and rate to switch between traits,
`switch_rate` with default value of 0.5.
"""
function DiscreteEvolve(numTraits::Int64, tree::BinaryTree, switch_rate = 0.5)
    # Create traits and assign to tips
    trts = DataFrame(trait1 = 1:numTraits)
    assign_traits!(tree, 0.5, trts)
    # Get traits from tree
    return DiscreteTrait(Array(get_traits(tree, true)[:, 1]))
end

"""
    ContinuousEvolve(val::Union{Float64, Unitful.Quantity{Float64}}, var::Union{Float64, Unitful.Quantity{Float64}}, tree::BinaryTree)

Evolve a continuous trait along a BinaryTree, `tree` via Brownian motion. Takes
in a starting value, `val` and a variance, `var`.
"""
function ContinuousEvolve(val::Union{Float64, Unitful.Quantity{Float64}},
                          var::Union{Float64, Unitful.Quantity{Float64}},
                          tree::BinaryTree)
    # Create traits and assign to tips
    numspecies = length(getleafnames(tree))
    trts = DataFrame(start = ustrip(val), σ² = ustrip(var))
    assign_traits!(tree, trts)
    # Get traits from tree
    newtrts = get_traits(tree, true)
    newtrts[!, :start] = newtrts[!, :start] .* unit(val)
    return GaussTrait(newtrts[!, :start], newtrts[!, :σ²])
end

"""
    ContinuousTrait{C <: Number} <: AbstractTraits{T}

Abstract trait type that holds information on a single continuous trait for each
species, of any Number type `C`.
"""
abstract type ContinuousTrait{C <: Number} <: AbstractTraits{C} end

"""
    GaussTrait{C <: Number} <: ContinuousTrait{C}

Trait type that holds a Gaussian habitat preference for each species, of any number
type `C`: the optimum `mean` and the standard deviation `sd` (the width of the
preference curve, which the [`Gauss`](@ref) relationship squares to a variance).
"""
struct GaussTrait{A <: NicheAxis, C <: Number} <: ContinuousTrait{C}
    mean::Vector{C}
    sd::Vector{C}
end

iscontinuous(trait::GaussTrait) = true

function GaussTrait(::Type{A}, mean::Vector{C},
                    sd::Vector{C}) where {A <: NicheAxis,
                                          C <: Unitful.Temperature}
    # `mean` is an absolute temperature (affine conversion), but `sd` is a temperature
    # interval (a width), so convert it via the scale only: subtracting the affine
    # unit's zero point (`0 * unit`, i.e. `0.0°F`/`0.0°C`/`0.0K`) turns the value into
    # its underlying interval before converting to K, so a σ given in °C/°F/K all land
    # on the correct K width (e.g. 9°F → 5K, not 9K or an offset).
    meanK = uconvert.(K, mean)
    sdK = uconvert.(K, sd .- 0 * unit(eltype(sd)))
    return GaussTrait{A, typeof(1.0K)}(meanK, sdK)
end
function GaussTrait(::Type{A}, mean::Vector{C},
                    sd::Vector{C}) where {A <: NicheAxis, C}
    return GaussTrait{A, C}(mean, sd)
end
# Back-compat: no axis given → `Unclassified` (matches positionally, not by axis).
function GaussTrait(mean::Vector{C}, sd::Vector{C}) where {C}
    return GaussTrait(Unclassified, mean, sd)
end

"""
    TempBin{C <: Int} <: ContinuousTrait{C}

Trait type that holds binned temperature preference information created through
ClimatePref. Holds an array of counts per temperature band (°C).
"""
struct TempBin{C <: Int} <: ContinuousTrait{C}
    dist::Matrix{C}
end

iscontinuous(trait::TempBin) = true

eltype(::TempBin) = typeof(1.0K)

"""
    RainBin{C <: Int} <: ContinuousTrait{C}

Trait type that holds binned rainfall preference information created through
ClimatePref. Holds an array of counts per rainfall band (mm).
"""
struct RainBin{C <: Int} <: ContinuousTrait{C}
    dist::Matrix{C}
end

iscontinuous(trait::RainBin) = true

eltype(::RainBin) = typeof(1.0mm)

"""
    TraitCollection2{T1, T2} <: AbstractTraits{Tuple{T1, T2}}

Trait collection that holds two trait types, `TR1` and `TR2`.
"""
struct TraitCollection2{T1, T2} <: AbstractTraits{Tuple{T1, T2}}
    t1::T1
    t2::T2
end

function iscontinuous(trait::TraitCollection2)
    return [iscontinuous(trait.t1), iscontinuous(trait.t2)]
end

eltype(trait::TraitCollection2) = [eltype(trait.t1), eltype(trait.t2)]

"""
    TraitCollection3{T1, T2, T3} <: AbstractTraits{Tuple{T1, T2, T3}}

Trait collection that holds three trait types, `TR1`, `TR2` and `TR3`.
"""
struct TraitCollection3{T1, T2, T3} <: AbstractTraits{Tuple{T1, T2, T3}}
    t1::T1
    t2::T2
    t3::T3
end

function iscontinuous(trait::TraitCollection3)
    return [
        iscontinuous(trait.t1),
        iscontinuous(trait.t2),
        iscontinuous(trait.t3)
    ]
end

function Base.eltype(trait::TraitCollection3)
    return [eltype(trait.t1), eltype(trait.t2), eltype(trait.t3)]
end
