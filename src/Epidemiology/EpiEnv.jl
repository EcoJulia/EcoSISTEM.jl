using Diversity
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref

using Diversity.API

"""
    AbstractEpiEnv{H <: AbstractHabitat, C <: AbstractControl} <:
    AbstractPartition{H}

Abstract supertype for all epi environment types and a subtype of AbstractPartition.
"""
abstract type AbstractEpiEnv{H <: AbstractHabitat, C <: AbstractControl} <:
   AbstractPartition{H} end

"""
    GridEpiEnv{H, C} <: AbstractAbiotic{H, C}

This epi environment type holds a habitat and control strategy, as well as a string of
subcommunity names, and initial susceptible population.
"""
mutable struct GridEpiEnv{H, C, A} <: AbstractEpiEnv{H, C}
    habitat::H
    active::A
    control::C
    names::Vector{String}
    function (::Type{GridEpiEnv{H, C}})(
        habitat::H,
        active::A,
        control::C,
        names::Vector{String}=map(x -> "$x", 1:countsubcommunities(habitat))
    ) where {H, C, A <: AbstractMatrix{Bool}}
        countsubcommunities(habitat) == length(names) ||
            error("Number of subcommunities must match subcommunity names")
        if (size(habitat, 1), size(habitat, 2)) != size(active)
            throw(DimensionMismatch(
                "size(habitat)=$(size(habitat, 1)), $(size(habitat, 2)) != " *
                "size(active)=$(size(active))"
            ))
        end
        return new{H, C, A}(habitat, active, control, names)
    end
end

import Diversity.API: _countsubcommunities
function _countsubcommunities(epienv::GridEpiEnv)
  return countsubcommunities(epienv.habitat)
end
import Diversity.API: _getsubcommunitynames
function _getsubcommunitynames(epienv::GridEpiEnv)
    return epienv.names
end

"""
    function simplehabitatAE(
        val::Union{Float64, Unitful.Quantity{Float64}},
        dimension::Tuple{Int64, Int64},
        area::Unitful.Area{Float64},
        active::AbstractMatrix{Bool},
        control::C,
    )

Function to create a simple `ContinuousHab` type epi environment. It creates a
`ContinuousHab` filled with a given value `val`, of dimensions `dimension` and specified
area `area`. If a Bool matrix `active` of active grid squares is included, this is used,
else one is created with all grid cells active.

!!! note
    The simulation grid will be shrunk so that it tightly wraps the active values
"""
function simplehabitatAE(
    val::Union{Float64, Unitful.Quantity{Float64}},
    dimension::Tuple{Int64, Int64},
    area::Unitful.Area{Float64},
    active::M,
    control::C,
) where {C <: AbstractControl, M <: AbstractMatrix{Bool}}
    if typeof(val) <: Unitful.Temperature
        val = uconvert(K, val)
    end
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))

    # Shrink to active region
    # This doesn't change the gridsquaresize
    active = shrink_to_active(active, active)
    dimension = size(active)

    hab = simplehabitat(val, gridsquaresize, dimension)
    return GridEpiEnv{typeof(hab), typeof(control)}(hab, active, control)
end

function simplehabitatAE(
    val::Union{Float64, Unitful.Quantity{Float64}},
    dimension::Tuple{Int64, Int64},
    area::Unitful.Area{Float64},
    control::C,
) where C <: AbstractControl
    active = fill(true, dimension)
    return simplehabitatAE(val, dimension, area, active, control)
end
