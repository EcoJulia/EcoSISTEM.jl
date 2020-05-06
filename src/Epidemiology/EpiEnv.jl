using Diversity
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref

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

This epi environment type holds a habitat and control strategy, as well as a string of subcommunity names.
"""
mutable struct GridEpiEnv{H, C} <: AbstractEpiEnv{H, C}
  habitat::H
  active::Array{Bool, 2}
  control::C
  names::Vector{String}
  function (::Type{GridEpiEnv{H, C}})(habitat::H, active::Array{Bool,2}, control::C, names::Vector{String} =
       map(x -> "$x", 1:countsubcommunities(habitat))) where {H, C}
    countsubcommunities(habitat) == length(names) ||
      error("Number of subcommunities must match subcommunity names")
    return new{H, C}(habitat, active, control, names)
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
    simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
    dimension::Tuple{Int64, Int64}, area::Unitful.Area{Float64},
    active::Array{Bool, 2}, control::C)

Function to create a simple `ContinuousHab` type epi
environment. It creates a `ContinuousHab` filled with a given value, `val`, dimensions (`dimension`) and a specified area (`area`). The rate of environmental change is specified using the parameter `rate`. If a Bool matrix of active grid squares is included, `active`, this is used, else one is created with all grid cells active.
"""
function simplehabitatAE(
    val::Union{Float64, Unitful.Quantity{Float64}},
    dimension::Tuple{Int64, Int64},
    area::Unitful.Area{Float64},
    active::Array{Bool, 2},
    control::C
) where C <: AbstractControl
    if typeof(val) <: Unitful.Temperature
        val = uconvert(K, val)
    end
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = simplehabitat(val, gridsquaresize, dimension)
    return GridEpiEnv{typeof(hab), typeof(control)}(hab, active, control)
end

function simplehabitatAE(
    val::Union{Float64, Unitful.Quantity{Float64}},
    dimension::Tuple{Int64, Int64},
    area::Unitful.Area{Float64},
    control::C
) where C <: AbstractControl

    active = fill(true, dimension)
    simplehabitatAE(val, dimension, area, active, control)
end
