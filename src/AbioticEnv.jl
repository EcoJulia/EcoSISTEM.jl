using Diversity
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using Compat
using ClimatePref

using Diversity.API
"""
    AbstractAbiotic{H <: AbstractHabitat, B <: AbstractBudget} <: AbstractPartition

Abstract supertype for all abiotic environment types and a subtype of
AbstractPartition
"""
abstract type AbstractAbiotic{H <: AbstractHabitat, B <: AbstractBudget} <:
   AbstractPartition end

"""
    GridAbioticEnv{H, B} <: AbstractAbiotic{H, B}

This abiotic environment type holds a habitat and budget, as well as a string
of subcommunity names.
"""
mutable struct GridAbioticEnv{H, B} <: AbstractAbiotic{H, B}
  habitat::H
  active::Array{Bool, 2}
  budget::B
  names::Vector{String}
  function (::Type{GridAbioticEnv{H, B}})(habitat::H, active::Array{Bool,2},
       budget::B, names::Vector{String} =
       map(x -> "$x", 1:countsubcommunities(habitat))) where {H, B}

    countsubcommunities(habitat) == countsubcommunities(budget) ||
      error("Habitat and budget must have same dimensions")
    countsubcommunities(habitat) == length(names) ||
      error("Number of subcommunities must match subcommunity names")
    return new{H, B}(habitat, active, budget, names)
  end
end
"""
    simplenicheAE(numniches::Int64, dimension::Tuple,
                        maxBud::Float64, area::Unitful.Area{Float64},
                        active::Array{Bool, 2})

Function to create a simple `DiscreteHab`, `SimpleBudget` type abiotic environment. Given a
number of niche types `numniches`, it creates a `DiscreteHab` environment with
dimensions `dimension` and a specified area `area`. It also creates a
`SimpleBudget` type filled with the maximum budget value `maxbud`. If a Bool
matrix of active grid squares is included, `active`, this is used, else one is
created with all grid cells active.
"""
function simplenicheAE(numniches::Int64, dimension::Tuple,
                        maxbud::Float64, area::Unitful.Area{Float64},
                        active::Array{Bool, 2})
  # Create niches
  niches = collect(1:numniches)
  area = uconvert(km^2, area)
  gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
  # Create niche-like environment
  hab = randomniches(dimension, niches, 0.5, fill(1.0/numniches, numniches),
  gridsquaresize)
  # Create empty budget and for now fill with one value
  bud = zeros(dimension)
  fill!(bud, maxbud/(dimension[1]*dimension[2]))
  return GridAbioticEnv{typeof(hab), SimpleBudget}(hab, active, SimpleBudget(bud))
end
function simplenicheAE(numniches::Int64, dimension::Tuple,
                        maxbud::Float64, area::Unitful.Area{Float64})
    active = Array{Bool,2}(Compat.undef, dimension)
    fill!(active, true)
    simplenicheAE(numniches, dimension, maxbud, area, active)
end
import Diversity.API: _countsubcommunities
function _countsubcommunities(gae::GridAbioticEnv)
  return countsubcommunities(gae.habitat)
end
import Diversity.API: _getsubcommunitynames
function _getsubcommunitynames(gae::GridAbioticEnv)
    return gae.names
end

function getavailableenergy(gae::GridAbioticEnv)
    return _getavailableenergy(gae.budget)
end
"""
    tempgradAE(min::Unitful.Temperature{Float64},
      max::Unitful.Temperature{Float64},
      dimension::Tuple{Int64, Int64}, maxbud::Float64,
      area::Unitful.Area{Float64}, rate::Quantity{Float64, typeof(𝚯*𝐓^-1)},
      active::Array{Bool, 2})

Function to create a temperature gradient `ContinuousHab`, `SimpleBudget` type abiotic
environment. Given a `min` and `max` temperature, it generates a
gradient from minimum at the bottom to maximum at the top. It creates a
`ContinuousHab` environment with dimensions `dimension` and a specified area
`area`. It also creates a `SimpleBudget` type filled with the maximum budget
value `maxbud`. The rate of temperature change is specified using the parameter
`rate`. If a Bool matrix of active grid squares is included, `active`,
this is used, else one is created with all grid cells active.
"""
function tempgradAE(minT::Unitful.Temperature{Float64},
  maxT::Unitful.Temperature{Float64},
  dimension::Tuple{Int64, Int64}, maxbud::Float64,
  area::Unitful.Area{Float64}, rate::Quantity{Float64, typeof(𝚯*𝐓^-1)},
  active::Array{Bool, 2})
  minT = uconvert(°C, minT); maxT = uconvert(°C, maxT)
  area = uconvert(km^2, area)
  gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
  hab = tempgrad(min, max, gridsquaresize, dimension, rate)
  bud = zeros(dimension)
  fill!(bud, maxbud/(dimension[1]*dimension[2]))

  return GridAbioticEnv{typeof(hab), SimpleBudget}(hab, active, SimpleBudget(bud))
end

function tempgradAE(minT::Unitful.Temperature{Float64},
  maxT::Unitful.Temperature{Float64},
  dimension::Tuple{Int64, Int64}, maxbud::Float64,
  area::Unitful.Area{Float64}, rate::Quantity{Float64, typeof(𝚯*𝐓^-1)})

  active = fill(true, dimension)
  tempgradAE(minT, maxT, dimension, maxbud, area, rate, active)
 end

 function eraAE(era::ERA, maxbud::Float64)
    dimension = size(era.array)[1:2]
    gridsquaresize = era.array.axes[1].val[2] - era.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    active = fill(true, dimension)
    active[isnan.(era.array[:,:,1])] .= false

    hab = ContinuousTimeHab(Array(era.array), 1, gridsquaresize,
        HabitatUpdate{Unitful.Dimensions{()}}(eraChange, 0.0/s))
    bud = zeros(dimension)
    fill!(bud, maxbud/(dimension[1]*dimension[2]))

     return GridAbioticEnv{typeof(hab), SimpleBudget}(hab, active, SimpleBudget(bud))
 end
 function eraAE(era::ERA, maxbud::Float64, active::Array{Bool, 2})
    dimension = size(era.array)[1:2]
    gridsquaresize = era.array.axes[1].val[2] - era.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = ContinuousTimeHab(Array(era.array), 1, gridsquaresize,
        HabitatUpdate{Unitful.Dimensions{()}}(eraChange, 0.0/s))
    bud = zeros(dimension)
    fill!(bud, maxbud/(dimension[1]*dimension[2]))

     return GridAbioticEnv{typeof(hab), SimpleBudget}(hab, active, SimpleBudget(bud))
 end
 function eraAE(era::ERA, bud::SolarBudget, active::Array{Bool, 2})
    dimension = size(era.array)[1:2]
    gridsquaresize = era.array.axes[1].val[2] - era.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = ContinuousTimeHab(Array(era.array), 1, gridsquaresize,
        HabitatUpdate{Unitful.Dimensions{()}}(eraChange, 0.0/s))

     return GridAbioticEnv{typeof(hab), SolarBudget}(hab, active, bud)
 end
 function worldclimAE(wc::Worldclim, maxbud::Float64)
    dimension = size(wc.array)[1:2]
    gridsquaresize = wc.array.axes[1].val[2] - wc.array.axes[1].val[1]

    active = fill(true, dimension)
    active[isnan.(wc.array[:,:,1])] = false

    hab = ContinuousTimeHab(Array(wc.array), 1, gridsquaresize,
        HabitatUpdate{Unitful.Dimensions{()}}(worldclimChange, 0.0/s))
    bud = zeros(dimension)
    fill!(bud, maxbud/(dimension[1]*dimension[2]))

     return GridAbioticEnv{typeof(hab), SimpleBudget}(hab, active, SimpleBudget(bud))
 end
 function worldclimAE(wc::Worldclim, maxbud::Float64, active::Array{Bool, 2})
    dimension = size(wc.array)[1:2]
    gridsquaresize = wc.array.axes[1].val[2] - wc.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = ContinuousTimeHab(Array(wc.array), 1, gridsquaresize,
        HabitatUpdate{Unitful.Dimensions{()}}(worldclimChange, 0.0/s))
    bud = zeros(dimension)
    fill!(bud, maxbud/(dimension[1]*dimension[2]))

     return GridAbioticEnv{typeof(hab), SimpleBudget}(hab, active, SimpleBudget(bud))
 end
 function worldclimAE(wc::Worldclim, bud::SolarBudget, active::Array{Bool, 2})
    dimension = size(wc.array)[1:2]
    gridsquaresize = wc.array.axes[1].val[2] - wc.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = ContinuousTimeHab(Array(wc.array), 1, gridsquaresize,
        HabitatUpdate{Unitful.Dimensions{()}}(worldclimChange, 0.0/s))

     return GridAbioticEnv{typeof(hab), SolarBudget}(hab, active, bud)
 end
 """
     simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
         dimension::Tuple{Int64, Int64}, maxbud::Float64, area::Unitful.Area{Float64},
         active::Array{Bool, 2})

 Function to create a simple `ContinuousHab`, `SimpleBudget` type abiotic
 environment. It creates a `ContinuousHab` filled with a given value, `val`,
 dimensions (`dimension`) and a specified area (`area`). It also creates a
 `SimpleBudget` type filled with the maximum budget value (`maxbud`).
 The rate of temperature change is specified using the parameter
 `rate`. If a Bool matrix of active grid squares is included, `active`,
 this is used, else one is created with all grid cells active.
 """
function simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
  dimension::Tuple{Int64, Int64}, maxbud::Float64, area::Unitful.Area{Float64},
  active::Array{Bool, 2})
  area = uconvert(km^2, area)
  gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
  hab = simplehabitat(val, gridsquaresize, dimension)
  bud = zeros(dimension)
  fill!(bud, maxbud/(dimension[1]*dimension[2]))

  return GridAbioticEnv{typeof(hab), SimpleBudget}(hab, active, SimpleBudget(bud))
end


function simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
  dimension::Tuple{Int64, Int64}, maxbud::Float64, area::Unitful.Area{Float64})

  active = fill(true, dimension)
  simplehabitatAE(val, dimension, maxbud, area, active)
end
