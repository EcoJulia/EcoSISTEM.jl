# SPDX-License-Identifier: LGPL-3.0-or-later

using Diversity
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref

using Diversity.API

matchdict = Dict(kJ => SolarBudget, mm => WaterBudget, NoUnits => SimpleBudget,
                 m^3 => VolWaterBudget)
checkbud(maxbud) = unit(maxbud) in keys(matchdict)
function cancel(a::Quantity{<:Real, 𝐌 * 𝐓^-2}, b::Quantity{<:Real, 𝐋^2})
    return uconvert(kJ, a * b)
end
function cancel(a::Quantity{<:Real, 𝐋 * 𝐋^-2}, b::Quantity{<:Real, 𝐋^2})
    return uconvert(mm, a * b)
end
function cancel(a::Quantity{<:Real, 𝐋^-2}, b::Quantity{<:Real, 𝐋^2})
    return uconvert(NoUnits, a * b)
end
function cancel(a::Quantity{<:Real, 𝐋^3 * 𝐋^-2}, b::Quantity{<:Real, 𝐋^2})
    return uconvert(m^3, a * b)
end
"""
    AbstractAbiotic{H <: AbstractHabitat, B <: AbstractBudget} <: AbstractPartition

Abstract supertype for all abiotic environment types and a subtype of
AbstractPartition
"""
abstract type AbstractAbiotic{H <: AbstractHabitat, B <: AbstractBudget} <:
              AbstractPartition{H} end

"""
    GridAbioticEnv{H, B} <: AbstractAbiotic{H, B}

This abiotic environment type holds a habitat and budget, as well as a string
of subcommunity names.
"""
mutable struct GridAbioticEnv{H, B} <: AbstractAbiotic{H, B}
    habitat::H
    active::Matrix{Bool}
    budget::B
    names::Vector{String}
    function (::Type{GridAbioticEnv{H, B}})(habitat::H, active::Matrix{Bool},
                                            budget::B,
                                            names::Vector{String} =
                                            map(x -> "$x",
                                                1:countsubcommunities(habitat))) where {
                                                                                        H,
                                                                                        B
                                                                                        }
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
                        active::Matrix{Bool})

Function to create a simple `DiscreteHab`, `SimpleBudget` type abiotic environment. Given a
number of niche types `numniches`, it creates a `DiscreteHab` environment with
dimensions `dimension` and a specified area `area`. It also creates a
`SimpleBudget` type filled with the maximum budget value `maxbud`. If a Bool
matrix of active grid squares is included, `active`, this is used, else one is
created with all grid cells active.
"""
function simplenicheAE(numniches::Int64, dimension::Tuple,
                       maxbud::Unitful.Quantity{Float64},
                       area::Unitful.Area{Float64},
                       active::Matrix{Bool})
    # Create niches
    niches = collect(1:numniches)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    # Create niche-like environment
    hab = randomniches(dimension, niches, 0.5, fill(1.0 / numniches, numniches),
                       gridsquaresize)
    # Create empty budget and for now fill with one value
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end
function simplenicheAE(numniches::Int64, dimension::Tuple,
                       maxbud::Unitful.Quantity{Float64},
                       area::Unitful.Area{Float64})
    active = Matrix{Bool}(undef, dimension)
    fill!(active, true)
    return simplenicheAE(numniches, dimension, maxbud, area, active)
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
      active::Matrix{Bool})

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
                    dimension::Tuple{Int64, Int64},
                    maxbud::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                    active::Matrix{Bool})
    minT = uconvert(K, minT)
    maxT = uconvert(K, maxT)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = tempgrad(minT, maxT, gridsquaresize, dimension, rate)
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end

function tempgradAE(minT::Unitful.Temperature{Float64},
                    maxT::Unitful.Temperature{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxbud::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝚯 * 𝐓^-1})
    active = fill(true, dimension)
    return tempgradAE(minT, maxT, dimension, maxbud, area, rate, active)
end

"""
    peakedgradAE(minT::Unitful.Temperature{Float64},
       maxT::Unitful.Temperature{Float64},
       dimension::Tuple{Int64, Int64}, maxbud::Unitful.Quantity{Float64},
       area::Unitful.Area{Float64}, rate::Quantity{Float64, 𝚯*𝐓^-1},
       active::Matrix{Bool})

Function to create a temperature gradient `ContinuousHab`, `SimpleBudget` type abiotic
environment. Given a `min` and `max` temperature, it generates a
gradient from minima at the top and bottom peaking to maximum in the middle. It creates a
`ContinuousHab` environment with dimensions `dimension` and a specified area
`area`. It also creates a `SimpleBudget` type filled with the maximum budget
value `maxbud`. The rate of temperature change is specified using the parameter
`rate`. If a Bool matrix of active grid squares is included, `active`,
this is used, else one is created with all grid cells active.
"""
function peakedgradAE(minT::Unitful.Temperature{Float64},
                      maxT::Unitful.Temperature{Float64},
                      dimension::Tuple{Int64, Int64},
                      maxbud::Unitful.Quantity{Float64},
                      area::Unitful.Area{Float64},
                      rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                      active::Matrix{Bool})
    minT = uconvert(K, minT)
    maxT = uconvert(K, maxT)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = tempgrad(minT, maxT + (maxT - minT), gridsquaresize, dimension, rate)
    hab.matrix[(ceil(Int, dimension[1] / 2) + 1):end, :] = hab.matrix[floor(Int,
                                                                            dimension[1] /
                                                                            2):-1:1,
                                                                      :]
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end

function peakedgradAE(minT::Unitful.Temperature{Float64},
                      maxT::Unitful.Temperature{Float64},
                      dimension::Tuple{Int64, Int64},
                      maxbud::Unitful.Quantity{Float64},
                      area::Unitful.Area{Float64},
                      rate::Quantity{Float64, 𝚯 * 𝐓^-1})
    active = fill(true, dimension)
    return peakedgradAE(minT, maxT, dimension, maxbud, area, rate, active)
end

"""
    raingradAE(min::Unitful.Temperature{Float64},
      max::Unitful.Temperature{Float64},
      dimension::Tuple{Int64, Int64}, maxbud::Float64,
      area::Unitful.Area{Float64}, rate::Quantity{Float64, typeof(𝚯*𝐓^-1)},
      active::Matrix{Bool})

Function to create a rain gradient `ContinuousHab`, `SimpleBudget` type abiotic environment. Given a `min` and `max` rainfall, it generates a
gradient from minimum at the bottom to maximum at the top. It creates a
`ContinuousHab` environment with dimensions `dimension` and a specified area `area`. It also creates a `SimpleBudget` type filled with the maximum budget value `maxbud`. The rate of rainfall change is specified using the parameter `rate`. If a Bool matrix of active grid squares is included, `active`, this is used, else one is created with all grid cells active.
"""
function raingradAE(minR::Unitful.Length{Float64},
                    maxR::Unitful.Length{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxbud::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝐋 * 𝐓^-1},
                    active::Matrix{Bool})
    minR = uconvert(mm, minR)
    maxR = uconvert(mm, maxR)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = raingrad(minR, maxR, gridsquaresize, dimension, rate)
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end

function raingradAE(minR::Unitful.Length{Float64},
                    maxR::Unitful.Length{Float64},
                    dimension::Tuple{Int64, Int64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝐋 * 𝐓^-1},
                    active::Matrix{Bool})
    minR = uconvert(mm, minR)
    maxR = uconvert(mm, maxR)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = raingrad(minR, maxR, gridsquaresize, dimension, rate)
    bud = WaterBudget(hab.matrix)
    return GridAbioticEnv{typeof(hab), typeof(bud)}(hab, active, bud)
end

function raingradAE(minR::Unitful.Length{Float64},
                    maxR::Unitful.Length{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxbud::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝐋 * 𝐓^-1})
    active = fill(true, dimension)
    return raingradAE(minR, maxR, dimension, maxbud, area, rate, active)
end
function raingradAE(minR::Unitful.Length{Float64},
                    maxR::Unitful.Length{Float64},
                    dimension::Tuple{Int64, Int64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝐋 * 𝐓^-1})
    active = fill(true, dimension)
    return raingradAE(minR, maxR, dimension, area, rate, active)
end
"""
   eraAE(era::ERA, maxbud::Unitful.Quantity{Float64})

Function to create a `ContinuousHab`, `SimpleBudget` type abiotic environment from an ERA type climate. It either creates a `SimpleBudget` type filled with the maximum budget value `maxbud` or uses a provided budget of type `SolarTimeBudget`. If a Bool matrix of active grid squares is included, `active`, this is used, else one is created with all grid cells active.
"""
function eraAE(era::ERA, maxbud::Unitful.Quantity{Float64},
               area::Unitful.Area{Float64})
    dimension = size(era.array)[1:2]
    gridsquaresize = era.array.axes[1].val[2] - era.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    active = fill(true, dimension)
    active[isnan.(era.array[:, :, 1])] .= false

    hab = ContinuousTimeHab(Array(era.array), 1, gridsquaresize,
                            HabitatUpdate(eraChange, 0.0 / s,
                                          Unitful.Dimensions{()}))
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end
function eraAE(era::ERA, maxbud::Unitful.Quantity{Float64},
               area::Unitful.Area{Float64}, active::Matrix{Bool})
    dimension = size(era.array)[1:2]
    gridsquaresize = era.array.axes[1].val[2] - era.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = ContinuousTimeHab(Array(era.array), 1, gridsquaresize,
                            HabitatUpdate(eraChange, 0.0 / s,
                                          Unitful.Dimensions{()}))
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end
function eraAE(era::ERA, bud::B,
               active::Matrix{Bool}) where {B <: AbstractTimeBudget}
    dimension = size(era.array)[1:2]
    gridsquaresize = era.array.axes[1].val[2] - era.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = ContinuousTimeHab(Array(era.array), 1, gridsquaresize,
                            HabitatUpdate(eraChange, 0.0 / s,
                                          Unitful.Dimensions{()}))

    return GridAbioticEnv{typeof(hab), typeof(bud)}(hab, active, bud)
end

"""
   worldclimAE(wc::Worldclim_monthly, maxbud::Unitful.Quantity{Float64})

Function to create a `ContinuousTimeHab`, `SimpleBudget` type abiotic environment from a Wordclim type climate. 
It either creates a `SimpleBudget` type filled with the maximum budget value `maxbud` or uses a provided budget of type `SolarTimeBudget`. 
If a Bool matrix of active grid squares is included, `active`, this is used, otherwise one is all grid cells are considered active.
"""
function worldclimAE(wc::Worldclim_monthly, maxbud::Unitful.Quantity{Float64},
                     area::Unitful.Area{Float64})
    dimension = size(wc.array)[1:2]
    gridsquaresize = wc.array.axes[1].val[2] - wc.array.axes[1].val[1]

    active = fill(true, dimension)
    active[isnan.(wc.array[:, :, 1])] .= false

    hab = ContinuousTimeHab(Array(wc.array), 1, gridsquaresize,
                            HabitatUpdate(worldclimChange, 0.0 / s,
                                          Unitful.Dimensions{()}))
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end
function worldclimAE(wc::Worldclim_monthly, maxbud::Unitful.Quantity{Float64},
                     area::Unitful.Area{Float64}, active::Matrix{Bool})
    dimension = size(wc.array)[1:2]
    gridsquaresize = wc.array.axes[1].val[2] - wc.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = ContinuousTimeHab(Array(wc.array), 1, gridsquaresize,
                            HabitatUpdate(worldclimChange, 0.0 / s,
                                          Unitful.Dimensions{()}))
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end
function worldclimAE(wc::Worldclim_monthly, bud::B,
                     active::Matrix{Bool}) where {B <: AbstractTimeBudget}
    gridsquaresize = wc.array.axes[1].val[2] - wc.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = ContinuousTimeHab(Array(wc.array), 1, gridsquaresize,
                            HabitatUpdate(worldclimChange, 0.0 / s,
                                          Unitful.Dimensions{()}))

    return GridAbioticEnv{typeof(hab), typeof(bud)}(hab, active, bud)
end

"""
  bioclimAE(bc::Worldclim_bioclim, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area{Float64})

Function to create a `ContinuousHab`, `SimpleBudget` type abiotic environment from an Wordclim type climate. 
It either creates a `SimpleBudget` type filled with the maximum budget value `maxbud` or uses a provided budget of type `SolarBudget`. 
If a Bool matrix of active grid squares is included, `active`, this is used, else one is created with all grid cells active.
"""
function bioclimAE(bc::Worldclim_bioclim, maxbud::Unitful.Quantity{Float64},
                   area::Unitful.Area{Float64})
    dimension = size(bc.array)[1:2]
    gridsquaresize = bc.array.axes[1].val[2] - bc.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km

    active = fill(true, dimension)
    active[isnan.(bc.array[:, :, 1])] .= false

    hab = ContinuousHab(Array(bc.array), gridsquaresize,
                        HabitatUpdate(NoChange, 0.0 / s,
                                      Unitful.Dimensions{()}))
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end
function bioclimAE(bc::Worldclim_bioclim, maxbud::Unitful.Quantity{Float64},
                   area::Unitful.Area{Float64}, active::Matrix{Bool})
    dimension = size(bc.array)[1:2]
    gridsquaresize = bc.array.axes[1].val[2] - bc.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = ContinuousHab(Array(bc.array), gridsquaresize,
                        HabitatUpdate(NoChange, 0.0 / s,
                                      Unitful.Dimensions{()}))
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end
function bioclimAE(bc::Worldclim_bioclim, bud::B,
                   active::Matrix{Bool}) where {B <: AbstractBudget}
    gridsquaresize = bc.array.axes[1].val[2] - bc.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = ContinuousHab(Array(bc.array), gridsquaresize,
                        HabitatUpdate(NoChange, 0.0 / s,
                                      Unitful.Dimensions{()}))

    return GridAbioticEnv{typeof(hab), typeof(bud)}(hab, active, bud)
end

"""
    simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
        dimension::Tuple{Int64, Int64}, maxbud::Float64, area::Unitful.Area{Float64},
        active::Matrix{Bool})

Function to create a simple `ContinuousHab`, `SimpleBudget` type abiotic
environment. It creates a `ContinuousHab` filled with a given value, `val`,
dimensions (`dimension`) and a specified area (`area`). It also creates a
`SimpleBudget` type filled with the maximum budget value (`maxbud`).
The rate of temperature change is specified using the parameter
`rate`. If a Bool matrix of active grid squares is included, `active`,
this is used, else one is created with all grid cells active.
"""
function simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
                         dimension::Tuple{Int64, Int64},
                         maxbud::Unitful.Quantity{Float64},
                         area::Unitful.Area{Float64},
                         active::Matrix{Bool})
    if typeof(val) <: Unitful.Temperature
        val = uconvert(K, val)
    end
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = simplehabitat(val, gridsquaresize, dimension)
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end

function simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
                         dimension::Tuple{Int64, Int64},
                         maxbud::Unitful.Quantity{Float64},
                         area::Unitful.Area{Float64})
    active = fill(true, dimension)
    return simplehabitatAE(val, dimension, maxbud, area, active)
end

import EcoBase.getcoords

getcoords(abenv::GridAbioticEnv) = abenv.habitat

function lcAE(lc::Landcover, maxbud::Unitful.Quantity{Float64},
              area::Unitful.Area)
    dimension = size(lc.array)
    gridsquaresize = lc.array.axes[1].val[2] - lc.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    active = fill(true, dimension)
    active[isnan.(lc.array[:, :, 1])] .= false

    hab = DiscreteHab(Array(lc.array), gridsquaresize,
                      HabitatUpdate(NoChange, 0.0 / s))
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end
function lcAE(lc::Landcover, maxbud::Unitful.Quantity{Float64},
              area::Unitful.Area, active::Matrix{Bool})
    dimension = size(lc.array)[1:2]
    gridsquaresize = lc.array.axes[1].val[2] - lc.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = DiscreteHab(Array(lc.array), gridsquaresize,
                      HabitatUpdate(NoChange, 0.0 / s))
    B = cancel(maxbud, area)
    bud = zeros(typeof(B), dimension)
    fill!(bud, B / (dimension[1] * dimension[2]))
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = matchdict[unit(B)]
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end
function lcAE(lc::Landcover, bud::B,
              active::Matrix{Bool}) where {B <: AbstractBudget}
    gridsquaresize = lc.array.axes[1].val[2] - lc.array.axes[1].val[1]
    gridsquaresize = ustrip.(gridsquaresize) * 111.32km
    hab = DiscreteHab(Array(lc.array), gridsquaresize,
                      HabitatUpdate(NoChange, 0.0 / s))

    return GridAbioticEnv{typeof(hab), B}(hab, active, bud)
end
