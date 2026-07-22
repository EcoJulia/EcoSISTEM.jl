# SPDX-License-Identifier: LGPL-3.0-or-later

using Diversity
using Unitful
using Unitful.DefaultSymbols
using Diversity.API
using RasterDataSources
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref

const RDS = RasterDataSources
const BUDGETDICT = Dict(kJ => SolarBudget, mm => WaterBudget,
                        NoUnits => SimpleBudget,
                        m^3 => VolWaterBudget)
checkbud(maxbud) = unit(maxbud) in keys(BUDGETDICT)
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

# Shared tail of the maximum-budget `*AE` constructors: total the per-area maximum budget
# `maxbud` over `area`, spread it uniformly across the grid of `active`, pick the budget
# type from its units, and assemble the `GridAbioticEnv` with habitat `hab`.
function _maxbudget_env(hab, active::Matrix{Bool},
                        maxbud::Unitful.Quantity{Float64}, area::Unitful.Area)
    B = cancel(maxbud, area)
    checkbud(B) || error("Unrecognised unit in budget")
    budtype = BUDGETDICT[unit(B)]
    bud = fill(B / length(active), size(active))
    return GridAbioticEnv{typeof(hab), budtype}(hab, active, budtype(bud))
end

# The time-varying (monthly) climate habitat shared by `eraAE`/`worldclimAE`: a
# `ContinuousTimeHab` starting at month 1 that advances via `cyclicChange`.
function _timeclimate_hab(array)
    return ContinuousTimeHab(Array(array), 1, _cellsize(array),
                             HabitatUpdate(cyclicChange, 0.0 / s,
                                           Unitful.Dimensions{()}))
end

# Static (no-dynamics) data habitats shared by `bioclimAE` (continuous) and `lcAE`
# (discrete land cover).
function _continuoushab(array)
    return ContinuousHab(Array(array), _cellsize(array),
                         HabitatUpdate(NoChange, 0.0 / s,
                                       Unitful.Dimensions{()}))
end
function _discretehab(array)
    return DiscreteHab(Array(array), _cellsize(array),
                       HabitatUpdate(NoChange, 0.0 / s))
end

# Default active mask for a data habitat: every cell active except those that are NaN in
# the first layer. Returns a `Matrix{Bool}` (matching `_maxbudget_env`).
function _nanactive(array)
    active = fill(true, size(array)[1:2])
    active[isnan.(array[:, :, 1])] .= false
    return active
end

# Length of a degree of latitude (≈ constant); a degree of longitude is this
# times cos(latitude).
const LONGITUDE_DEGREE_LENGTH = 111.32km / °

"""
    _cellsize(lats, longs)

Area-preserving physical side of a grid cell whose cell-centre latitudes and
longitudes are `lats` and `longs`. The behaviour depends on the units of the
coordinates (multiple dispatch):

  - **Geographic (angle) coordinates** (degrees): the north–south side is
    `Δlat × 111.32 km`; the east–west side uses the true length of a degree of
    longitude at the grid's centre latitude, `Δlong × 111.32 km × cos φ`. Returns
    the geometric mean of the two so that each cell's area is geodetically correct,
    and reports (`@info`) how much the east–west length varies from the top to the
    bottom of the grid.
  - **Projected (length) coordinates** (`m`, `km`, …): the grid is already metric,
    so there is no spherical adjustment — the cell side is simply
    `sqrt(Δlat × Δlong)`.
"""
function _cellsize(lats, longs)
    dlat = abs((lats[2] - lats[1]))
    dlong = abs((longs[2] - longs[1]))
    φtop, φbot = maximum(lats), minimum(lats)
    ns = dlat * LONGITUDE_DEGREE_LENGTH
    ew = dlong * LONGITUDE_DEGREE_LENGTH * cos((φtop + φbot) / 2)
    ewtop = dlong * LONGITUDE_DEGREE_LENGTH * cos(φtop)
    ewbot = dlong * LONGITUDE_DEGREE_LENGTH * cos(φbot)
    isapprox(ewtop, ewbot; rtol = 1.0e-2) ||
        @info "East–west cell length varies with latitude across this grid: " *
              "$(round(typeof(1.0km), ewtop, digits = 2)) at $(round(typeof(1.0°), φtop, digits = 1)) (top), " *
              "$(round(typeof(1.0km), ew, digits = 2)) at the centre, " *
              "$(round(typeof(1.0km), ewbot, digits = 2)) at $(round(typeof(1.0°), φbot, digits = 1)) (bottom); " *
              "using the area-preserving cell size $(round(typeof(1.0km), sqrt(ns * ew), digits = 2))."
    return sqrt(ns * ew)
end
# Projected (length) coordinates: already metric, so no spherical adjustment. The
# result keeps its native length unit; `genlookups` makes the cell-size/dispersal
# ratio dimensionless explicitly, so no unit normalisation is needed here.
function _cellsize(lats::AbstractVector{<:Unitful.Length},
                   longs::AbstractVector{<:Unitful.Length})
    return sqrt(abs(lats[2] - lats[1]) * abs(longs[2] - longs[1]))
end
_cellsize(A) = _cellsize(A.axes[1].val, A.axes[2].val)

"""
    AbstractAbiotic{H <: AbstractHabitat, B <: AbstractBudget} <: AbstractPartition

Abstract supertype for all abiotic environment types and a subtype of
AbstractPartition
"""
abstract type AbstractAbiotic{H <: AbstractHabitat, B <: AbstractBudget} <:
              AbstractPartition{H} end

"""
    GridAbioticEnv{H, B} <: AbstractAbiotic{H, B}

Abiotic environment type holding a `habitat` of type `H`, a boolean `active`
matrix indicating which grid cells are accessible, a `budget` of type `B`
representing available resources, and a vector of `names` for each subcommunity.
"""
struct GridAbioticEnv{H, B} <: AbstractAbiotic{H, B}
    habitat::H
    active::Matrix{Bool}
    budget::B
    names::Vector{String}
    function (::Type{GridAbioticEnv{H, B}})(habitat::H,
                                            active::Matrix{Bool},
                                            budget::B,
                                            names::Vector{String} = map(x -> "$x",
                                                                        1:countsubcommunities(habitat))) where {H,
                                                                                                                B}
        countsubcommunities(habitat) == countsubcommunities(budget) ||
            error("Habitat and budget must have same dimensions")
        countsubcommunities(habitat) == length(names) ||
            error("Number of subcommunities must match subcommunity names")
        return new{H, B}(habitat, active, budget, names)
    end
end

"""
    simplenicheAE(numniches::Int64, dimension::Tuple,
                        maxbud::Unitful.Quantity{Float64}, area::Unitful.Area{Float64},
                        active::Matrix{Bool})

Create a simple [`DiscreteHab`](@ref), [`SimpleBudget`](@ref) type abiotic
environment. Given a number of niche types `numniches`, it creates a
[`DiscreteHab`](@ref) environment with dimensions `dimension` and a specified
area `area`. It also creates a [`SimpleBudget`](@ref) type filled with the
maximum budget value `maxbud`. If a Bool matrix of active grid squares is
included, `active`, this is used, else one is created with all grid cells
active.
"""
function simplenicheAE(numniches::Int64,
                       dimension::Tuple,
                       maxbud::Unitful.Quantity{Float64},
                       area::Unitful.Area{Float64},
                       active::Matrix{Bool};
                       axis::Type{<:NicheAxis} = Unclassified)
    # Create niches
    niches = collect(1:numniches)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    # Create niche-like environment, tagged with its niche axis
    hab = _reaxis(randomniches(dimension,
                               niches,
                               0.5,
                               fill(1.0 / numniches, numniches),
                               gridsquaresize), axis)
    # Create empty budget and for now fill with one value
    return _maxbudget_env(hab, active, maxbud, area)
end
function simplenicheAE(numniches::Int64,
                       dimension::Tuple,
                       maxbud::Unitful.Quantity{Float64},
                       area::Unitful.Area{Float64};
                       axis::Type{<:NicheAxis} = Unclassified)
    active = Matrix{Bool}(undef, dimension)
    fill!(active, true)
    return simplenicheAE(numniches, dimension, maxbud, area, active;
                         axis = axis)
end
@doc (@doc simplenicheAE) simplenicheAE(::Int64,
                                        ::Tuple,
                                        ::Unitful.Quantity{Float64},
                                        ::Unitful.Area{Float64})

import Diversity.API: _countsubcommunities
function _countsubcommunities(gae::GridAbioticEnv)
    return countsubcommunities(gae.habitat)
end
import Diversity.API: _getsubcommunitynames
function _getsubcommunitynames(gae::GridAbioticEnv)
    return gae.names
end

"""
    getavailableenergy(gae::GridAbioticEnv)

Return the available energy budget from a [`GridAbioticEnv`](@ref).
"""
function getavailableenergy(gae::GridAbioticEnv)
    return _getavailableenergy(gae.budget)
end
"""
    tempgradAE(minT::Unitful.Temperature{Float64},
      maxT::Unitful.Temperature{Float64},
      dimension::Tuple{Int64, Int64}, maxbud::Unitful.Quantity{Float64},
      area::Unitful.Area{Float64}, rate::Quantity{Float64, 𝚯*𝐓^-1},
      active::Matrix{Bool})

Create a temperature gradient [`ContinuousHab`](@ref), [`SimpleBudget`](@ref)
type abiotic environment. Given a `minT` and `maxT` temperature, it generates a
gradient from minimum at the bottom to maximum at the top. It creates a
[`ContinuousHab`](@ref) environment with dimensions `dimension` and a specified
area `area`. It also creates a [`SimpleBudget`](@ref) type filled with the
maximum budget value `maxbud`. The rate of temperature change is specified using
the parameter `rate`. If a Bool matrix of active grid squares is included,
`active`, this is used, else one is created with all grid cells active.
"""
function tempgradAE(minT::Unitful.Temperature{Float64},
                    maxT::Unitful.Temperature{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxbud::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                    active::Matrix{Bool};
                    axis::Type{<:NicheAxis} = MeanTemperature)
    minT = _canonical(minT, axis)
    maxT = _canonical(maxT, axis)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = _reaxis(tempgrad(minT, maxT, gridsquaresize, dimension, rate), axis)
    return _maxbudget_env(hab, active, maxbud, area)
end

function tempgradAE(minT::Unitful.Temperature{Float64},
                    maxT::Unitful.Temperature{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxbud::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝚯 * 𝐓^-1};
                    axis::Type{<:NicheAxis} = MeanTemperature)
    active = fill(true, dimension)
    return tempgradAE(minT, maxT, dimension, maxbud, area, rate, active;
                      axis = axis)
end
@doc (@doc tempgradAE) tempgradAE(::Unitful.Temperature{Float64},
                                  ::Unitful.Temperature{Float64},
                                  ::Tuple{Int64, Int64},
                                  ::Unitful.Quantity{Float64},
                                  ::Unitful.Area{Float64},
                                  ::Quantity{Float64, 𝚯 * 𝐓^-1})

"""
    peakedgradAE(minT::Unitful.Temperature{Float64},
       maxT::Unitful.Temperature{Float64},
       dimension::Tuple{Int64, Int64}, maxbud::Unitful.Quantity{Float64},
       area::Unitful.Area{Float64}, rate::Quantity{Float64, 𝚯*𝐓^-1},
       active::Matrix{Bool})

Create a peaked temperature gradient [`ContinuousHab`](@ref),
[`SimpleBudget`](@ref) type abiotic environment. Given a `minT` and `maxT`
temperature, it generates a gradient with minima at the top and bottom, peaking
at `maxT` in the middle. It creates a [`ContinuousHab`](@ref) environment with
dimensions `dimension` and a specified area `area`. It also creates a
[`SimpleBudget`](@ref) type filled with the maximum budget value `maxbud`. The
rate of temperature change is specified using the parameter `rate`. If a Bool
matrix of active grid squares is included, `active`, this is used, else one is
created with all grid cells active.
"""
function peakedgradAE(minT::Unitful.Temperature{Float64},
                      maxT::Unitful.Temperature{Float64},
                      dimension::Tuple{Int64, Int64},
                      maxbud::Unitful.Quantity{Float64},
                      area::Unitful.Area{Float64},
                      rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                      active::Matrix{Bool};
                      axis::Type{<:NicheAxis} = MeanTemperature)
    minT = _canonical(minT, axis)
    maxT = _canonical(maxT, axis)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = tempgrad(minT, maxT + (maxT - minT), gridsquaresize, dimension, rate)
    hab.matrix[(ceil(Int, dimension[1] / 2) + 1):end, :] = hab.matrix[floor(Int,
                                                                            dimension[1] / 2):-1:1,
                                                                      :]
    return _maxbudget_env(_reaxis(hab, axis), active, maxbud, area)
end

function peakedgradAE(minT::Unitful.Temperature{Float64},
                      maxT::Unitful.Temperature{Float64},
                      dimension::Tuple{Int64, Int64},
                      maxbud::Unitful.Quantity{Float64},
                      area::Unitful.Area{Float64},
                      rate::Quantity{Float64, 𝚯 * 𝐓^-1};
                      axis::Type{<:NicheAxis} = MeanTemperature)
    active = fill(true, dimension)
    return peakedgradAE(minT, maxT, dimension, maxbud, area, rate, active;
                        axis = axis)
end
@doc (@doc peakedgradAE) peakedgradAE(::Unitful.Temperature{Float64},
                                      ::Unitful.Temperature{Float64},
                                      ::Tuple{Int64, Int64},
                                      ::Unitful.Quantity{Float64},
                                      ::Unitful.Area{Float64},
                                      ::Quantity{Float64, 𝚯 * 𝐓^-1})

"""
    raingradAE(minR::Unitful.Length{Float64},
      maxR::Unitful.Length{Float64},
      dimension::Tuple{Int64, Int64}, maxbud::Unitful.Quantity{Float64},
      area::Unitful.Area{Float64}, rate::Quantity{Float64, 𝐋*𝐓^-1},
      active::Matrix{Bool})

Create a rainfall gradient [`ContinuousHab`](@ref), [`SimpleBudget`](@ref) type
abiotic environment. Given a `minR` and `maxR` rainfall, it generates a gradient
from minimum at the bottom to maximum at the top. It creates a
[`ContinuousHab`](@ref) environment with dimensions `dimension` and a specified
area `area`. It also creates a [`SimpleBudget`](@ref) type filled with the
maximum budget value `maxbud`. The rate of rainfall change is specified using
the parameter `rate`. If a Bool matrix of active grid squares is included,
`active`, this is used, else one is created with all grid cells active.
"""
function raingradAE(minR::Unitful.Length{Float64},
                    maxR::Unitful.Length{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxbud::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝐋 * 𝐓^-1},
                    active::Matrix{Bool};
                    axis::Type{<:NicheAxis} = Precipitation)
    minR = _canonical(minR, axis)
    maxR = _canonical(maxR, axis)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = _reaxis(raingrad(minR, maxR, gridsquaresize, dimension, rate), axis)
    return _maxbudget_env(hab, active, maxbud, area)
end

"""
    raingradAE(minR::Unitful.Length{Float64}, maxR::Unitful.Length{Float64},
      dimension::Tuple{Int64, Int64}, area::Unitful.Area{Float64},
      rate::Quantity{Float64, 𝐋*𝐓^-1}, active::Matrix{Bool})

As [`raingradAE`](@ref) but uses the rainfall values from the gradient directly
as a [`WaterBudget`](@ref), rather than computing a budget from a separate
`maxbud` value.
"""
function raingradAE(minR::Unitful.Length{Float64},
                    maxR::Unitful.Length{Float64},
                    dimension::Tuple{Int64, Int64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝐋 * 𝐓^-1},
                    active::Matrix{Bool};
                    axis::Type{<:NicheAxis} = Precipitation)
    minR = _canonical(minR, axis)
    maxR = _canonical(maxR, axis)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = _reaxis(raingrad(minR, maxR, gridsquaresize, dimension, rate), axis)
    bud = WaterBudget(hab.matrix)
    return GridAbioticEnv{typeof(hab), typeof(bud)}(hab, active, bud)
end

function raingradAE(minR::Unitful.Length{Float64},
                    maxR::Unitful.Length{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxbud::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝐋 * 𝐓^-1};
                    axis::Type{<:NicheAxis} = Precipitation)
    active = fill(true, dimension)
    return raingradAE(minR, maxR, dimension, maxbud, area, rate, active;
                      axis = axis)
end
@doc (@doc raingradAE) raingradAE(::Unitful.Length{Float64},
                                  ::Unitful.Length{Float64},
                                  ::Tuple{Int64, Int64},
                                  ::Unitful.Quantity{Float64},
                                  ::Unitful.Area{Float64},
                                  ::Quantity{Float64, 𝐋 * 𝐓^-1})
function raingradAE(minR::Unitful.Length{Float64},
                    maxR::Unitful.Length{Float64},
                    dimension::Tuple{Int64, Int64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝐋 * 𝐓^-1};
                    axis::Type{<:NicheAxis} = Precipitation)
    active = fill(true, dimension)
    return raingradAE(minR, maxR, dimension, area, rate, active; axis = axis)
end
@doc (@doc raingradAE) raingradAE(::Unitful.Length{Float64},
                                  ::Unitful.Length{Float64},
                                  ::Tuple{Int64, Int64},
                                  ::Unitful.Area{Float64},
                                  ::Quantity{Float64, 𝐋 * 𝐓^-1})
"""
   eraAE(era::ERA, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area{Float64})

Create a [`ContinuousHab`](@ref), [`SimpleBudget`](@ref) type abiotic
environment from an ERA type climate. It either creates a [`SimpleBudget`](@ref)
type filled with the maximum budget value `maxbud` or uses a provided budget of
type [`SolarTimeBudget`](@ref). If a Bool matrix of active grid squares is
included, `active`, this is used, else one is created with all grid cells
active.
"""
function eraAE(era::ERA, maxbud::Unitful.Quantity{Float64},
               area::Unitful.Area{Float64})
    return _maxbudget_env(_timeclimate_hab(era.array), _nanactive(era.array),
                          maxbud, area)
end
"""
    eraAE(era::ERA, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area{Float64},
      active::Matrix{Bool})

As [`eraAE`](@ref) with an explicit `active` matrix of grid squares, rather than
inferring active cells from NaN values in the ERA data.
"""
function eraAE(era::ERA,
               maxbud::Unitful.Quantity{Float64},
               area::Unitful.Area{Float64},
               active::Matrix{Bool})
    return _maxbudget_env(_timeclimate_hab(era.array), active, maxbud, area)
end
"""
    eraAE(era::ERA, bud::B, active::Matrix{Bool}) where B <: AbstractTimeBudget

As [`eraAE`](@ref) but accepts a pre-constructed `AbstractTimeBudget` object
`bud` rather than computing a budget from a maximum value.
"""
function eraAE(era::ERA, bud::B,
               active::Matrix{Bool}) where {B <: AbstractTimeBudget}
    hab = _timeclimate_hab(era.array)
    return GridAbioticEnv{typeof(hab), typeof(bud)}(hab, active, bud)
end

"""
   worldclimAE(wc::ClimateRaster{WorldClim{Climate}}, maxbud::Unitful.Quantity{Float64},
     area::Unitful.Area{Float64})

Create a [`ContinuousTimeHab`](@ref), [`SimpleBudget`](@ref) type abiotic
environment from a Worldclim type climate. It either creates a
[`SimpleBudget`](@ref) type filled with the maximum budget value `maxbud` or
uses a provided budget of type [`SolarTimeBudget`](@ref). If a Bool matrix of
active grid squares is included, `active`, this is used, otherwise all grid
cells are considered active.
"""
function worldclimAE(wc::ClimateRaster{WorldClim{Climate}},
                     maxbud::Unitful.Quantity{Float64},
                     area::Unitful.Area{Float64})
    return _maxbudget_env(_timeclimate_hab(wc.array), _nanactive(wc.array),
                          maxbud, area)
end
"""
    worldclimAE(wc::ClimateRaster{WorldClim{Climate}}, maxbud::Unitful.Quantity{Float64},
      area::Unitful.Area{Float64}, active::Matrix{Bool})

As [`worldclimAE`](@ref) with an explicit `active` matrix of grid squares,
rather than inferring active cells from NaN values in the Worldclim data.
"""
function worldclimAE(wc::ClimateRaster{WorldClim{Climate}},
                     maxbud::Unitful.Quantity{Float64},
                     area::Unitful.Area{Float64},
                     active::Matrix{Bool})
    return _maxbudget_env(_timeclimate_hab(wc.array), active, maxbud, area)
end
"""
    worldclimAE(wc::ClimateRaster{WorldClim{Climate}}, bud::B, active::Matrix{Bool}) where B <: AbstractTimeBudget

As [`worldclimAE`](@ref) but accepts a pre-constructed `AbstractTimeBudget`
object `bud` rather than computing a budget from a maximum value.
"""
function worldclimAE(wc::ClimateRaster{WorldClim{Climate}},
                     bud::B,
                     active::Matrix{Bool}) where {B <: AbstractTimeBudget}
    hab = _timeclimate_hab(wc.array)
    return GridAbioticEnv{typeof(hab), typeof(bud)}(hab, active, bud)
end

"""
    bioclimAE(bc::ClimateRaster{WorldClim{BioClim}}, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area{Float64})

Create a [`ContinuousHab`](@ref), [`SimpleBudget`](@ref) type abiotic
environment from a Worldclim type climate. It either creates a
[`SimpleBudget`](@ref) type filled with the maximum budget value `maxbud` or
uses a provided budget of type [`SolarBudget`](@ref). If a Bool matrix of active
grid squares is included, `active`, this is used, else one is created with all
grid cells active.
"""
function bioclimAE(bc::ClimateRaster{WorldClim{BioClim}},
                   maxbud::Unitful.Quantity{Float64},
                   area::Unitful.Area{Float64})
    return _maxbudget_env(_continuoushab(bc.array), _nanactive(bc.array),
                          maxbud, area)
end
"""
    bioclimAE(bc::ClimateRaster{WorldClim{BioClim}}, maxbud::Unitful.Quantity{Float64},
      area::Unitful.Area{Float64}, active::Matrix{Bool})

As [`bioclimAE`](@ref) with an explicit `active` matrix of grid squares, rather
than inferring active cells from NaN values in the bioclim data.
"""
function bioclimAE(bc::ClimateRaster{WorldClim{BioClim}},
                   maxbud::Unitful.Quantity{Float64},
                   area::Unitful.Area{Float64}, active::Matrix{Bool})
    return _maxbudget_env(_continuoushab(bc.array), active, maxbud, area)
end
"""
    bioclimAE(bc::ClimateRaster{WorldClim{BioClim}}, bud::B, active::Matrix{Bool}) where B <: AbstractBudget

As [`bioclimAE`](@ref) but accepts a pre-constructed [`AbstractBudget`](@ref)
object `bud` rather than computing a budget from a maximum value.
"""
function bioclimAE(bc::ClimateRaster{WorldClim{BioClim}}, bud::B,
                   active::Matrix{Bool}) where {B <: AbstractBudget}
    hab = _continuoushab(bc.array)
    return GridAbioticEnv{typeof(hab), typeof(bud)}(hab, active, bud)
end

"""
    simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
        dimension::Tuple{Int64, Int64}, maxbud::Float64, area::Unitful.Area{Float64},
        active::Matrix{Bool})

Create a simple [`ContinuousHab`](@ref), [`SimpleBudget`](@ref) type abiotic
environment. It creates a [`ContinuousHab`](@ref) filled with a given value,
`val`, dimensions (`dimension`) and a specified area (`area`). It also creates a
[`SimpleBudget`](@ref) type filled with the maximum budget value (`maxbud`). If
a Bool matrix of active grid squares is included, `active`, this is used, else
one is created with all grid cells active.
"""
function simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
                         dimension::Tuple{Int64, Int64},
                         maxbud::Unitful.Quantity{Float64},
                         area::Unitful.Area{Float64},
                         active::Matrix{Bool};
                         axis::Type{<:NicheAxis} = Unclassified)
    val = _canonical(val, axis)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    hab = _reaxis(simplehabitat(val, gridsquaresize, dimension), axis)
    return _maxbudget_env(hab, active, maxbud, area)
end

function simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
                         dimension::Tuple{Int64, Int64},
                         maxbud::Unitful.Quantity{Float64},
                         area::Unitful.Area{Float64};
                         axis::Type{<:NicheAxis} = Unclassified)
    active = fill(true, dimension)
    return simplehabitatAE(val, dimension, maxbud, area, active; axis = axis)
end
@doc (@doc simplehabitatAE) simplehabitatAE(::Union{Float64,
                                                    Unitful.Quantity{Float64}},
                                            ::Tuple{Int64, Int64},
                                            ::Unitful.Quantity{Float64},
                                            ::Unitful.Area{Float64})

import EcoBase.getcoords

getcoords(abenv::GridAbioticEnv) = abenv.habitat

"""
    lcAE(lc::ClimateRaster{<:EarthEnv{<:LandCover}}, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area)

Create a [`DiscreteHab`](@ref), [`SimpleBudget`](@ref) type abiotic environment
from a land cover [`ClimateRaster`](@ref) dataset. It creates a
[`DiscreteHab`](@ref) habitat from the land cover array and a
[`SimpleBudget`](@ref) type filled with the maximum budget value `maxbud`, scaled
to the given `area`. If a Bool matrix of active grid squares is included,
`active`, this is used, else one is created with all grid cells active.
"""
function lcAE(lc::ClimateRaster{T, A}, maxbud::Unitful.Quantity{Float64},
              area::Unitful.Area) where {T <: EarthEnv{<:LandCover}, A}
    return _maxbudget_env(_discretehab(lc.array), _nanactive(lc.array),
                          maxbud, area)
end
"""
    lcAE(lc::ClimateRaster{<:EarthEnv{<:LandCover}}, maxbud::Unitful.Quantity{Float64}, area::Unitful.Area,
      active::Matrix{Bool})

As [`lcAE`](@ref) with an explicit `active` matrix of grid squares, rather than
setting all cells active.
"""
function lcAE(lc::ClimateRaster{T, A}, maxbud::Unitful.Quantity{Float64},
              area::Unitful.Area,
              active::Matrix{Bool}) where {T <: EarthEnv{<:LandCover}, A}
    return _maxbudget_env(_discretehab(lc.array), active, maxbud, area)
end
"""
    lcAE(lc::ClimateRaster{<:EarthEnv{<:LandCover}}, bud::B, active::Matrix{Bool}) where B <: AbstractBudget

As [`lcAE`](@ref) but accepts a pre-constructed [`AbstractBudget`](@ref) object
`bud` rather than computing a budget from a maximum value.
"""
function lcAE(lc::ClimateRaster{T, A}, bud::B,
              active::Matrix{Bool}) where {T <: EarthEnv{<:LandCover},
                                           A, B <: AbstractBudget}
    hab = _discretehab(lc.array)
    return GridAbioticEnv{typeof(hab), B}(hab, active, bud)
end
