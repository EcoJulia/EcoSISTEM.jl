# SPDX-License-Identifier: LGPL-3.0-or-later

using Diversity
using Unitful
using Unitful.DefaultSymbols
using Diversity.API
using RasterDataSources
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref

const RDS = RasterDataSources
const SUPPLYDICT = Dict(kJ => SolarSupply, mm => WaterSupply,
                        NoUnits => SimpleSupply,
                        m^3 => VolWaterSupply)
checksupply(maxsupply) = unit(maxsupply) in keys(SUPPLYDICT)
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

# Shared tail of the maximum-supply `*AE` constructors: total the per-area maximum supply
# `maxsupply` over `area`, spread it uniformly across the grid of `active`, pick the supply
# type from its units, and assemble the `GridHabitat` with regime `regime`.
function _maxsupply_env(regime, active::Matrix{Bool},
                        maxsupply::Unitful.Quantity{Float64},
                        area::Unitful.Area)
    B = cancel(maxsupply, area)
    checksupply(B) || error("Unrecognised unit in supply")
    supplytype = SUPPLYDICT[unit(B)]
    supply = fill(B / length(active), size(active))
    return GridHabitat{typeof(regime), supplytype}(regime, active,
                                                   supplytype(supply))
end

# The time-varying (monthly) climate regime shared by `eraAE`/`worldclimAE`: a
# `ContinuousTimeRegime` starting at month 1 that advances via `cyclicChange`.
function _timeclimate_regime(array)
    return ContinuousTimeRegime(Array(array), 1, _cellsize(array),
                                LayerUpdate(cyclicChange, 0.0 / s,
                                            Unitful.Dimensions{()}))
end

# Static (no-dynamics) data regimes shared by `bioclimAE` (continuous) and `lcAE`
# (discrete land cover).
function _continuousregime(array)
    return ContinuousRegime(Array(array), _cellsize(array),
                            LayerUpdate(NoChange, 0.0 / s,
                                        Unitful.Dimensions{()}))
end
function _discreteregime(array)
    return DiscreteRegime(Array(array), _cellsize(array),
                          LayerUpdate(NoChange, 0.0 / s))
end

# Default active mask for a data regime: every cell active except those that are NaN in
# the first layer. Returns a `Matrix{Bool}` (matching `_maxsupply_env`).
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
    AbstractHabitat{H <: AbstractRegime, B <: AbstractSupply} <: AbstractPartition

Abstract supertype for all abiotic environment types and a subtype of
AbstractPartition
"""
abstract type AbstractHabitat{H <: AbstractRegime, B <: AbstractSupply} <:
              AbstractPartition{H} end

"""
    GridHabitat{H, B} <: AbstractHabitat{H, B}

Abiotic environment type holding a `regime` of type `H`, a boolean `active`
matrix indicating which grid cells are accessible, a `supply` of type `B`
representing available resources, and a vector of `names` for each subcommunity.
"""
struct GridHabitat{H, B} <: AbstractHabitat{H, B}
    regime::H
    active::Matrix{Bool}
    supply::B
    names::Vector{String}
    function (::Type{GridHabitat{H, B}})(regime::H,
                                         active::Matrix{Bool},
                                         supply::B,
                                         names::Vector{String} = map(x -> "$x",
                                                                     1:countsubcommunities(regime))) where {H,
                                                                                                            B}
        countsubcommunities(regime) == countsubcommunities(supply) ||
            error("Condition and supply must have same dimensions")
        countsubcommunities(regime) == length(names) ||
            error("Number of subcommunities must match subcommunity names")
        return new{H, B}(regime, active, supply, names)
    end
end

"""
    simplenicheAE(numniches::Int64, dimension::Tuple,
                        maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area{Float64},
                        active::Matrix{Bool})

Create a simple [`DiscreteRegime`](@ref), [`SimpleSupply`](@ref) type abiotic
environment. Given a number of niche types `numniches`, it creates a
[`DiscreteRegime`](@ref) environment with dimensions `dimension` and a specified
area `area`. It also creates a [`SimpleSupply`](@ref) type filled with the
maximum supply value `maxsupply`. If a Bool matrix of active grid squares is
included, `active`, this is used, else one is created with all grid cells
active.
"""
function simplenicheAE(numniches::Int64,
                       dimension::Tuple,
                       maxsupply::Unitful.Quantity{Float64},
                       area::Unitful.Area{Float64},
                       active::Matrix{Bool};
                       axis::Type{<:NicheAxis} = Unclassified)
    # Create niches
    niches = collect(1:numniches)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    # Create niche-like environment, tagged with its niche axis
    regime = _reaxis(randomniches(dimension,
                                  niches,
                                  0.5,
                                  fill(1.0 / numniches, numniches),
                                  gridsquaresize), axis)
    # Create empty supply and for now fill with one value
    return _maxsupply_env(regime, active, maxsupply, area)
end
function simplenicheAE(numniches::Int64,
                       dimension::Tuple,
                       maxsupply::Unitful.Quantity{Float64},
                       area::Unitful.Area{Float64};
                       axis::Type{<:NicheAxis} = Unclassified)
    active = Matrix{Bool}(undef, dimension)
    fill!(active, true)
    return simplenicheAE(numniches, dimension, maxsupply, area, active;
                         axis = axis)
end
@doc (@doc simplenicheAE) simplenicheAE(::Int64,
                                        ::Tuple,
                                        ::Unitful.Quantity{Float64},
                                        ::Unitful.Area{Float64})

import Diversity.API: _countsubcommunities
function _countsubcommunities(gae::GridHabitat)
    return countsubcommunities(gae.regime)
end
import Diversity.API: _getsubcommunitynames
function _getsubcommunitynames(gae::GridHabitat)
    return gae.names
end

"""
    getavailablesupply(gae::GridHabitat)

Return the available resource supply from a [`GridHabitat`](@ref).
"""
function getavailablesupply(gae::GridHabitat)
    return _getavailablesupply(gae.supply)
end
"""
    tempgradAE(minT::Unitful.Temperature{Float64},
      maxT::Unitful.Temperature{Float64},
      dimension::Tuple{Int64, Int64}, maxsupply::Unitful.Quantity{Float64},
      area::Unitful.Area{Float64}, rate::Quantity{Float64, 𝚯*𝐓^-1},
      active::Matrix{Bool})

Create a temperature gradient [`ContinuousRegime`](@ref), [`SimpleSupply`](@ref)
type abiotic environment. Given a `minT` and `maxT` temperature, it generates a
gradient from minimum at the bottom to maximum at the top. It creates a
[`ContinuousRegime`](@ref) environment with dimensions `dimension` and a specified
area `area`. It also creates a [`SimpleSupply`](@ref) type filled with the
maximum supply value `maxsupply`. The rate of temperature change is specified using
the parameter `rate`. If a Bool matrix of active grid squares is included,
`active`, this is used, else one is created with all grid cells active.
"""
function tempgradAE(minT::Unitful.Temperature{Float64},
                    maxT::Unitful.Temperature{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxsupply::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                    active::Matrix{Bool};
                    axis::Type{<:NicheAxis} = MeanTemperature)
    minT = _canonical(minT, axis)
    maxT = _canonical(maxT, axis)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    regime = _reaxis(tempgrad(minT, maxT, gridsquaresize, dimension, rate),
                     axis)
    return _maxsupply_env(regime, active, maxsupply, area)
end

function tempgradAE(minT::Unitful.Temperature{Float64},
                    maxT::Unitful.Temperature{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxsupply::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝚯 * 𝐓^-1};
                    axis::Type{<:NicheAxis} = MeanTemperature)
    active = fill(true, dimension)
    return tempgradAE(minT, maxT, dimension, maxsupply, area, rate, active;
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
       dimension::Tuple{Int64, Int64}, maxsupply::Unitful.Quantity{Float64},
       area::Unitful.Area{Float64}, rate::Quantity{Float64, 𝚯*𝐓^-1},
       active::Matrix{Bool})

Create a peaked temperature gradient [`ContinuousRegime`](@ref),
[`SimpleSupply`](@ref) type abiotic environment. Given a `minT` and `maxT`
temperature, it generates a gradient with minima at the top and bottom, peaking
at `maxT` in the middle. It creates a [`ContinuousRegime`](@ref) environment with
dimensions `dimension` and a specified area `area`. It also creates a
[`SimpleSupply`](@ref) type filled with the maximum supply value `maxsupply`. The
rate of temperature change is specified using the parameter `rate`. If a Bool
matrix of active grid squares is included, `active`, this is used, else one is
created with all grid cells active.
"""
function peakedgradAE(minT::Unitful.Temperature{Float64},
                      maxT::Unitful.Temperature{Float64},
                      dimension::Tuple{Int64, Int64},
                      maxsupply::Unitful.Quantity{Float64},
                      area::Unitful.Area{Float64},
                      rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                      active::Matrix{Bool};
                      axis::Type{<:NicheAxis} = MeanTemperature)
    minT = _canonical(minT, axis)
    maxT = _canonical(maxT, axis)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    regime = tempgrad(minT, maxT + (maxT - minT), gridsquaresize, dimension,
                      rate)
    regime.matrix[(ceil(Int, dimension[1] / 2) + 1):end,
                  :] = regime.matrix[floor(Int,
                                           dimension[1] / 2):-1:1,
                                     :]
    return _maxsupply_env(_reaxis(regime, axis), active, maxsupply, area)
end

function peakedgradAE(minT::Unitful.Temperature{Float64},
                      maxT::Unitful.Temperature{Float64},
                      dimension::Tuple{Int64, Int64},
                      maxsupply::Unitful.Quantity{Float64},
                      area::Unitful.Area{Float64},
                      rate::Quantity{Float64, 𝚯 * 𝐓^-1};
                      axis::Type{<:NicheAxis} = MeanTemperature)
    active = fill(true, dimension)
    return peakedgradAE(minT, maxT, dimension, maxsupply, area, rate, active;
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
      dimension::Tuple{Int64, Int64}, maxsupply::Unitful.Quantity{Float64},
      area::Unitful.Area{Float64}, rate::Quantity{Float64, 𝐋*𝐓^-1},
      active::Matrix{Bool})

Create a rainfall gradient [`ContinuousRegime`](@ref), [`SimpleSupply`](@ref) type
abiotic environment. Given a `minR` and `maxR` rainfall, it generates a gradient
from minimum at the bottom to maximum at the top. It creates a
[`ContinuousRegime`](@ref) environment with dimensions `dimension` and a specified
area `area`. It also creates a [`SimpleSupply`](@ref) type filled with the
maximum supply value `maxsupply`. The rate of rainfall change is specified using
the parameter `rate`. If a Bool matrix of active grid squares is included,
`active`, this is used, else one is created with all grid cells active.
"""
function raingradAE(minR::Unitful.Length{Float64},
                    maxR::Unitful.Length{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxsupply::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝐋 * 𝐓^-1},
                    active::Matrix{Bool};
                    axis::Type{<:NicheAxis} = Precipitation)
    minR = _canonical(minR, axis)
    maxR = _canonical(maxR, axis)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    regime = _reaxis(raingrad(minR, maxR, gridsquaresize, dimension, rate),
                     axis)
    return _maxsupply_env(regime, active, maxsupply, area)
end

"""
    raingradAE(minR::Unitful.Length{Float64}, maxR::Unitful.Length{Float64},
      dimension::Tuple{Int64, Int64}, area::Unitful.Area{Float64},
      rate::Quantity{Float64, 𝐋*𝐓^-1}, active::Matrix{Bool})

As [`raingradAE`](@ref) but uses the rainfall values from the gradient directly
as a [`WaterSupply`](@ref), rather than computing a supply from a separate
`maxsupply` value.
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
    regime = _reaxis(raingrad(minR, maxR, gridsquaresize, dimension, rate),
                     axis)
    supply = WaterSupply(regime.matrix)
    return GridHabitat{typeof(regime), typeof(supply)}(regime, active, supply)
end

function raingradAE(minR::Unitful.Length{Float64},
                    maxR::Unitful.Length{Float64},
                    dimension::Tuple{Int64, Int64},
                    maxsupply::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    rate::Quantity{Float64, 𝐋 * 𝐓^-1};
                    axis::Type{<:NicheAxis} = Precipitation)
    active = fill(true, dimension)
    return raingradAE(minR, maxR, dimension, maxsupply, area, rate, active;
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
   eraAE(era::ERA, maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area{Float64})

Create a [`ContinuousRegime`](@ref), [`SimpleSupply`](@ref) type abiotic
environment from an ERA type climate. It either creates a [`SimpleSupply`](@ref)
type filled with the maximum supply value `maxsupply` or uses a provided supply of
type [`SolarTimeSupply`](@ref). If a Bool matrix of active grid squares is
included, `active`, this is used, else one is created with all grid cells
active.
"""
function eraAE(era::ERA, maxsupply::Unitful.Quantity{Float64},
               area::Unitful.Area{Float64})
    return _maxsupply_env(_timeclimate_regime(era.array), _nanactive(era.array),
                          maxsupply, area)
end
"""
    eraAE(era::ERA, maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area{Float64},
      active::Matrix{Bool})

As [`eraAE`](@ref) with an explicit `active` matrix of grid squares, rather than
inferring active cells from NaN values in the ERA data.
"""
function eraAE(era::ERA,
               maxsupply::Unitful.Quantity{Float64},
               area::Unitful.Area{Float64},
               active::Matrix{Bool})
    return _maxsupply_env(_timeclimate_regime(era.array), active, maxsupply,
                          area)
end
"""
    eraAE(era::ERA, supply::B, active::Matrix{Bool}) where B <: AbstractTimeSupply

As [`eraAE`](@ref) but accepts a pre-constructed `AbstractTimeSupply` object
`supply` rather than computing a supply from a maximum value.
"""
function eraAE(era::ERA, supply::B,
               active::Matrix{Bool}) where {B <: AbstractTimeSupply}
    regime = _timeclimate_regime(era.array)
    return GridHabitat{typeof(regime), typeof(supply)}(regime, active, supply)
end

"""
   worldclimAE(wc::ClimateRaster{WorldClim{Climate}}, maxsupply::Unitful.Quantity{Float64},
     area::Unitful.Area{Float64})

Create a [`ContinuousTimeRegime`](@ref), [`SimpleSupply`](@ref) type abiotic
environment from a Worldclim type climate. It either creates a
[`SimpleSupply`](@ref) type filled with the maximum supply value `maxsupply` or
uses a provided supply of type [`SolarTimeSupply`](@ref). If a Bool matrix of
active grid squares is included, `active`, this is used, otherwise all grid
cells are considered active.
"""
function worldclimAE(wc::ClimateRaster{WorldClim{Climate}},
                     maxsupply::Unitful.Quantity{Float64},
                     area::Unitful.Area{Float64})
    return _maxsupply_env(_timeclimate_regime(wc.array), _nanactive(wc.array),
                          maxsupply, area)
end
"""
    worldclimAE(wc::ClimateRaster{WorldClim{Climate}}, maxsupply::Unitful.Quantity{Float64},
      area::Unitful.Area{Float64}, active::Matrix{Bool})

As [`worldclimAE`](@ref) with an explicit `active` matrix of grid squares,
rather than inferring active cells from NaN values in the Worldclim data.
"""
function worldclimAE(wc::ClimateRaster{WorldClim{Climate}},
                     maxsupply::Unitful.Quantity{Float64},
                     area::Unitful.Area{Float64},
                     active::Matrix{Bool})
    return _maxsupply_env(_timeclimate_regime(wc.array), active, maxsupply,
                          area)
end
"""
    worldclimAE(wc::ClimateRaster{WorldClim{Climate}}, supply::B, active::Matrix{Bool}) where B <: AbstractTimeSupply

As [`worldclimAE`](@ref) but accepts a pre-constructed `AbstractTimeSupply`
object `supply` rather than computing a supply from a maximum value.
"""
function worldclimAE(wc::ClimateRaster{WorldClim{Climate}},
                     supply::B,
                     active::Matrix{Bool}) where {B <: AbstractTimeSupply}
    regime = _timeclimate_regime(wc.array)
    return GridHabitat{typeof(regime), typeof(supply)}(regime, active, supply)
end

"""
    bioclimAE(bc::ClimateRaster{WorldClim{BioClim}}, maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area{Float64})

Create a [`ContinuousRegime`](@ref), [`SimpleSupply`](@ref) type abiotic
environment from a Worldclim type climate. It either creates a
[`SimpleSupply`](@ref) type filled with the maximum supply value `maxsupply` or
uses a provided supply of type [`SolarSupply`](@ref). If a Bool matrix of active
grid squares is included, `active`, this is used, else one is created with all
grid cells active.
"""
function bioclimAE(bc::ClimateRaster{WorldClim{BioClim}},
                   maxsupply::Unitful.Quantity{Float64},
                   area::Unitful.Area{Float64})
    return _maxsupply_env(_continuousregime(bc.array), _nanactive(bc.array),
                          maxsupply, area)
end
"""
    bioclimAE(bc::ClimateRaster{WorldClim{BioClim}}, maxsupply::Unitful.Quantity{Float64},
      area::Unitful.Area{Float64}, active::Matrix{Bool})

As [`bioclimAE`](@ref) with an explicit `active` matrix of grid squares, rather
than inferring active cells from NaN values in the bioclim data.
"""
function bioclimAE(bc::ClimateRaster{WorldClim{BioClim}},
                   maxsupply::Unitful.Quantity{Float64},
                   area::Unitful.Area{Float64}, active::Matrix{Bool})
    return _maxsupply_env(_continuousregime(bc.array), active, maxsupply, area)
end
"""
    bioclimAE(bc::ClimateRaster{WorldClim{BioClim}}, supply::B, active::Matrix{Bool}) where B <: AbstractSupply

As [`bioclimAE`](@ref) but accepts a pre-constructed [`AbstractSupply`](@ref)
object `supply` rather than computing a supply from a maximum value.
"""
function bioclimAE(bc::ClimateRaster{WorldClim{BioClim}}, supply::B,
                   active::Matrix{Bool}) where {B <: AbstractSupply}
    regime = _continuousregime(bc.array)
    return GridHabitat{typeof(regime), typeof(supply)}(regime, active, supply)
end

"""
    simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
        dimension::Tuple{Int64, Int64}, maxsupply::Float64, area::Unitful.Area{Float64},
        active::Matrix{Bool})

Create a simple [`ContinuousRegime`](@ref), [`SimpleSupply`](@ref) type abiotic
environment. It creates a [`ContinuousRegime`](@ref) filled with a given value,
`val`, dimensions (`dimension`) and a specified area (`area`). It also creates a
[`SimpleSupply`](@ref) type filled with the maximum supply value (`maxsupply`). If
a Bool matrix of active grid squares is included, `active`, this is used, else
one is created with all grid cells active.
"""
function simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
                         dimension::Tuple{Int64, Int64},
                         maxsupply::Unitful.Quantity{Float64},
                         area::Unitful.Area{Float64},
                         active::Matrix{Bool};
                         axis::Type{<:NicheAxis} = Unclassified)
    val = _canonical(val, axis)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    regime = _reaxis(simpleregime(val, gridsquaresize, dimension, axis), axis)
    return _maxsupply_env(regime, active, maxsupply, area)
end

function simplehabitatAE(val::Union{Float64, Unitful.Quantity{Float64}},
                         dimension::Tuple{Int64, Int64},
                         maxsupply::Unitful.Quantity{Float64},
                         area::Unitful.Area{Float64};
                         axis::Type{<:NicheAxis} = Unclassified)
    active = fill(true, dimension)
    return simplehabitatAE(val, dimension, maxsupply, area, active; axis = axis)
end
@doc (@doc simplehabitatAE) simplehabitatAE(::Union{Float64,
                                                    Unitful.Quantity{Float64}},
                                            ::Tuple{Int64, Int64},
                                            ::Unitful.Quantity{Float64},
                                            ::Unitful.Area{Float64})

import EcoBase.getcoords

getcoords(habitat::GridHabitat) = habitat.regime

"""
    lcAE(lc::ClimateRaster{<:EarthEnv{<:LandCover}}, maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area)

Create a [`DiscreteRegime`](@ref), [`SimpleSupply`](@ref) type abiotic environment
from a land cover [`ClimateRaster`](@ref) dataset. It creates a
[`DiscreteRegime`](@ref) regime from the land cover array and a
[`SimpleSupply`](@ref) type filled with the maximum supply value `maxsupply`, scaled
to the given `area`. If a Bool matrix of active grid squares is included,
`active`, this is used, else one is created with all grid cells active.
"""
function lcAE(lc::ClimateRaster{T, A}, maxsupply::Unitful.Quantity{Float64},
              area::Unitful.Area) where {T <: EarthEnv{<:LandCover}, A}
    return _maxsupply_env(_discreteregime(lc.array), _nanactive(lc.array),
                          maxsupply, area)
end
"""
    lcAE(lc::ClimateRaster{<:EarthEnv{<:LandCover}}, maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area,
      active::Matrix{Bool})

As [`lcAE`](@ref) with an explicit `active` matrix of grid squares, rather than
setting all cells active.
"""
function lcAE(lc::ClimateRaster{T, A}, maxsupply::Unitful.Quantity{Float64},
              area::Unitful.Area,
              active::Matrix{Bool}) where {T <: EarthEnv{<:LandCover}, A}
    return _maxsupply_env(_discreteregime(lc.array), active, maxsupply, area)
end
"""
    lcAE(lc::ClimateRaster{<:EarthEnv{<:LandCover}}, supply::B, active::Matrix{Bool}) where B <: AbstractSupply

As [`lcAE`](@ref) but accepts a pre-constructed [`AbstractSupply`](@ref) object
`supply` rather than computing a supply from a maximum value.
"""
function lcAE(lc::ClimateRaster{T, A}, supply::B,
              active::Matrix{Bool}) where {T <: EarthEnv{<:LandCover},
                                           A, B <: AbstractSupply}
    regime = _discreteregime(lc.array)
    return GridHabitat{typeof(regime), B}(regime, active, supply)
end
