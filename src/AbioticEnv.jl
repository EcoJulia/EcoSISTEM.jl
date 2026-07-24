# SPDX-License-Identifier: LGPL-3.0-or-later

using Diversity
using Unitful
using Unitful.DefaultSymbols
using Diversity.API
using RasterDataSources
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref

const RDS = RasterDataSources

# Turn a per-area *rate* (an areal flux — energy/water per unit area per unit time) into an
# absolute Resource quantity (per cell, not per area) by multiplying by that cell's `area` —
# dispatched on the input's own dimension, not a `Dict`. Solar's `kJ/m²/day` (𝐌𝐓⁻³) × m²
# (𝐋²) → `kJ/day` (𝐋²𝐌𝐓⁻³, `Unitful.Power`); water's `L/m²/day` (𝐋𝐓⁻¹) × m² → `L/day`
# (𝐋³𝐓⁻¹, `VolumeFlow`); the free/simple case is dimensionless-per-area (𝐋⁻²𝐓⁻¹) × m² →
# a bare rate (𝐓⁻¹). All three verified directly against `dimension(...)` before wiring.
function cancel(a::Quantity{<:Real, 𝐌 * 𝐓^-3}, b::Quantity{<:Real, 𝐋^2})
    return uconvert(unit(_SolarRate), a * b)
end
function cancel(a::Quantity{<:Real, 𝐋 * 𝐓^-1}, b::Quantity{<:Real, 𝐋^2})
    return uconvert(unit(_WaterRate), a * b)
end
function cancel(a::Quantity{<:Real, 𝐋^-2 * 𝐓^-1}, b::Quantity{<:Real, 𝐋^2})
    return uconvert(unit(_SimpleRate), a * b)
end

# Shared tail of the maximum-supply `*AE` constructors: total the per-area maximum supply
# `maxsupply` (an areal rate) over `area`, spread it uniformly across the grid of `active`,
# pick the supply type from its dimension (`supplytype`, dispatched — see `NicheInfo.jl`),
# and assemble the `GridHabitat` with regime `regime`.
function _maxsupply_env(regime, active::Matrix{Bool},
                        maxsupply::Unitful.Quantity{Float64},
                        area::Unitful.Area)
    B = cancel(maxsupply, area)
    T = supplytype(typeof(B))
    supply = fill(B / length(active), size(active))
    return GridHabitat{typeof(regime), T}(regime, active, T(supply))
end

# The time-varying (monthly) climate regime shared by `erahabitat`/`worldclimhabitat`: a
# `ContinuousTimeRegime` starting at month 1 that advances via `cyclicChange`, tagged with niche `axis`.
function _timeclimate_regime(array, axis::Type{<:NicheAxis} = Unclassified)
    return _reaxis(ContinuousTimeRegime(Array(array), 1, _cellsize(array),
                                        LayerUpdate(cyclicChange, 0.0 / s,
                                                    Unitful.Dimensions{()})),
                   axis)
end

# Static (no-dynamics) data regimes shared by `bioclimhabitat` (continuous) and `landcoverhabitat`
# (discrete land cover), each tagged with niche `axis`.
function _continuousregime(array, axis::Type{<:NicheAxis} = Unclassified)
    return _reaxis(ContinuousRegime(Array(array), _cellsize(array),
                                    LayerUpdate(NoChange, 0.0 / s,
                                                Unitful.Dimensions{()})), axis)
end
function _discreteregime(array, axis::Type{<:NicheAxis} = Unclassified)
    return _reaxis(DiscreteRegime(Array(array), _cellsize(array),
                                  LayerUpdate(NoChange, 0.0 / s)), axis)
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
    simplenichehabitat(numniches::Int64, dimension::Tuple,
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
function simplenichehabitat(numniches::Int64,
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
function simplenichehabitat(numniches::Int64,
                            dimension::Tuple,
                            maxsupply::Unitful.Quantity{Float64},
                            area::Unitful.Area{Float64};
                            axis::Type{<:NicheAxis} = Unclassified)
    active = Matrix{Bool}(undef, dimension)
    fill!(active, true)
    return simplenichehabitat(numniches, dimension, maxsupply, area, active;
                              axis = axis)
end
@doc (@doc simplenichehabitat) simplenichehabitat(::Int64,
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
    tempgradhabitat(minT::Unitful.Temperature{Float64},
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
function tempgradhabitat(minT::Unitful.Temperature{Float64},
                         maxT::Unitful.Temperature{Float64},
                         dimension::Tuple{Int64, Int64},
                         maxsupply::Unitful.Quantity{Float64},
                         area::Unitful.Area{Float64},
                         rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                         active::Matrix{Bool};
                         axis::Type{<:NicheAxis} = Temperature)
    minT = _canonical(minT, axis)
    maxT = _canonical(maxT, axis)
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))
    regime = _reaxis(tempgrad(minT, maxT, gridsquaresize, dimension, rate),
                     axis)
    return _maxsupply_env(regime, active, maxsupply, area)
end

function tempgradhabitat(minT::Unitful.Temperature{Float64},
                         maxT::Unitful.Temperature{Float64},
                         dimension::Tuple{Int64, Int64},
                         maxsupply::Unitful.Quantity{Float64},
                         area::Unitful.Area{Float64},
                         rate::Quantity{Float64, 𝚯 * 𝐓^-1};
                         axis::Type{<:NicheAxis} = Temperature)
    active = fill(true, dimension)
    return tempgradhabitat(minT, maxT, dimension, maxsupply, area, rate, active;
                           axis = axis)
end
@doc (@doc tempgradhabitat) tempgradhabitat(::Unitful.Temperature{Float64},
                                            ::Unitful.Temperature{Float64},
                                            ::Tuple{Int64, Int64},
                                            ::Unitful.Quantity{Float64},
                                            ::Unitful.Area{Float64},
                                            ::Quantity{Float64, 𝚯 * 𝐓^-1})

"""
    peakedgradhabitat(minT::Unitful.Temperature{Float64},
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
function peakedgradhabitat(minT::Unitful.Temperature{Float64},
                           maxT::Unitful.Temperature{Float64},
                           dimension::Tuple{Int64, Int64},
                           maxsupply::Unitful.Quantity{Float64},
                           area::Unitful.Area{Float64},
                           rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                           active::Matrix{Bool};
                           axis::Type{<:NicheAxis} = Temperature)
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

function peakedgradhabitat(minT::Unitful.Temperature{Float64},
                           maxT::Unitful.Temperature{Float64},
                           dimension::Tuple{Int64, Int64},
                           maxsupply::Unitful.Quantity{Float64},
                           area::Unitful.Area{Float64},
                           rate::Quantity{Float64, 𝚯 * 𝐓^-1};
                           axis::Type{<:NicheAxis} = Temperature)
    active = fill(true, dimension)
    return peakedgradhabitat(minT, maxT, dimension, maxsupply, area, rate,
                             active;
                             axis = axis)
end
@doc (@doc peakedgradhabitat) peakedgradhabitat(::Unitful.Temperature{Float64},
                                                ::Unitful.Temperature{Float64},
                                                ::Tuple{Int64, Int64},
                                                ::Unitful.Quantity{Float64},
                                                ::Unitful.Area{Float64},
                                                ::Quantity{Float64, 𝚯 * 𝐓^-1})

"""
    raingradhabitat(minR::Unitful.Length{Float64},
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
function raingradhabitat(minR::Unitful.Length{Float64},
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
    raingradhabitat(minR::Unitful.Length{Float64}, maxR::Unitful.Length{Float64},
      dimension::Tuple{Int64, Int64}, area::Unitful.Area{Float64},
      rate::Quantity{Float64, 𝐋*𝐓^-1}, active::Matrix{Bool})

As [`raingradhabitat`](@ref) but uses the rainfall values from the gradient directly
as a [`WaterSupply`](@ref), rather than computing a supply from a separate
`maxsupply` value.
"""
function raingradhabitat(minR::Unitful.Length{Float64},
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

function raingradhabitat(minR::Unitful.Length{Float64},
                         maxR::Unitful.Length{Float64},
                         dimension::Tuple{Int64, Int64},
                         maxsupply::Unitful.Quantity{Float64},
                         area::Unitful.Area{Float64},
                         rate::Quantity{Float64, 𝐋 * 𝐓^-1};
                         axis::Type{<:NicheAxis} = Precipitation)
    active = fill(true, dimension)
    return raingradhabitat(minR, maxR, dimension, maxsupply, area, rate, active;
                           axis = axis)
end
@doc (@doc raingradhabitat) raingradhabitat(::Unitful.Length{Float64},
                                            ::Unitful.Length{Float64},
                                            ::Tuple{Int64, Int64},
                                            ::Unitful.Quantity{Float64},
                                            ::Unitful.Area{Float64},
                                            ::Quantity{Float64, 𝐋 * 𝐓^-1})
function raingradhabitat(minR::Unitful.Length{Float64},
                         maxR::Unitful.Length{Float64},
                         dimension::Tuple{Int64, Int64},
                         area::Unitful.Area{Float64},
                         rate::Quantity{Float64, 𝐋 * 𝐓^-1};
                         axis::Type{<:NicheAxis} = Precipitation)
    active = fill(true, dimension)
    return raingradhabitat(minR, maxR, dimension, area, rate, active;
                           axis = axis)
end
@doc (@doc raingradhabitat) raingradhabitat(::Unitful.Length{Float64},
                                            ::Unitful.Length{Float64},
                                            ::Tuple{Int64, Int64},
                                            ::Unitful.Area{Float64},
                                            ::Quantity{Float64, 𝐋 * 𝐓^-1})
"""
   erahabitat(era::ERA, maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area{Float64})

Create a [`ContinuousRegime`](@ref), [`SimpleSupply`](@ref) type abiotic
environment from an ERA type climate. It either creates a [`SimpleSupply`](@ref)
type filled with the maximum supply value `maxsupply` or uses a provided supply of
type [`SolarTimeSupply`](@ref). If a Bool matrix of active grid squares is
included, `active`, this is used, else one is created with all grid cells
active.
"""
function erahabitat(era::ERA, maxsupply::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64};
                    axis::Type{<:NicheAxis} = Unclassified)
    return _maxsupply_env(_timeclimate_regime(era.array, axis),
                          _nanactive(era.array), maxsupply, area)
end
"""
    erahabitat(era::ERA, maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area{Float64},
      active::Matrix{Bool})

As [`erahabitat`](@ref) with an explicit `active` matrix of grid squares, rather than
inferring active cells from NaN values in the ERA data.
"""
function erahabitat(era::ERA,
                    maxsupply::Unitful.Quantity{Float64},
                    area::Unitful.Area{Float64},
                    active::Matrix{Bool};
                    axis::Type{<:NicheAxis} = Unclassified)
    return _maxsupply_env(_timeclimate_regime(era.array, axis), active,
                          maxsupply, area)
end
"""
    erahabitat(era::ERA, supply::B, active::Matrix{Bool}) where B <: AbstractTimeSupply

As [`erahabitat`](@ref) but accepts a pre-constructed `AbstractTimeSupply` object
`supply` rather than computing a supply from a maximum value.
"""
function erahabitat(era::ERA, supply::B,
                    active::Matrix{Bool};
                    axis::Type{<:NicheAxis} = Unclassified) where {B <:
                                                                   AbstractTimeSupply}
    regime = _timeclimate_regime(era.array, axis)
    return GridHabitat{typeof(regime), typeof(supply)}(regime, active, supply)
end

"""
   worldclimhabitat(worldclim::ClimateRaster{WorldClim{Climate}}, maxsupply::Unitful.Quantity{Float64},
     area::Unitful.Area{Float64})

Create a [`ContinuousTimeRegime`](@ref), [`SimpleSupply`](@ref) type abiotic
environment from a Worldclim type climate. It either creates a
[`SimpleSupply`](@ref) type filled with the maximum supply value `maxsupply` or
uses a provided supply of type [`SolarTimeSupply`](@ref). If a Bool matrix of
active grid squares is included, `active`, this is used, otherwise all grid
cells are considered active.
"""
function worldclimhabitat(worldclim::ClimateRaster{WorldClim{Climate}},
                          maxsupply::Unitful.Quantity{Float64},
                          area::Unitful.Area{Float64};
                          axis::Type{<:NicheAxis} = Unclassified)
    return _maxsupply_env(_timeclimate_regime(worldclim.array, axis),
                          _nanactive(worldclim.array), maxsupply, area)
end
"""
    worldclimhabitat(worldclim::ClimateRaster{WorldClim{Climate}}, maxsupply::Unitful.Quantity{Float64},
      area::Unitful.Area{Float64}, active::Matrix{Bool})

As [`worldclimhabitat`](@ref) with an explicit `active` matrix of grid squares,
rather than inferring active cells from NaN values in the Worldclim data.
"""
function worldclimhabitat(worldclim::ClimateRaster{WorldClim{Climate}},
                          maxsupply::Unitful.Quantity{Float64},
                          area::Unitful.Area{Float64},
                          active::Matrix{Bool};
                          axis::Type{<:NicheAxis} = Unclassified)
    return _maxsupply_env(_timeclimate_regime(worldclim.array, axis), active,
                          maxsupply, area)
end
"""
    worldclimhabitat(worldclim::ClimateRaster{WorldClim{Climate}}, supply::B, active::Matrix{Bool}) where B <: AbstractTimeSupply

As [`worldclimhabitat`](@ref) but accepts a pre-constructed `AbstractTimeSupply`
object `supply` rather than computing a supply from a maximum value.
"""
function worldclimhabitat(worldclim::ClimateRaster{WorldClim{Climate}},
                          supply::B,
                          active::Matrix{Bool};
                          axis::Type{<:NicheAxis} = Unclassified) where {B <:
                                                                         AbstractTimeSupply}
    regime = _timeclimate_regime(worldclim.array, axis)
    return GridHabitat{typeof(regime), typeof(supply)}(regime, active, supply)
end

"""
    bioclimhabitat(bioclim::ClimateRaster{WorldClim{BioClim}}, maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area{Float64})

Create a [`ContinuousRegime`](@ref), [`SimpleSupply`](@ref) type abiotic
environment from a Worldclim type climate. It either creates a
[`SimpleSupply`](@ref) type filled with the maximum supply value `maxsupply` or
uses a provided supply of type [`SolarSupply`](@ref). If a Bool matrix of active
grid squares is included, `active`, this is used, else one is created with all
grid cells active.
"""
function bioclimhabitat(bioclim::ClimateRaster{WorldClim{BioClim}},
                        maxsupply::Unitful.Quantity{Float64},
                        area::Unitful.Area{Float64};
                        axis::Type{<:NicheAxis} = Unclassified)
    return _maxsupply_env(_continuousregime(bioclim.array, axis),
                          _nanactive(bioclim.array), maxsupply, area)
end
"""
    bioclimhabitat(bioclim::ClimateRaster{WorldClim{BioClim}}, maxsupply::Unitful.Quantity{Float64},
      area::Unitful.Area{Float64}, active::Matrix{Bool})

As [`bioclimhabitat`](@ref) with an explicit `active` matrix of grid squares, rather
than inferring active cells from NaN values in the bioclim data.
"""
function bioclimhabitat(bioclim::ClimateRaster{WorldClim{BioClim}},
                        maxsupply::Unitful.Quantity{Float64},
                        area::Unitful.Area{Float64}, active::Matrix{Bool};
                        axis::Type{<:NicheAxis} = Unclassified)
    return _maxsupply_env(_continuousregime(bioclim.array, axis), active,
                          maxsupply,
                          area)
end
"""
    bioclimhabitat(bioclim::ClimateRaster{WorldClim{BioClim}}, supply::B, active::Matrix{Bool}) where B <: AbstractSupply

As [`bioclimhabitat`](@ref) but accepts a pre-constructed [`AbstractSupply`](@ref)
object `supply` rather than computing a supply from a maximum value.
"""
function bioclimhabitat(bioclim::ClimateRaster{WorldClim{BioClim}}, supply::B,
                        active::Matrix{Bool};
                        axis::Type{<:NicheAxis} = Unclassified) where {B <:
                                                                       AbstractSupply}
    regime = _continuousregime(bioclim.array, axis)
    return GridHabitat{typeof(regime), typeof(supply)}(regime, active, supply)
end

"""
    simplehabitat(val::Union{Float64, Unitful.Quantity{Float64}},
        dimension::Tuple{Int64, Int64}, maxsupply::Float64, area::Unitful.Area{Float64},
        active::Matrix{Bool})

Create a simple [`ContinuousRegime`](@ref), [`SimpleSupply`](@ref) type abiotic
environment. It creates a [`ContinuousRegime`](@ref) filled with a given value,
`val`, dimensions (`dimension`) and a specified area (`area`). It also creates a
[`SimpleSupply`](@ref) type filled with the maximum supply value (`maxsupply`). If
a Bool matrix of active grid squares is included, `active`, this is used, else
one is created with all grid cells active.
"""
function simplehabitat(val::Union{Float64, Unitful.Quantity{Float64}},
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

function simplehabitat(val::Union{Float64, Unitful.Quantity{Float64}},
                       dimension::Tuple{Int64, Int64},
                       maxsupply::Unitful.Quantity{Float64},
                       area::Unitful.Area{Float64};
                       axis::Type{<:NicheAxis} = Unclassified)
    active = fill(true, dimension)
    return simplehabitat(val, dimension, maxsupply, area, active; axis = axis)
end
@doc (@doc simplehabitat) simplehabitat(::Union{Float64,
                                                Unitful.Quantity{Float64}},
                                        ::Tuple{Int64, Int64},
                                        ::Unitful.Quantity{Float64},
                                        ::Unitful.Area{Float64})

import EcoBase.getcoords

getcoords(habitat::GridHabitat) = habitat.regime

"""
    landcoverhabitat(landcover::ClimateRaster{<:EarthEnv{<:LandCover}}, maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area)

Create a [`DiscreteRegime`](@ref), [`SimpleSupply`](@ref) type abiotic environment
from a land cover [`ClimateRaster`](@ref) dataset. It creates a
[`DiscreteRegime`](@ref) regime from the land cover array and a
[`SimpleSupply`](@ref) type filled with the maximum supply value `maxsupply`, scaled
to the given `area`. If a Bool matrix of active grid squares is included,
`active`, this is used, else one is created with all grid cells active.
"""
function landcoverhabitat(landcover::ClimateRaster{T, A},
                          maxsupply::Unitful.Quantity{Float64},
                          area::Unitful.Area;
                          axis::Type{<:NicheAxis} = Unclassified) where {T <:
                                                                         EarthEnv{<:LandCover},
                                                                         A}
    return _maxsupply_env(_discreteregime(landcover.array, axis),
                          _nanactive(landcover.array), maxsupply, area)
end
"""
    landcoverhabitat(landcover::ClimateRaster{<:EarthEnv{<:LandCover}}, maxsupply::Unitful.Quantity{Float64}, area::Unitful.Area,
      active::Matrix{Bool})

As [`landcoverhabitat`](@ref) with an explicit `active` matrix of grid squares, rather than
setting all cells active.
"""
function landcoverhabitat(landcover::ClimateRaster{T, A},
                          maxsupply::Unitful.Quantity{Float64},
                          area::Unitful.Area,
                          active::Matrix{Bool};
                          axis::Type{<:NicheAxis} = Unclassified) where {T <:
                                                                         EarthEnv{<:LandCover},
                                                                         A}
    return _maxsupply_env(_discreteregime(landcover.array, axis), active,
                          maxsupply, area)
end
"""
    landcoverhabitat(landcover::ClimateRaster{<:EarthEnv{<:LandCover}}, supply::B, active::Matrix{Bool}) where B <: AbstractSupply

As [`landcoverhabitat`](@ref) but accepts a pre-constructed [`AbstractSupply`](@ref) object
`supply` rather than computing a supply from a maximum value.
"""
function landcoverhabitat(landcover::ClimateRaster{T, A}, supply::B,
                          active::Matrix{Bool};
                          axis::Type{<:NicheAxis} = Unclassified) where {T <:
                                                                         EarthEnv{<:LandCover},
                                                                         A,
                                                                         B <:
                                                                         AbstractSupply}
    regime = _discreteregime(landcover.array, axis)
    return GridHabitat{typeof(regime), B}(regime, active, supply)
end
