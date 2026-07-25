# SPDX-License-Identifier: LGPL-3.0-or-later

using EcoSISTEM.ClimatePref
using RecipesBase
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources

import Base: eltype, length
"""
    Abstract1Demand{R}

Abstract supertype for all species resource demand types, parameterised by
the type(s) of resource required `R`.
"""
abstract type AbstractDemand{R} end
abstract type Abstract1Demand{R} <: AbstractDemand{R} end
abstract type Abstract2Demands{R} <: AbstractDemand{R} end

numdemands(::Type{<:Abstract1Demand}) = 1
numdemands(::Type{<:Abstract2Demands}) = 2

function Base.eltype(::Abstract1Demand{R}) where {R}
    return R
end

"""
    SimpleDemand <: Abstract1Demand{typeof(1.0/day)}

A simple resource demand — a per-species free-resource rate for each species, `/day`. Accepts
a plain `Vector{<:Real}` (auto-tagged `/day`) so the "no units hassle" case stays easy.
"""
struct SimpleDemand <: Abstract1Demand{_SimpleRate}
    resource::Vector{_SimpleRate}
    exchange_rate::typeof(1.0 * day)

    function SimpleDemand(resource::Vector{<:Real},
                          exchange_rate::Unitful.Quantity{Float64} = 1.0 * day)
        return new(resource .* unit(_SimpleRate), uconvert(day, exchange_rate))
    end
end

Base.length(demand::SimpleDemand) = length(demand.resource)

function _getdemand(abun::Vector{Int64}, demand::SimpleDemand)
    return sum(abun .* demand.resource)
end

"""
    SizeDemand <: Abstract1Demand{typeof(1.0/day)}

A resource demand scaled by species size — a per-species free-resource rate for each
species, `/day` (see [`SimpleDemand`](@ref); `pop_mass_rel`/`area` are used only once, at
initial-abundance construction time in `SpeciesList`, not by the ongoing demand:supply
balance).
"""
struct SizeDemand <: Abstract1Demand{_SimpleRate}
    resource::Vector{_SimpleRate}
    pop_mass_rel::Float64
    area::Unitful.Area
    exchange_rate::typeof(1.0 * day)

    function SizeDemand(resource::Vector{<:Real},
                        pop_mass_rel::Float64,
                        area::Unitful.Area,
                        exchange_rate::Unitful.Quantity{Float64} = 1.0 * day)
        return new(resource .* unit(_SimpleRate), pop_mass_rel, area,
                   uconvert(day, exchange_rate))
    end
end

Base.length(demand::SizeDemand) = length(demand.resource)
function _getdemand(abun::Vector{Int64}, demand::SizeDemand)
    return sum(abun .* demand.resource)
end

"""
    SolarDemand <: Abstract1Demand{typeof(1.0*kJ/day)}

A vector of per-species solar resource demands, `kJ/day`.
"""
struct SolarDemand <: Abstract1Demand{_SolarRate}
    resource::Vector{_SolarRate}
    exchange_rate::typeof(inv(oneunit(_SolarRate)))

    function SolarDemand(resource::Vector{<:Unitful.Power{Float64}},
                         exchange_rate::Unitful.Quantity{Float64} = 1.0 /
                                                                    mean(resource))
        return new(uconvert.(unit(_SolarRate), resource),
                   uconvert(inv(unit(_SolarRate)), exchange_rate))
    end
end

Base.length(demand::SolarDemand) = length(demand.resource)
function _getdemand(abun::Vector{Int64}, demand::SolarDemand)
    return sum(abun .* demand.resource)
end

"""
    WaterDemand <: Abstract1Demand{typeof(1.0*L/day)}

A vector of per-species water demands, `L/day`.
"""
struct WaterDemand <: Abstract1Demand{_WaterRate}
    resource::Vector{_WaterRate}
    exchange_rate::typeof(inv(oneunit(_WaterRate)))

    function WaterDemand(resource::Vector{<:VolumeFlow{Float64}},
                         exchange_rate::Unitful.Quantity{Float64} = 1.0 /
                                                                    mean(resource))
        return new(uconvert.(unit(_WaterRate), resource),
                   uconvert(inv(unit(_WaterRate)), exchange_rate))
    end
end
Base.length(demand::WaterDemand) = length(demand.resource)
function _getdemand(abun::Vector{Int64}, demand::WaterDemand)
    return sum(abun .* demand.resource)
end

"""
    DemandCollection2{R1, R2}

A pair of species resource demands (e.g. solar resource and water) consumed
together, to match a two-resource supply.
"""
struct DemandCollection2{R1, R2} <: Abstract2Demands{Tuple{R1, R2}}
    one::R1
    two::R2
end
Base.length(demand::DemandCollection2) = length(demand.one.resource)
function Base.eltype(demand::DemandCollection2)
    return [eltype(demand.one), eltype(demand.two)]
end
function _getdemand(abun::Vector{Int64}, demand::DemandCollection2)
    return [_getdemand(abun, demand.one), _getdemand(abun, demand.two)]
end

# A plot title for a Resource layer's element type, dispatched on its dimension (mirrors
# `supplytype`/`demandtype` in `NicheInfo.jl`) rather than looked up by unit in a `Dict`.
_resourcetitle(::Type{<:Unitful.Power}) = "Solar Radiation (kJ/day)"
_resourcetitle(::Type{<:VolumeFlow}) = "Available water (L/day)"
_resourcetitle(::Type{<:Real}) = "Free resource"

# ---------------------------------------------------------------------------
# Supplies — `Resource`-role layers (folded onto `ContinuousLayer{Resource}`)
# ---------------------------------------------------------------------------
# The concrete supply names (`SimpleSupply`/`SolarSupply`/… and the time-varying
# `*TimeSupply`s) are aliases (in Layer.jl) over `ContinuousLayer{Resource, axis, V, Arr}`:
# a static supply is 2-D (`Matrix`), a time supply 3-D (`Array`, indexed by `time`). The
# constructors below fill the (unused) `size` and the per-timestep `dynamics` rule and zero
# NaNs, reproducing the old supply structs. A supply's `size` is never read
# (geometry/dispersal use the regime), so a placeholder is stored; `NoChange`/`cyclicChange`
# live in HabitatUpdate.jl (included later) and resolve at call time.
const _SUPPLY_SIZE = 1.0m
_supply_static() = LayerUpdate(NoChange, 0.0 / s, Unitful.Dimensions{()})
_supply_cyclic() = LayerUpdate(cyclicChange, 0.0 / s, Unitful.Dimensions{()})

Base.eltype(::ContinuousLayer{Resource, A, V}) where {A, V} = V
function Base.eltype(supply::LayerCollection2{Resource})
    return [eltype(supply.one), eltype(supply.two)]
end

countsubcommunities(ab::AbstractSupply) = _countsubcommunities(ab)
function _countsubcommunities(supply::ContinuousLayer{Resource, A, V, Arr}) where {A,
                                                                                   V,
                                                                                   Arr <:
                                                                                   AbstractMatrix{V}}
    return length(supply.matrix)
end
function _countsubcommunities(supply::ContinuousLayer{Resource, A, V, Arr}) where {A,
                                                                                   V,
                                                                                   Arr <:
                                                                                   AbstractArray{V,
                                                                                                 3}}
    return length(@view supply.matrix[:, :, 1])
end
function _countsubcommunities(supply::LayerCollection2{Resource})
    return _countsubcommunities(supply.one)
end

# The resource available in each cell: the full matrix (static), or the current time slice.
function _getsupply(supply::ContinuousLayer{Resource, A, V, Arr}) where {A, V,
                                                                         Arr <:
                                                                         AbstractMatrix{V}}
    return supply.matrix
end
function _getsupply(supply::ContinuousLayer{Resource, A, V, Arr}) where {A, V,
                                                                         Arr <:
                                                                         AbstractArray{V,
                                                                                       3}}
    return @view supply.matrix[:, :, supply.time]
end
function _getsupply(supply::LayerCollection2{Resource}, field::Symbol)
    return _getsupply(getfield(supply, field))
end

function _getavailablesupply(supply::ContinuousLayer{Resource})
    return sum(supply.matrix[.!isnan.(supply.matrix)])
end
function _getavailablesupply(supply::LayerCollection2{Resource})
    return [_getavailablesupply(supply.one), _getavailablesupply(supply.two)]
end

# --- Constructors reproducing the old per-type supply structs -------------------------
# Accepts a plain `Matrix{<:Real}` and auto-tags with `/day`, matching `SimpleDemand`'s same
# "no units hassle" ergonomics.
function SimpleSupply(mat::Matrix{<:Real})
    return SimpleSupply(mat .* unit(_SimpleRate))
end
function SimpleSupply(mat::Matrix{_SimpleRate})
    return ContinuousLayer{Resource, Unclassified, _SimpleRate,
                           Matrix{_SimpleRate}}(mat, 1, _SUPPLY_SIZE,
                                                _supply_static())
end
function SolarSupply(mat::Matrix{_SolarRate})
    mat[isnan.(mat)] .= zero(_SolarRate)
    return ContinuousLayer{Resource, SolarRadiation, _SolarRate,
                           Matrix{_SolarRate}}(mat, 1, _SUPPLY_SIZE,
                                               _supply_static())
end
function WaterSupply(mat::Matrix{_WaterRate})
    mat[isnan.(mat)] .= zero(_WaterRate)
    return ContinuousLayer{Resource, Precipitation, _WaterRate,
                           Matrix{_WaterRate}}(mat, 1, _SUPPLY_SIZE,
                                               _supply_static())
end
# The raw raster is a per-area rate at its own native reporting period (e.g. `L*m^-2*yr^-1`
# for `bio12`) — `cancel` dispatches on dimension alone (any native time unit) and converts
# straight to the canonical absolute `L/day`, against the grid's (area-preserving, see
# `_cellsize`) cell area; no separate pre-`uconvert` step needed.
function WaterSupply(bioclim::ClimateRaster{WorldClim{BioClim}})
    mat = Matrix(bioclim.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    cellarea = _cellsize(bioclim.array)^2
    return WaterSupply(cancel.(mat, Ref(cellarea)))
end
function SolarTimeSupply(mat::Array{_SolarRate, 3}, time::Int64)
    mat[isnan.(mat)] .= zero(_SolarRate)
    return ContinuousLayer{Resource, SolarRadiation, _SolarRate,
                           Array{_SolarRate, 3}}(mat, time, _SUPPLY_SIZE,
                                                 _supply_cyclic())
end
# Same per-cell-area treatment as `WaterSupply(::ClimateRaster{WorldClim{BioClim}})` above.
function SolarTimeSupply(worldclim::ClimateRaster{WorldClim{Climate}},
                         time::Int64)
    mat = Array(worldclim.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    cellarea = _cellsize(worldclim.array)^2
    return SolarTimeSupply(cancel.(mat, Ref(cellarea)), time)
end
function WaterTimeSupply(mat::Array{_WaterRate, 3}, time::Int64)
    mat[isnan.(mat)] .= zero(_WaterRate)
    return ContinuousLayer{Resource, Precipitation, _WaterRate,
                           Array{_WaterRate, 3}}(mat, time, _SUPPLY_SIZE,
                                                 _supply_cyclic())
end
# Same per-cell-area treatment as `WaterSupply(::ClimateRaster{WorldClim{BioClim}})` above.
function WaterTimeSupply(worldclim::ClimateRaster{WorldClim{Climate}},
                         time::Int64)
    mat = Array(worldclim.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    cellarea = _cellsize(worldclim.array)^2
    return WaterTimeSupply(cancel.(mat, Ref(cellarea)), time)
end
# --- Recipes --------------------------------------------------------------------------
@recipe function f(B::ContinuousLayer{Resource, A, V}) where {A, V}
    b = ustrip.(B.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> _resourcetitle(V)
    clims --> (minimum(b) * 0.99, maximum(b) * 1.01)
    return b
end
@recipe function f(B::ContinuousLayer{Resource, A, V, Arr},
                   time::Int64) where {A, V, Arr <: AbstractArray{V, 3}}
    b = ustrip.(B.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> _resourcetitle(V)
    clims --> (minimum(b[:, :, time]) * 0.99, maximum(b[:, :, time]) * 1.01)
    return b[:, :, time]
end
@recipe function f(H::LayerCollection2{Resource, B1, B2}) where {B1, B2}
    x, y = H.one, H.two
    layout := 2
    @series begin
        subplot := 1
        x
    end
    @series begin
        subplot := 2
        y
    end
end
