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
    SimpleDemand <: Abstract1Demand{Float64}

A simple resource demand is a single float for each species.
"""
struct SimpleDemand <: Abstract1Demand{Float64}
    resource::Vector{Float64}
    exchange_rate::Float64

    function SimpleDemand(resource::Vector{Float64})
        return new(resource, 1.0)
    end
end

Base.length(demand::SimpleDemand) = length(demand.resource)

function _getdemand(abun::Vector{Int64}, demand::SimpleDemand)
    return sum(abun .* demand.resource)
end

"""
    SizeDemand <: Abstract1Demand{Float64}

A simple resource demand is a single float for each species.
"""
struct SizeDemand <: Abstract1Demand{Float64}
    resource::Vector{Float64}
    pop_mass_rel::Float64
    area::Unitful.Area
    exchange_rate::Float64

    function SizeDemand(resource::Vector{Float64},
                        pop_mass_rel::Float64,
                        area::Unitful.Area,
                        exchange_rate::Float64 = 1.0)
        return new(resource, pop_mass_rel, area, exchange_rate)
    end
end

Base.length(demand::SizeDemand) = length(demand.resource)
function _getdemand(abun::Vector{Int64}, demand::SizeDemand)
    return sum(abun .* demand.resource)
end

"""
    SolarDemand <: Abstract1Demand{typeof(1.0*kJ)}

A vector of solar resource demands (kJ) for each species.
"""
struct SolarDemand <: Abstract1Demand{typeof(1.0 * kJ)}
    resource::Vector{typeof(1.0 * kJ)}
    exchange_rate::typeof(1.0 / kJ)

    function SolarDemand(resource::Vector{<:Unitful.Energy{Float64}},
                         exchange_rate::Unitful.Quantity{Float64} = 1.0 /
                                                                    mean(resource))
        return new(uconvert.(kJ, resource), uconvert(kJ^-1, exchange_rate))
    end
end

Base.length(demand::SolarDemand) = length(demand.resource)
function _getdemand(abun::Vector{Int64}, demand::SolarDemand)
    return sum(abun .* demand.resource)
end

"""
    WaterDemand <: Abstract1Demand{typeof(1.0*mm)}

A vector of water demands (mm) for each species.
"""
struct WaterDemand <: Abstract1Demand{typeof(1.0 * mm)}
    resource::Vector{typeof(1.0 * mm)}
    exchange_rate::typeof(1.0 / mm)

    function WaterDemand(resource::Vector{<:Unitful.Length{Float64}},
                         exchange_rate::Unitful.Quantity{Float64} = 1.0 /
                                                                    mean(resource))
        return new(uconvert.(mm, resource), uconvert.(mm^-1, exchange_rate))
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

unitdict = Dict(kJ => "Solar Radiation (kJ)",
                NoUnits => "Free resource",
                mm => "Available water (mm)")

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
function SimpleSupply(mat::Matrix{Float64})
    return ContinuousLayer{Resource, Unclassified, Float64, Matrix{Float64}}(mat,
                                                                             1,
                                                                             _SUPPLY_SIZE,
                                                                             _supply_static())
end
function SolarSupply(mat::Matrix{typeof(1.0 * kJ)})
    mat[isnan.(mat)] .= 0 * kJ
    return ContinuousLayer{Resource, SolarRadiation, typeof(1.0 * kJ),
                           Matrix{typeof(1.0 * kJ)}}(mat, 1, _SUPPLY_SIZE,
                                                     _supply_static())
end
function WaterSupply(mat::Matrix{typeof(1.0 * mm)})
    mat[isnan.(mat)] .= 0.0mm
    return ContinuousLayer{Resource, Precipitation, typeof(1.0 * mm),
                           Matrix{typeof(1.0 * mm)}}(mat, 1, _SUPPLY_SIZE,
                                                     _supply_static())
end
function WaterSupply(bioclim::ClimateRaster{WorldClim{BioClim}})
    mat = Matrix(bioclim.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    return WaterSupply(mat)
end
function SolarTimeSupply(mat::Array{typeof(1.0 * kJ), 3}, time::Int64)
    mat[isnan.(mat)] .= 0 * kJ
    return ContinuousLayer{Resource, SolarRadiation, typeof(1.0 * kJ),
                           Array{typeof(1.0 * kJ), 3}}(mat, time, _SUPPLY_SIZE,
                                                       _supply_cyclic())
end
function SolarTimeSupply(worldclim::ClimateRaster{WorldClim{Climate}},
                         time::Int64)
    mat = Array(worldclim.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    return SolarTimeSupply(mat, time)
end
function WaterTimeSupply(mat::Array{typeof(1.0 * mm), 3}, time::Int64)
    mat[isnan.(mat)] .= 0 * mm
    return ContinuousLayer{Resource, Precipitation, typeof(1.0 * mm),
                           Array{typeof(1.0 * mm), 3}}(mat, time, _SUPPLY_SIZE,
                                                       _supply_cyclic())
end
function WaterTimeSupply(worldclim::ClimateRaster{WorldClim{Climate}},
                         time::Int64)
    mat = Array(worldclim.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    return WaterTimeSupply(mat, time)
end
# --- Recipes --------------------------------------------------------------------------
@recipe function f(B::ContinuousLayer{Resource, A, V}) where {A, V}
    b = ustrip.(B.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(V)]
    clims --> (minimum(b) * 0.99, maximum(b) * 1.01)
    return b
end
@recipe function f(B::ContinuousLayer{Resource, A, V, Arr},
                   time::Int64) where {A, V, Arr <: AbstractArray{V, 3}}
    b = ustrip.(B.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(V)]
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
