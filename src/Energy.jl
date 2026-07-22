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

Base.length(dem::SimpleDemand) = length(dem.resource)

function _getdemand(abun::Vector{Int64}, dem::SimpleDemand)
    return sum(abun .* dem.resource)
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

Base.length(dem::SizeDemand) = length(dem.resource)
function _getdemand(abun::Vector{Int64}, dem::SizeDemand)
    return sum(abun .* dem.resource)
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

Base.length(dem::SolarDemand) = length(dem.resource)
function _getdemand(abun::Vector{Int64}, dem::SolarDemand)
    return sum(abun .* dem.resource)
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
Base.length(dem::WaterDemand) = length(dem.resource)
function _getdemand(abun::Vector{Int64}, dem::WaterDemand)
    return sum(abun .* dem.resource)
end

"""
    VolWaterDemand <: Abstract1Demand{typeof(1.0*mm)}
A vector of soil water volume demands (m^3) for each species.
"""
struct VolWaterDemand <: Abstract1Demand{typeof(1.0 * m^3)}
    resource::Vector{typeof(1.0 * m^3)}
    exchange_rate::typeof(1.0 / m^3)

    function VolWaterDemand(resource::Vector{<:Unitful.Volume{Float64}},
                            exchange_rate::Unitful.Quantity{Float64} = 1.0 /
                                                                       mean(resource))
        return new(uconvert.(m^3, resource), uconvert.(m^-3, exchange_rate))
    end
end
Base.length(dem::VolWaterDemand) = length(dem.resource)
function _getdemand(abun::Vector{Int64}, dem::VolWaterDemand)
    return sum(abun .* dem.resource)
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
Base.length(dem::DemandCollection2) = length(dem.one.resource)
function Base.eltype(dem::DemandCollection2)
    return [eltype(dem.one), eltype(dem.two)]
end
function _getdemand(abun::Vector{Int64}, dem::DemandCollection2)
    return [_getdemand(abun, dem.one), _getdemand(abun, dem.two)]
end

unitdict = Dict(kJ => "Solar Radiation (kJ)",
                NoUnits => "Free resource",
                mm => "Available water (mm)")

# ---------------------------------------------------------------------------
# Supplies — `Budget`-role layers (folded onto `ContinuousLayer{Budget}`)
# ---------------------------------------------------------------------------
# The concrete supply names (`SimpleSupply`/`SolarSupply`/… and the time-varying
# `*TimeSupply`s) are aliases (in Layer.jl) over `ContinuousLayer{Budget, axis, V, Arr}`:
# a static supply is 2-D (`Matrix`), a time supply 3-D (`Array`, indexed by `time`). The
# constructors below fill the (unused) `size` and the per-timestep `dynamics` rule and zero
# NaNs, reproducing the old supply structs. A supply's `size` is never read
# (geometry/dispersal use the habitat), so a placeholder is stored; `NoChange`/`cyclicChange`
# live in HabitatUpdate.jl (included later) and resolve at call time.
const _SUPPLY_SIZE = 1.0m
_supply_static() = HabitatUpdate(NoChange, 0.0 / s, Unitful.Dimensions{()})
_supply_cyclic() = HabitatUpdate(cyclicChange, 0.0 / s, Unitful.Dimensions{()})

Base.eltype(::ContinuousLayer{Budget, A, V}) where {A, V} = V
function Base.eltype(sup::LayerCollection2{Budget})
    return [eltype(sup.one), eltype(sup.two)]
end

countsubcommunities(ab::AbstractSupply) = _countsubcommunities(ab)
function _countsubcommunities(sup::ContinuousLayer{Budget, A, V, Arr}) where {A,
                                                                              V,
                                                                              Arr <:
                                                                              AbstractMatrix{V}}
    return length(sup.matrix)
end
function _countsubcommunities(sup::ContinuousLayer{Budget, A, V, Arr}) where {A,
                                                                              V,
                                                                              Arr <:
                                                                              AbstractArray{V,
                                                                                            3}}
    return length(@view sup.matrix[:, :, 1])
end
function _countsubcommunities(sup::LayerCollection2{Budget})
    return _countsubcommunities(sup.one)
end

# The resource available in each cell: the full matrix (static), or the current time slice.
function _getsupply(sup::ContinuousLayer{Budget, A, V, Arr}) where {A, V,
                                                                    Arr <:
                                                                    AbstractMatrix{V}}
    return sup.matrix
end
function _getsupply(sup::ContinuousLayer{Budget, A, V, Arr}) where {A, V,
                                                                    Arr <:
                                                                    AbstractArray{V,
                                                                                  3}}
    return @view sup.matrix[:, :, sup.time]
end
function _getsupply(sup::LayerCollection2{Budget}, field::Symbol)
    return _getsupply(getfield(sup, field))
end

function _getavailablesupply(sup::ContinuousLayer{Budget})
    return sum(sup.matrix[.!isnan.(sup.matrix)])
end
function _getavailablesupply(sup::LayerCollection2{Budget})
    return [_getavailablesupply(sup.one), _getavailablesupply(sup.two)]
end

# --- Constructors reproducing the old per-type supply structs -------------------------
function SimpleSupply(mat::Matrix{Float64})
    return ContinuousLayer{Budget, Unclassified, Float64, Matrix{Float64}}(mat,
                                                                           1,
                                                                           _SUPPLY_SIZE,
                                                                           _supply_static())
end
function SolarSupply(mat::Matrix{typeof(1.0 * kJ)})
    mat[isnan.(mat)] .= 0 * kJ
    return ContinuousLayer{Budget, SolarRadiation, typeof(1.0 * kJ),
                           Matrix{typeof(1.0 * kJ)}}(mat, 1, _SUPPLY_SIZE,
                                                     _supply_static())
end
function WaterSupply(mat::Matrix{typeof(1.0 * mm)})
    mat[isnan.(mat)] .= 0.0mm
    return ContinuousLayer{Budget, Precipitation, typeof(1.0 * mm),
                           Matrix{typeof(1.0 * mm)}}(mat, 1, _SUPPLY_SIZE,
                                                     _supply_static())
end
function WaterSupply(bc::ClimateRaster{WorldClim{BioClim}})
    mat = Matrix(bc.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    return WaterSupply(mat)
end
function VolWaterSupply(mat::Matrix{typeof(1.0 * m^3)})
    mat[isnan.(mat)] .= 0 * m^3
    return ContinuousLayer{Budget, VolumetricWater, typeof(1.0 * m^3),
                           Matrix{typeof(1.0 * m^3)}}(mat, 1, _SUPPLY_SIZE,
                                                      _supply_static())
end
function SolarTimeSupply(mat::Array{typeof(1.0 * kJ), 3}, time::Int64)
    mat[isnan.(mat)] .= 0 * kJ
    return ContinuousLayer{Budget, SolarRadiation, typeof(1.0 * kJ),
                           Array{typeof(1.0 * kJ), 3}}(mat, time, _SUPPLY_SIZE,
                                                       _supply_cyclic())
end
function SolarTimeSupply(wc::ClimateRaster{WorldClim{Climate}}, time::Int64)
    mat = Array(wc.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    return SolarTimeSupply(mat, time)
end
function WaterTimeSupply(mat::Array{typeof(1.0 * mm), 3}, time::Int64)
    mat[isnan.(mat)] .= 0 * mm
    return ContinuousLayer{Budget, Precipitation, typeof(1.0 * mm),
                           Array{typeof(1.0 * mm), 3}}(mat, time, _SUPPLY_SIZE,
                                                       _supply_cyclic())
end
function WaterTimeSupply(wc::ClimateRaster{WorldClim{Climate}}, time::Int64)
    mat = Array(wc.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    return WaterTimeSupply(mat, time)
end
function VolWaterTimeSupply(mat::Array{typeof(1.0 * m^3), 3}, time::Int64)
    mat[isnan.(mat)] .= 0 * m^3
    return ContinuousLayer{Budget, VolumetricWater, typeof(1.0 * m^3),
                           Array{typeof(1.0 * m^3), 3}}(mat, time, _SUPPLY_SIZE,
                                                        _supply_cyclic())
end

# --- Recipes --------------------------------------------------------------------------
@recipe function f(B::ContinuousLayer{Budget, A, V}) where {A, V}
    b = ustrip.(B.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(V)]
    clims --> (minimum(b) * 0.99, maximum(b) * 1.01)
    return b
end
@recipe function f(B::ContinuousLayer{Budget, A, V, Arr},
                   time::Int64) where {A, V, Arr <: AbstractArray{V, 3}}
    b = ustrip.(B.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(V)]
    clims --> (minimum(b[:, :, time]) * 0.99, maximum(b[:, :, time]) * 1.01)
    return b[:, :, time]
end
@recipe function f(H::LayerCollection2{Budget, B1, B2}) where {B1, B2}
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
