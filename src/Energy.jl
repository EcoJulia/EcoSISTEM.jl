# SPDX-License-Identifier: LGPL-3.0-or-later

using EcoSISTEM.ClimatePref
using RecipesBase
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources

import Base: eltype, length
"""
    Abstract1Requirement{Energy}

Abstract supertype for all species energy requirement types, parameterised by
the type(s) of energy required `Energy`.
"""
abstract type AbstractRequirement{Energy} end
abstract type Abstract1Requirement{Energy} <: AbstractRequirement{Energy} end
abstract type Abstract2Requirements{Energy} <: AbstractRequirement{Energy} end

numrequirements(::Type{<:Abstract1Requirement}) = 1
numrequirements(::Type{<:Abstract2Requirements}) = 2

function Base.eltype(::Abstract1Requirement{Energy}) where {Energy}
    return Energy
end

"""
    SimpleRequirement <: Abstract1Requirement{Float64}

A simple energy requirement is a single float for each species.
"""
struct SimpleRequirement <: Abstract1Requirement{Float64}
    energy::Vector{Float64}
    exchange_rate::Float64

    function SimpleRequirement(energy::Vector{Float64})
        return new(energy, 1.0)
    end
end

Base.length(req::SimpleRequirement) = length(req.energy)

function _getenergyusage(abun::Vector{Int64}, req::SimpleRequirement)
    return sum(abun .* req.energy)
end

"""
    SizeRequirement <: Abstract1Requirement{Float64}

A simple energy requirement is a single float for each species.
"""
struct SizeRequirement <: Abstract1Requirement{Float64}
    energy::Vector{Float64}
    pop_mass_rel::Float64
    area::Unitful.Area
    exchange_rate::Float64

    function SizeRequirement(energy::Vector{Float64},
                             pop_mass_rel::Float64,
                             area::Unitful.Area,
                             exchange_rate::Float64 = 1.0)
        return new(energy, pop_mass_rel, area, exchange_rate)
    end
end

Base.length(req::SizeRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::SizeRequirement)
    return sum(abun .* req.energy)
end

"""
    SolarRequirement <: Abstract1Requirement{typeof(1.0*kJ)}

A vector of solar energy requirements (kJ) for each species.
"""
struct SolarRequirement <: Abstract1Requirement{typeof(1.0 * kJ)}
    energy::Vector{typeof(1.0 * kJ)}
    exchange_rate::typeof(1.0 / kJ)

    function SolarRequirement(energy::Vector{<:Unitful.Energy{Float64}},
                              exchange_rate::Unitful.Quantity{Float64} = 1.0 /
                                                                         mean(energy))
        return new(uconvert.(kJ, energy), uconvert(kJ^-1, exchange_rate))
    end
end

Base.length(req::SolarRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::SolarRequirement)
    return sum(abun .* req.energy)
end

"""
    WaterRequirement <: Abstract1Requirement{typeof(1.0*mm)}

A vector of water requirements (mm) for each species.
"""
struct WaterRequirement <: Abstract1Requirement{typeof(1.0 * mm)}
    energy::Vector{typeof(1.0 * mm)}
    exchange_rate::typeof(1.0 / mm)

    function WaterRequirement(energy::Vector{<:Unitful.Length{Float64}},
                              exchange_rate::Unitful.Quantity{Float64} = 1.0 /
                                                                         mean(energy))
        return new(uconvert.(mm, energy), uconvert.(mm^-1, exchange_rate))
    end
end
Base.length(req::WaterRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::WaterRequirement)
    return sum(abun .* req.energy)
end

"""
    VolWaterRequirement <: Abstract1Requirement{typeof(1.0*mm)}
A vector of soil water volume requirements (m^3) for each species.
"""
struct VolWaterRequirement <: Abstract1Requirement{typeof(1.0 * m^3)}
    energy::Vector{typeof(1.0 * m^3)}
    exchange_rate::typeof(1.0 / m^3)

    function VolWaterRequirement(energy::Vector{<:Unitful.Volume{Float64}},
                                 exchange_rate::Unitful.Quantity{Float64} = 1.0 /
                                                                            mean(energy))
        return new(uconvert.(m^3, energy), uconvert.(m^-3, exchange_rate))
    end
end
Base.length(req::VolWaterRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::VolWaterRequirement)
    return sum(abun .* req.energy)
end

struct ReqCollection2{R1, R2} <: Abstract2Requirements{Tuple{R1, R2}}
    one::R1
    two::R2
end
Base.length(req::ReqCollection2) = length(req.one.energy)
function Base.eltype(req::ReqCollection2)
    return [eltype(req.one), eltype(req.two)]
end
function _getenergyusage(abun::Vector{Int64}, req::ReqCollection2)
    return [_getenergyusage(abun, req.one), _getenergyusage(abun, req.two)]
end

unitdict = Dict(kJ => "Solar Radiation (kJ)",
                NoUnits => "Free energy",
                mm => "Available water (mm)")

# ---------------------------------------------------------------------------
# Budgets — `Budget`-role layers (folded onto `ContinuousLayer{Budget}`)
# ---------------------------------------------------------------------------
# The concrete budget names (`SimpleBudget`/`SolarBudget`/… and the time-varying
# `*TimeBudget`s) are aliases (in Layer.jl) over `ContinuousLayer{Budget, axis, V, Arr}`:
# a static budget is 2-D (`Matrix`), a time budget 3-D (`Array`, indexed by `time`). The
# constructors below fill the (unused) `size` and the per-timestep `change` rule and zero
# NaNs, reproducing the old budget structs. A budget's `size` is never read
# (geometry/dispersal use the habitat), so a placeholder is stored; `NoChange`/`cyclicChange`
# live in HabitatUpdate.jl (included later) and resolve at call time.
const _BUDGET_SIZE = 1.0m
_budget_static() = HabitatUpdate(NoChange, 0.0 / s, Unitful.Dimensions{()})
_budget_cyclic() = HabitatUpdate(cyclicChange, 0.0 / s, Unitful.Dimensions{()})

Base.eltype(::ContinuousLayer{Budget, A, V}) where {A, V} = V
function Base.eltype(bud::LayerCollection2{Budget})
    return [eltype(bud.one), eltype(bud.two)]
end

countsubcommunities(ab::AbstractBudget) = _countsubcommunities(ab)
function _countsubcommunities(bud::ContinuousLayer{Budget, A, V, Arr}) where {A,
                                                                              V,
                                                                              Arr <:
                                                                              AbstractMatrix{V}}
    return length(bud.matrix)
end
function _countsubcommunities(bud::ContinuousLayer{Budget, A, V, Arr}) where {A,
                                                                              V,
                                                                              Arr <:
                                                                              AbstractArray{V,
                                                                                            3}}
    return length(@view bud.matrix[:, :, 1])
end
function _countsubcommunities(bud::LayerCollection2{Budget})
    return _countsubcommunities(bud.one)
end

# The resource available in each cell: the full matrix (static), or the current time slice.
function _getbudget(bud::ContinuousLayer{Budget, A, V, Arr}) where {A, V,
                                                                    Arr <:
                                                                    AbstractMatrix{V}}
    return bud.matrix
end
function _getbudget(bud::ContinuousLayer{Budget, A, V, Arr}) where {A, V,
                                                                    Arr <:
                                                                    AbstractArray{V,
                                                                                  3}}
    return @view bud.matrix[:, :, bud.time]
end
function _getbudget(bud::LayerCollection2{Budget}, field::Symbol)
    return _getbudget(getfield(bud, field))
end

function _getavailableenergy(bud::ContinuousLayer{Budget})
    return sum(bud.matrix[.!isnan.(bud.matrix)])
end
function _getavailableenergy(bud::LayerCollection2{Budget})
    return [_getavailableenergy(bud.one), _getavailableenergy(bud.two)]
end

# --- Constructors reproducing the old per-type budget structs -------------------------
function SimpleBudget(mat::Matrix{Float64})
    return ContinuousLayer{Budget, Unclassified, Float64, Matrix{Float64}}(mat,
                                                                           1,
                                                                           _BUDGET_SIZE,
                                                                           _budget_static())
end
function SolarBudget(mat::Matrix{typeof(1.0 * kJ)})
    mat[isnan.(mat)] .= 0 * kJ
    return ContinuousLayer{Budget, SolarRadiation, typeof(1.0 * kJ),
                           Matrix{typeof(1.0 * kJ)}}(mat, 1, _BUDGET_SIZE,
                                                     _budget_static())
end
function WaterBudget(mat::Matrix{typeof(1.0 * mm)})
    mat[isnan.(mat)] .= 0.0mm
    return ContinuousLayer{Budget, Precipitation, typeof(1.0 * mm),
                           Matrix{typeof(1.0 * mm)}}(mat, 1, _BUDGET_SIZE,
                                                     _budget_static())
end
function WaterBudget(bc::ClimateRaster{WorldClim{BioClim}})
    mat = Matrix(bc.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    return WaterBudget(mat)
end
function VolWaterBudget(mat::Matrix{typeof(1.0 * m^3)})
    mat[isnan.(mat)] .= 0 * m^3
    return ContinuousLayer{Budget, VolumetricWater, typeof(1.0 * m^3),
                           Matrix{typeof(1.0 * m^3)}}(mat, 1, _BUDGET_SIZE,
                                                      _budget_static())
end
function SolarTimeBudget(mat::Array{typeof(1.0 * kJ), 3}, time::Int64)
    mat[isnan.(mat)] .= 0 * kJ
    return ContinuousLayer{Budget, SolarRadiation, typeof(1.0 * kJ),
                           Array{typeof(1.0 * kJ), 3}}(mat, time, _BUDGET_SIZE,
                                                       _budget_cyclic())
end
function SolarTimeBudget(wc::Worldclim_monthly, time::Int64)
    mat = Array(wc.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    return SolarTimeBudget(mat, time)
end
function WaterTimeBudget(mat::Array{typeof(1.0 * mm), 3}, time::Int64)
    mat[isnan.(mat)] .= 0 * mm
    return ContinuousLayer{Budget, Precipitation, typeof(1.0 * mm),
                           Array{typeof(1.0 * mm), 3}}(mat, time, _BUDGET_SIZE,
                                                       _budget_cyclic())
end
function WaterTimeBudget(wc::Worldclim_monthly, time::Int64)
    mat = Array(wc.array)
    mat[isnan.(mat)] .= zero(eltype(mat))
    return WaterTimeBudget(mat, time)
end
function VolWaterTimeBudget(mat::Array{typeof(1.0 * m^3), 3}, time::Int64)
    mat[isnan.(mat)] .= 0 * m^3
    return ContinuousLayer{Budget, VolumetricWater, typeof(1.0 * m^3),
                           Array{typeof(1.0 * m^3), 3}}(mat, time, _BUDGET_SIZE,
                                                        _budget_cyclic())
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
