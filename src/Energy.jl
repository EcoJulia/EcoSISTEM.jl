using RecipesBase
using Unitful
using Unitful.DefaultSymbols
import Base: eltype, length
"""
    Abstract1Requirement{Energy}

Abstract supertype for all species energy requirement types, parameterised by
the type(s) of energy required `Energy`.
"""
abstract type AbstractRequirement{Energy} end
abstract type Abstract1Requirement{Energy} <: AbstractRequirement{Energy} end
abstract type Abstract2Requirements{Energy} <: AbstractRequirement{Energy} end

numrequirements(::Type{<: Abstract1Requirement}) = 1
numrequirements(::Type{<: Abstract2Requirements}) = 2


function eltype(::Abstract1Requirement{Energy}) where Energy
  return Energy
end

"""
    SimpleRequirement <: Abstract1Requirement{Float64}

A simple energy requirement is a single float for each species.
"""
mutable struct SimpleRequirement <: Abstract1Requirement{Float64}
  energy::Vector{Float64}
  exchange_rate::Float64

  function SimpleRequirement(energy::Vector{Float64})
      return new(energy, 1.0)
  end
end



length(req::SimpleRequirement) = length(req.energy)

function _getenergyusage(abun::Vector{Int64}, req::SimpleRequirement)
    sum(abun .* req.energy)
end
GLOBAL_typedict["SimpleRequirement"] = SimpleRequirement
"""
    SizeRequirement <: Abstract1Requirement{Float64}

A simple energy requirement is a single float for each species.
"""
mutable struct SizeRequirement <: Abstract1Requirement{Float64}
  energy::Vector{Float64}
  pop_mass_rel::Float64
  area::Unitful.Area
  exchange_rate::Float64

  function SizeRequirement(energy::Vector{Float64}, pop_mass_rel::Float64, area::Unitful.Area, exchange_rate::Float64 = 1.0)
      return new(energy, pop_mass_rel, area, exchange_rate)
  end
end

length(req::SizeRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::SizeRequirement)
    sum(abun .* req.energy)
end
GLOBAL_typedict["SizeRequirement"] = SizeRequirement
"""
    SolarRequirement <: Abstract1Requirement{typeof(1.0*day^-1*kJ*m^-2)}

"""
mutable struct SolarRequirement <: Abstract1Requirement{typeof(1.0*kJ)}
  energy::Vector{typeof(1.0*kJ)}
  exchange_rate::typeof(1.0/kJ)

  function SolarRequirement(energy::Vector{typeof(1.0*kJ)}, exchange_rate::typeof(1.0/kJ) = 1.0/kJ)
      return new(energy, exchange_rate)
  end
end

length(req::SolarRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::SolarRequirement)
    sum(abun .* req.energy)
end
GLOBAL_typedict["SolarRequirement"] = SolarRequirement
"""
    WaterRequirement <: Abstract1Requirement{typeof(1.0*mm)}

"""
mutable struct WaterRequirement <: Abstract1Requirement{typeof(1.0*mm)}
  energy::Vector{typeof(1.0*mm)}
  exchange_rate::typeof(1.0/mm)

  function WaterRequirement(energy::Vector{typeof(1.0*mm)}, exchange_rate::typeof(1.0/mm) = 1.0/mm)
      return new(energy, exchange_rate)
  end
end
length(req::WaterRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::WaterRequirement)
    sum(abun .* req.energy)
end
GLOBAL_typedict["WaterRequirement"] = WaterRequirement

mutable struct ReqCollection2{R1, R2} <: Abstract2Requirements{Tuple{R1, R2}}
    r1::R1
    r2::R2
end
length(req::ReqCollection2) = length(req.r1.energy)
function eltype(req::ReqCollection2)
    return [eltype(req.r1), eltype(req.r2)]
end
function _getenergyusage(abun::Vector{Int64}, req::ReqCollection2)
    [_getenergyusage(abun, req.r1), _getenergyusage(abun, req.r2)]
end

GLOBAL_typedict["ReqCollection2"] = ReqCollection2
unitdict= Dict(kJ => "Solar Radiation (kJ)",NoUnits => "Free energy",
    mm => "Available water (mm)")
"""
    AbstractBudget

Abstract supertype for all budget types
"""
abstract type AbstractBudget{Requirement} end
abstract type AbstractTimeBudget{Requirement} <: AbstractBudget{Requirement} end

function eltype(::AbstractBudget{Energy}) where Energy
  return Energy
end

function countsubcommunities(ab::AbstractBudget)
  return _countsubcommunities(ab)
end
@recipe function f(B::AbstractBudget{R}) where R
    b = ustrip.(B.matrix)
    seriestype  :=  :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(R)]
    clims --> (minimum(b) * 0.99, maximum(b) * 1.01)
    b
end
"""
    SimpleBudget <: AbstractBudget{Float64}

This budget type has a matrix of floats, representing the energy budget of each
subcommunity in the abiotic environment.
"""
mutable struct SimpleBudget <: AbstractBudget{Unitful.Quantity{Float64, 𝐋^-2}}
  matrix::Matrix{Unitful.Quantity{Float64, 𝐋^-2}}

end

function _countsubcommunities(bud::SimpleBudget)
  return length(bud.matrix)
end
function _getbudget(bud::SimpleBudget)
    return bud.matrix
end
function _getavailableenergy(bud::SimpleBudget)
    return sum(bud.matrix[.!isnan.(bud.matrix)])
end
GLOBAL_typedict["SimpleBudget"] = SimpleBudget

"""
    SolarBudget <: AbstractBudget{typeof(1.0*kJ)}

This budget type has a matrix of solar energy units, representing the energy budget of each
subcommunity in the abiotic environment at a fixed point in time.
"""
mutable struct SolarBudget <: AbstractBudget{typeof(1.0*kJ)}
  matrix::Array{typeof(1.0*kJ), 2}
  function SolarBudget(mat::Array{typeof(1.0*kJ), 2})
    mat[isnan.(mat)] .=  0*kJ
    return new(mat)
  end
end
function _countsubcommunities(bud::SolarBudget)
  return length(bud.matrix)
end

function _getbudget(bud::SolarBudget)
    return bud.matrix
end
function _getavailableenergy(bud::SolarBudget)
    return sum(bud.matrix[.!isnan.(bud.matrix)])
end
GLOBAL_typedict["SolarBudget"] = SolarBudget

"""
    SolarTimeBudget <: AbstractBudget{typeof(1.0*kJ)}

This budget type has a matrix of solar energy units, representing the energy budget of each
subcommunity in the abiotic environment along with which time dimension we are interested in.
"""
mutable struct SolarTimeBudget <: AbstractTimeBudget{typeof(1.0*kJ)}
  matrix::Array{typeof(1.0*kJ), 3}
  time::Int64
  function SolarTimeBudget(mat::Array{typeof(1.0*kJ), 3}, time::Int64)
    mat[isnan.(mat)] .=  0*kJ
    return new(mat, time)
  end
end
function _countsubcommunities(bud::SolarTimeBudget)
  return length(bud.matrix[:,:,1])
end

function _getbudget(bud::SolarTimeBudget)
    return bud.matrix[:, :, bud.time]
end
function _getavailableenergy(bud::SolarTimeBudget)
    return sum(bud.matrix[.!isnan.(bud.matrix)])
end
@recipe function f(B::SolarTimeBudget, time::Int64)
    b = ustrip.(B.matrix)
    seriestype  :=  :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[kJ]
    clims --> (minimum(b[:, :, time]) * 0.99, maximum(b[:, :, time]) * 1.01)
    b[:, :, time]
end
GLOBAL_typedict["SolarTimeBudget"] = SolarTimeBudget

"""
    WaterBudget <: AbstractBudget{typeof(1.0*mm)}

This budget type has a matrix of solar energy units, representing the energy budget of each
subcommunity in the abiotic environment at a fixed point in time.
"""
mutable struct WaterBudget <: AbstractBudget{typeof(1.0*mm)}
  matrix::Array{typeof(1.0*mm), 2}
  function WaterBudget(mat::Array{typeof(1.0*mm), 2})
    mat[isnan.(mat)] .=  0*mm
    return new(mat)
  end
end
function _countsubcommunities(bud::WaterBudget)
  return length(bud.matrix)
end

function _getbudget(bud::WaterBudget)
    return bud.matrix
end
function _getavailableenergy(bud::WaterBudget)
    return sum(bud.matrix[.!isnan.(bud.matrix)])
end
GLOBAL_typedict["WaterBudget"] = WaterBudget

"""
    WaterTimeBudget <: AbstractBudget{typeof(1.0*mm)}

This budget type has a matrix of solar energy units, representing the energy budget of each
subcommunity in the abiotic environment along with which time dimension we are interested in.
"""
mutable struct WaterTimeBudget <: AbstractTimeBudget{typeof(1.0*mm)}
  matrix::Array{typeof(1.0*mm), 3}
  time::Int64
  function WaterTimeBudget(mat::Array{typeof(1.0*mm), 3}, time::Int64)
    mat[isnan.(mat)] .=  0*mm
    return new(mat, time)
  end
end

function _countsubcommunities(bud::WaterTimeBudget)
  return length(bud.matrix[:,:,1])
end

function _getbudget(bud::WaterTimeBudget)
    return bud.matrix[:, :, bud.time]
end
function _getavailableenergy(bud::WaterTimeBudget)
    return sum(bud.matrix[.!isnan.(bud.matrix)])
end
@recipe function f(B::WaterTimeBudget, time::Int64)
    b = ustrip.(B.matrix)
    seriestype  :=  :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[mm]
    clims --> (minimum(b[:, :, time]) * 0.99, maximum(b[:, :, time]) * 1.01)
    b[:, :, time]
end
GLOBAL_typedict["WaterTimeBudget"] = WaterTimeBudget

mutable struct BudgetCollection2{B1, B2} <: AbstractBudget{Tuple{B1, B2}}
    b1::B1
    b2::B2
end

function eltype(bud::BudgetCollection2)
    return [eltype(bud.b1), eltype(bud.b2)]
end
function _countsubcommunities(bud::BudgetCollection2)
  return length(bud.b1.matrix[:,:,1])
end

function _getbudget(bud::BudgetCollection2, field::Symbol)
    B = getfield(bud, field)
    return _getbudget(B)
end

function _getbudget(bud::Bud, field::Symbol) where Bud <: AbstractTimeBudget
    B = getfield(bud, field)
    return B.matrix[:, :, B.time]
end


function _getavailableenergy(bud::BudgetCollection2)
    return [_getavailableenergy(bud.b1), _getavailableenergy(bud.b2)]
end

@recipe function f(H::BudgetCollection2{B1, B2}) where {B1, B2}
    x, y = B.b1, B.b2
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

GLOBAL_typedict["BudgetCollection2"] = BudgetCollection2
