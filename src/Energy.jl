import Base: eltype, length
"""
    AbstractRequirement{Energy}

Abstract supertype for all species energy requirement types, parameterised by
the type(s) of energy required `Energy`.
"""
abstract type AbstractRequirement{Energy} end

function eltype(::AbstractRequirement{Energy}) where Energy
  return Energy
end

"""
    SimpleRequirement <: AbstractRequirement{Float64}

A simple energy requirement is a single float for each species.
"""
mutable struct SimpleRequirement <: AbstractRequirement{Float64}
  energy::Vector{Float64}
end
length(req::SimpleRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::SimpleRequirement)
    sum(abun .* req.energy)
end

"""
    SizeRequirement <: AbstractRequirement{Float64}

A simple energy requirement is a single float for each species.
"""
mutable struct SizeRequirement <: AbstractRequirement{Float64}
  energy::Vector{Float64}
  pop_mass_rel::Float64
  area::Unitful.Area
end
length(req::SizeRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::SizeRequirement)
    sum(abun .* req.energy)
end

"""
    SolarRequirement <: AbstractRequirement{typeof(1.0*day^-1*kJ*m^-2)}

"""
mutable struct SolarRequirement <: AbstractRequirement{typeof(1.0*day^-1*kJ*m^-2)}
  energy::Vector{typeof(1.0*day^-1*kJ*m^-2)}
end
length(req::SolarRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::SolarRequirement)
    sum(abun .* req.energy)
end
"""
    WaterRequirement <: AbstractRequirement{typeof(1.0*mm)}

"""
mutable struct WaterRequirement <: AbstractRequirement{typeof(1.0*mm)}
  energy::Vector{typeof(1.0*mm)}
end
length(req::WaterRequirement) = length(req.energy)
function _getenergyusage(abun::Vector{Int64}, req::WaterRequirement)
    sum(abun .* req.energy)
end

mutable struct ReqCollection2{R1, R2} <: AbstractRequirement{Tuple{R1, R2}}
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

"""
    AbstractBudget

Abstract supertype for all budget types
"""
abstract type AbstractBudget{Requirement} end

function eltype(::AbstractBudget{Energy}) where Energy
  return Energy
end

function countsubcommunities(ab::AbstractBudget)
  return _countsubcommunities(ab)
end

"""
    SimpleBudget <: AbstractBudget{Float64}

This budget type has a matrix of floats, representing the energy budget of each
subcommunity in the abiotic environment.
"""
mutable struct SimpleBudget <: AbstractBudget{Float64}
  matrix::Matrix{Float64}
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

"""
    SolarBudget <: AbstractBudget{typeof(1.0*day^-1*kJ*m^-2)}

This budget type has a matrix of solar energy units, representing the energy budget of each
subcommunity in the abiotic environment along with which time dimension we are interested in.
"""
mutable struct SolarBudget <: AbstractBudget{typeof(1.0*day^-1*kJ*m^-2)}
  matrix::Array{typeof(1.0*day^-1*kJ*m^-2), 3}
  time::Int64
  function SolarBudget(mat::Array{typeof(1.0*day^-1*kJ*m^-2), 3}, time::Int64)
    mat[isnan.(mat)] .=  0*day^-1*kJ*m^-2
    return new(mat, time)
  end
end

function _countsubcommunities(bud::SolarBudget)
  return length(bud.matrix[:,:,1])
end

function _getbudget(bud::SolarBudget)
    return bud.matrix[:, :, bud.time]
end
function _getavailableenergy(bud::SolarBudget)
    return sum(bud.matrix[.!isnan.(bud.matrix)])
end

"""
    WaterBudget <: AbstractBudget{typeof(1.0*mm)}

This budget type has a matrix of solar energy units, representing the energy budget of each
subcommunity in the abiotic environment along with which time dimension we are interested in.
"""
mutable struct WaterBudget <: AbstractBudget{typeof(1.0*mm)}
  matrix::Array{typeof(1.0*mm), 3}
  time::Int64
  function WaterBudget(mat::Array{typeof(1.0*mm), 3}, time::Int64)
    mat[isnan.(mat)] .=  0*mm
    return new(mat, time)
  end
end

function _countsubcommunities(bud::WaterBudget)
  return length(bud.matrix[:,:,1])
end

function _getbudget(bud::WaterBudget)
    return bud.matrix[:, :, bud.time]
end
function _getavailableenergy(bud::WaterBudget)
    return sum(bud.matrix[.!isnan.(bud.matrix)])
end

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
    return B.matrix[:, :, B.time]
end

function _getavailableenergy(bud::BudgetCollection2)
    return [_getavailableenergy(bud.b1), _getavailableenergy(bud.b2)]
end
