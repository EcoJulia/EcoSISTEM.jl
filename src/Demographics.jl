"""
    AbstractParams

Abstract supertype for all simulation parameter types
"""
abstract type AbstractParams end
"""
    PopGrowth <: AbstractParams

Basic parameter type that holds information on a population's birth and death
rates, `birth` and `death`, as well as how these are altered by energy
availability
"""
mutable struct PopGrowth{U <: Rates} <: AbstractParams
  birth::Vector{Quantity{Float64, typeof(𝐓^-1), U}}
  death::Vector{Quantity{Float64, typeof(𝐓^-1), U}}
  l::Float64
  s::Float64
  boost::Float64

  function PopGrowth{U}(birth::Vector{Quantity{Float64, typeof(𝐓^-1), U}},
    death::Vector{Quantity{Float64, typeof(𝐓^-1), U}},
    l::Float64, s::Float64, boost::Float64) where {U <: Rates}
    l > s || error("l must be greater than s")
    l >= 0 && s >= 0 || error("l and s must be greater than zero")
    new{U}(birth, death, l, s, boost)
  end
end

mutable struct EqualPop <: AbstractParams
  birth::Quantity{Float64, typeof(𝐓^-1)}
  death::Quantity{Float64, typeof(𝐓^-1)}
  l::Float64
  s::Float64
  boost::Float64
end

function equalpop(params::EqualPop, numspp)
  u = unit(params.birth)
  PopGrowth{typeof(u)}(repmat([params.birth], numspp), repmat([params.death], numspp),
  params.l, params.s, params.boost)
end
