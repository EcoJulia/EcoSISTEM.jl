using myunitful
"""
    AbstractParams

Abstract supertype for all simulation parameter types
"""
abstract type AbstractParams end
"""
    PopGrowth <: AbstractParams

Basic parameter type that holds information on a population's birth and death
rates, `birth` and `death`, as well as how these are altered by energy
availability. `l` represents the longevity of species based on their energy
requirements and `s` is the survival of species dependent on how well their
traits reflect the environment.
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
"""
    EqualPop <: AbstractParams

Parameter type that holds information on a population's birth and death
rates, `birth` and `death`, specifically populations where all species have the
same information. `l` represents the longevity of species based on their energy
requirements and `s` is the survival of species dependent on how well their
traits reflect the environment. Finally `boost` is used to manipulate how much
of a boost the species get from being in an environment with lots of available
energy.
"""
mutable struct EqualPop <: AbstractParams
  birth::Quantity{Float64, typeof(𝐓^-1)}
  death::Quantity{Float64, typeof(𝐓^-1)}
  l::Float64
  s::Float64
  boost::Float64
end
"""
    equalpop(params::EqualPop, numspp)

Function that takes demographic parameters from type `EqualPop` and converts
them into type `PopGrowth` based on the number of species (`numspp`).
"""
function equalpop(params::EqualPop, numspp)
  u = unit(params.birth)
  PopGrowth{typeof(u)}(repmat([params.birth], numspp), repmat([params.death], numspp),
  params.l, params.s, params.boost)
end
function equalpop(params::PopGrowth, numspp)
  return params
end
