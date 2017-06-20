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
mutable struct PopGrowth <: AbstractParams
  birth::Vector{Float64}
  death::Vector{Float64}
  l::Float64
  s::Float64

  function PopGrowth(birth::Vector{Float64}, death::Vector{Float64},
    l::Float64, s::Float64)
    l > s || error("l must be greater than s")
    l >= 0 && s >= 0 || error("l and s must be greater than zero")
    new(birth, death, l, s)
  end
end

mutable struct EqualPop <: AbstractParams
  birth::Float64
  death::Float64
  l::Float64
  s::Float64
end

function equalpop(params::EqualPop, numspp)
  PopGrowth(repmat([params.birth], numspp), repmat([params.death], numspp),
  params.l, params.s)
end
