using MyUnitful
"""
    AbstractParams

Abstract supertype for all simulation parameter types
"""
abstract type AbstractParams end

@static if VERSION >= v"0.7"
    const TimeUnitType{U} = Quantity{Float64, 𝐓^-1, U}
else
    const TimeUnitType{U} = Quantity{Float64, typeof(𝐓^-1), U}
end

    """
        PopGrowth <: AbstractParams

    Basic parameter type that holds information on a population's birth and death
    rates, `birth` and `death`, as well as how these are altered by energy
    availability. `l` represents the longevity of species based on their energy
    requirements and `s` is the survival of species dependent on how well their
    traits reflect the environment.
    """
mutable struct PopGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      longevity::Float64
      survival::Float64
      boost::Float64
    function PopGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}},
        longevity::Float64, survival::Float64, boost::Float64) where {U <: Unitful.Units}
        longevity > survival || error("l must be greater than s")
        longevity >= 0 && survival >= 0 || error("l and s must be greater than zero")
        new{U}(birth, death, longevity, survival, boost)
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
mutable struct EqualPop{U <: Unitful.Units} <: AbstractParams
  birth::TimeUnitType{U}
  death::TimeUnitType{U}
  longevity::Float64
  survival::Float64
  boost::Float64
end

mutable struct NoGrowth{U <: Unitful.Units} <: AbstractParams
    birth::Vector{TimeUnitType{U}}
    death::Vector{TimeUnitType{U}}
    longevity::Float64
    survival::Float64
    boost::Float64

    function NoGrowth{U}(birth::Vector{TimeUnitType{U}},
      death::Vector{TimeUnitType{U}},
      longevity::Float64, survival::Float64, boost::Float64) where {U <: Unitful.Units}
      #longevity > survival || error("longevity must be greater than survival")
      longevity >= 0 && survival >= 0 || error("longevity and survival must be greater than zero")
      new{U}(birth, death, longevity, survival, boost)
    end
end


GLOBAL_typedict["PopGrowth"] = PopGrowth
GLOBAL_typedict["EqualPop"] = EqualPop
GLOBAL_typedict["NoGrowth"] = NoGrowth
"""
    equalpop(params::EqualPop, numspp)

Function that takes demographic parameters from type `EqualPop` and converts
them into type `PopGrowth` based on the number of species (`numspp`).
"""
function equalpop(params::EqualPop, numspp)
  u = unit(params.birth)
  PopGrowth{typeof(u)}(fill(params.birth, numspp), fill(params.death, numspp),
  params.longevity, params.survival, params.boost)
end
function equalpop(params::PopGrowth, numspp)
  return params
end

function equalpop(params::NoGrowth, numspp)
  return params
end

GLOBAL_funcdict["equalpop"] = equalpop
