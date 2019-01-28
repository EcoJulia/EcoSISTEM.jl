using MyUnitful
"""
    AbstractParams

Abstract supertype for all simulation parameter types
"""
abstract type AbstractParams end

if VERSION >= v"0.7"
    """
        PopGrowth <: AbstractParams

    Basic parameter type that holds information on a population's birth and death
    rates, `birth` and `death`, as well as how these are altered by energy
    availability. `l` represents the longevity of species based on their energy
    requirements and `s` is the survival of species dependent on how well their
    traits reflect the environment.
    """
    mutable struct PopGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{Quantity{Float64, 𝐓^-1, U}}
      death::Vector{Quantity{Float64, 𝐓^-1, U}}
      l::Float64
      s::Float64
      boost::Float64

      function PopGrowth{U}(birth::Vector{Quantity{Float64, 𝐓^-1, U}},
        death::Vector{Quantity{Float64, 𝐓^-1, U}},
        l::Float64, s::Float64, boost::Float64) where {U <: Unitful.Units}
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
      birth::Quantity{Float64, 𝐓^-1}
      death::Quantity{Float64, 𝐓^-1}
      l::Float64
      s::Float64
      boost::Float64
    end

    mutable struct NoGrowth{U <: Unitful.Units} <: AbstractParams
        birth::Vector{Quantity{Float64, 𝐓^-1, U}}
        death::Vector{Quantity{Float64, 𝐓^-1, U}}
        l::Float64
        s::Float64
        boost::Float64

        function NoGrowth{U}(birth::Vector{Quantity{Float64, 𝐓^-1, U}},
          death::Vector{Quantity{Float64, 𝐓^-1, U}},
          l::Float64, s::Float64, boost::Float64) where {U <: Unitful.Units}
          l > s || error("l must be greater than s")
          l >= 0 && s >= 0 || error("l and s must be greater than zero")
          new{U}(birth, death, l, s, boost)
        end
    end
else
    """
        PopGrowth <: AbstractParams

    Basic parameter type that holds information on a population's birth and death
    rates, `birth` and `death`, as well as how these are altered by energy
    availability. `l` represents the longevity of species based on their energy
    requirements and `s` is the survival of species dependent on how well their
    traits reflect the environment.
    """
    mutable struct PopGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{Quantity{Float64, typeof(𝐓^-1), U}}
      death::Vector{Quantity{Float64, typeof(𝐓^-1), U}}
      l::Float64
      s::Float64
      boost::Float64

      function PopGrowth{U}(birth::Vector{Quantity{Float64, typeof(𝐓^-1), U}},
        death::Vector{Quantity{Float64, typeof(𝐓^-1), U}},
        l::Float64, s::Float64, boost::Float64) where {U <: Unitful.Units}
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

    mutable struct NoGrowth{U <: Unitful.Units} <: AbstractParams
        birth::Vector{Quantity{Float64, typeof(𝐓^-1), U}}
        death::Vector{Quantity{Float64, typeof(𝐓^-1), U}}
        l::Float64
        s::Float64
        boost::Float64

        function NoGrowth{U}(birth::Vector{Quantity{Float64, typeof(𝐓^-1), U}},
          death::Vector{Quantity{Float64, typeof(𝐓^-1), U}},
          l::Float64, s::Float64, boost::Float64) where {U <: Unitful.Units}
          l > s || error("l must be greater than s")
          l >= 0 && s >= 0 || error("l and s must be greater than zero")
          new{U}(birth, death, l, s, boost)
        end
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
  params.l, params.s, params.boost)
end
function equalpop(params::PopGrowth, numspp)
  return params
end

function equalpop(params::NoGrowth, numspp)
  return params
end

GLOBAL_funcdict["equalpop"] = equalpop
