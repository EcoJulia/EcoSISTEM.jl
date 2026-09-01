# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The per-species demographic rates - birth and death - as the simulation reads them.

using Unitful
using Unitful.DefaultSymbols

"""
    AbstractParams

The demographic rates a simulation runs on: a per-species birth rate and death rate, plus the three
exponents and ceilings that say how those rates respond to the environment. [`EqualPop`](@ref),
[`PopGrowth`](@ref) and [`NoGrowth`](@ref) are the concrete forms.

Every subtype carries the same five fields, and two of them are what make the model a niche model:

  - `longevity` is the exponent on a species' own resource demand, applied to birth **and** death
    alike. It therefore sets the *tempo* of turnover - a body-size proxy - and cancels out of the
    birth/death ratio, so it does not move where a species can persist.
  - `survival` is the exponent on niche suitability, applied with **opposite** signs to birth and
    death. It does move the ratio, which is what makes it the niche.

Rates are `Unitful` quantities per unit time, so a run is timestep-independent: twelve one-month
steps and one twelve-month step see the same rates.
"""
abstract type AbstractParams end

# A rate: a `Float64` per unit time, with the time unit left free so that a set of parameters given
# in per-month and one given in per-year are different types rather than silently mixed.
const TimeUnitType{U} = Quantity{Float64, 𝐓^-1, U}

"""
    EqualPop{U <: Unitful.Units} <: AbstractParams

Demographic rates shared by every species: one birth rate and one death rate for the whole
assemblage. [`equalpop`](@ref) expands it into a [`PopGrowth`](@ref) once the number of species is
known.

# Fields

  - `birth`, `death`: the rates, per unit time, the same for every species.
  - `longevity`, `survival`: the two exponents described in [`AbstractParams`](@ref).
  - `boost`: the ceiling on the birth multiplier where resources are abundant.
"""
struct EqualPop{U <: Unitful.Units} <: AbstractParams
    birth::TimeUnitType{U}
    death::TimeUnitType{U}
    longevity::Float64
    survival::Float64
    boost::Float64
end

"""
    PopGrowth{U <: Unitful.Units} <: AbstractParams

Demographic rates given per species: the general form, and what the hot loop reads.

# Fields

  - `birth`, `death`: one rate per species, per unit time.
  - `longevity`, `survival`: the two exponents described in [`AbstractParams`](@ref).
  - `boost`: the ceiling on the birth multiplier where resources are abundant.
"""
struct PopGrowth{U <: Unitful.Units} <: AbstractParams
    birth::Vector{TimeUnitType{U}}
    death::Vector{TimeUnitType{U}}
    longevity::Float64
    survival::Float64
    boost::Float64
    function PopGrowth{U}(birth::Vector{TimeUnitType{U}},
                          death::Vector{TimeUnitType{U}},
                          longevity::Float64,
                          survival::Float64,
                          boost::Float64) where {U <: Unitful.Units}
        longevity >= 0 && survival >= 0 ||
            error("longevity and survival must not be negative")
        return new{U}(birth, death, longevity, survival, boost)
    end
end

"""
    NoGrowth{U <: Unitful.Units} <: AbstractParams

Demographic rates for a population that does not grow. The same fields as [`PopGrowth`](@ref), but
the resource adjustment is short-circuited to zero, so births and deaths are taken at their stated
rates and the environment does not modulate them.

# Fields

  - `birth`, `death`: one rate per species, per unit time.
  - `longevity`, `survival`: the two exponents described in [`AbstractParams`](@ref). Carried so
    that the same parameters can be reused with a growing model, but not read while growth is off.
  - `boost`: the ceiling on the birth multiplier where resources are abundant, likewise unread.
"""
struct NoGrowth{U <: Unitful.Units} <: AbstractParams
    birth::Vector{TimeUnitType{U}}
    death::Vector{TimeUnitType{U}}
    longevity::Float64
    survival::Float64
    boost::Float64

    function NoGrowth{U}(birth::Vector{TimeUnitType{U}},
                         death::Vector{TimeUnitType{U}},
                         longevity::Float64,
                         survival::Float64,
                         boost::Float64) where {U <: Unitful.Units}
        longevity >= 0 && survival >= 0 ||
            error("longevity and survival must not be negative")
        return new{U}(birth, death, longevity, survival, boost)
    end
end

# == Functions ==================================================================================

"""
    equalpop(params::AbstractParams, numspp)

Expand [`EqualPop`](@ref)'s single pair of rates into the `numspp` per-species rates of a
[`PopGrowth`](@ref). Parameters that are already per species are returned unchanged, so a caller can
apply this without first asking which form it holds.
"""
function equalpop(params::EqualPop, numspp)
    u = unit(params.birth)
    return PopGrowth{typeof(u)}(fill(params.birth, numspp),
                                fill(params.death, numspp),
                                params.longevity,
                                params.survival,
                                params.boost)
end
function equalpop(params::PopGrowth, numspp)
    return params
end

function equalpop(params::NoGrowth, numspp)
    return params
end

# Choose the parameter type from the shape of what was passed: scalar rates give an `EqualPop`, a
# vector of either gives a `PopGrowth` of `n` species.
function _params(birth, death, longevity, survival, boost, n::Integer)
    if birth isa AbstractVector || death isa AbstractVector
        b = float.(_tofield(birth, n, "birth"))
        d = float.(_tofield(death, n, "death"))
        u = unit(first(b))
        return PopGrowth{typeof(u)}(b, d, Float64(longevity), Float64(survival),
                                    Float64(boost))
    else
        return EqualPop(float(birth), float(death), Float64(longevity),
                        Float64(survival), Float64(boost))
    end
end
