# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units

"""
    AbstractScenario

Abstract supertype for all whole ecosystem change scenarios
"""
abstract type AbstractScenario end

RateType = typeof(1.0 / year)

"""
    SimpleScenario <: AbstractScenario

Scenario type that applies a uniform rate of change across the entire ecosystem
at each timestep. `fun` is the change function called as `fun(eco, timestep,
rate)`, and `rate` is the magnitude of change per unit time.
"""
struct SimpleScenario{F <: Function} <: AbstractScenario
    fun::F
    rate::Union{Quantity{Float64, 𝐓^-1},
                Quantity{Float64, 𝚯 * 𝐓^-1},
                Quantity{Float64, 𝐋 * 𝐓^-1}}
end

"""
    FluctScenario <: AbstractScenario

Scenario type that fluctuates the environment periodically. `fun` is the
fluctuation function, `rate` is the rate of temperature change, and `startarray`
holds the baseline regime values at the start of the simulation.
"""
struct FluctScenario{F <: Function} <: AbstractScenario
    fun::F
    rate::Quantity{Float64, 𝚯 * 𝐓^-1}
    startarray::Matrix{typeof(1.0K)}
end

"""
    MultiScenario{S1 <: AbstractScenario, S2 <: AbstractScenario} <: AbstractScenario

Scenario type that composes two scenarios, applying `sc1` and then `sc2` in
sequence at each timestep.
"""
struct MultiScenario{S1 <: AbstractScenario, S2 <: AbstractScenario} <:
       AbstractScenario
    sc1::S1
    sc2::S2
end

"""
    runscenario!(eco::Ecosystem, timestep::Unitful.Time, scenario::S, currentstep::Unitful.Time) where S <: AbstractScenario

This function runs any scenario type for one timestep.
"""
function runscenario!(eco::Ecosystem,
                      timestep::Unitful.Time,
                      scenario::SimpleScenario,
                      currentstep::Unitful.Time)
    return scenario.fun(eco, timestep, scenario.rate)
end

function runscenario!(eco::Ecosystem,
                      timestep::Unitful.Time,
                      scenario::FluctScenario,
                      currentstep::Unitful.Time)
    return scenario.fun(eco, timestep, scenario.rate, currentstep,
                        scenario.startarray)
end

function runscenario!(eco::Ecosystem,
                      timestep::Unitful.Time,
                      scenario::MultiScenario,
                      currentstep::Unitful.Time)
    scenario.sc1.fun(eco, timestep, scenario.sc1.rate)
    return scenario.sc2.fun(eco,
                            timestep,
                            scenario.sc2.rate,
                            currentstep,
                            scenario.sc2.startarray)
end
