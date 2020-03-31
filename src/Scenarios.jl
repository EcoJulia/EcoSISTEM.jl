
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

"""
    AbstractScenario

Abstract supertype for all whole ecosystem change scenarios
"""
abstract type AbstractScenario end

RateType = typeof(1.0/year)

"""
    SimpleScenario <: AbstractScenario

This scenario type holds a function that acts to change the entire ecosystem.
"""
mutable struct SimpleScenario <: AbstractScenario
    fun::Function
    rate::Union{Quantity{Float64, 𝐓^-1}, Quantity{Float64, 𝚯*𝐓^-1}, Quantity{Float64, 𝐋*𝐓^-1}}
end

mutable struct FluctScenario <: AbstractScenario
        fun::Function
        rate::Quantity{Float64, 𝚯*𝐓^-1}
        startarray::Array{Unitful.Temperature{Float64}, 2}
end

mutable struct MultiScenario{S1 <: AbstractScenario, S2 <: AbstractScenario} <: AbstractScenario
        sc1::S1
        sc2::S2
end


function runscenario!(eco::Ecosystem, timestep::Unitful.Time, scenario::SimpleScenario, currentstep::Unitful.Time)
    scenario.fun(eco, timestep, scenario.rate)
end

function runscenario!(eco::Ecosystem, timestep::Unitful.Time, scenario::FluctScenario, currentstep::Unitful.Time)
    scenario.fun(eco, timestep, scenario.rate, currentstep, scenario.startarray)
end

function runscenario!(eco::Ecosystem, timestep::Unitful.Time, scenario::MultiScenario, currentstep::Unitful.Time)
    scenario.sc1.fun(eco, timestep, scenario.sc1.rate)
    scenario.sc2.fun(eco, timestep, scenario.sc2.rate, currentstep, scenario.sc2.startarray)
end
