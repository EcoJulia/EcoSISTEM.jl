
using Unitful
using Unitful.DefaultSymbols

"""
    AbstractScenario

Abstract supertype for all whole ecosystem change scenarios
"""
abstract type AbstractScenario end

"""
    SimpleScenario <: AbstractScenario

This scenario type holds a function that acts to change the entire ecosystem.
"""
mutable struct SimpleScenario <: AbstractScenario
    fun::Function
    rate::Quantity{Float64, typeof(𝐓^-1)}
end

function RandHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time, rate::Quantity{Float64, typeof(𝐓^-1)})
    v = uconvert(unit(timestep)^-1, rate)
    pos = find(eco.abenv.active)
    smp = sample(pos, jbinom(1, length(pos), ustrip(v))[1])
    eco.abenv.budget.matrix[smp] = 0.0
end

function ClustHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time, rate::Quantity{Float64, typeof(𝐓^-1)})
    v = uconvert(unit(timestep)^-1, rate)
    if all(eco.abenv.budget.matrix .> 0.0)
        pos = find(eco.abenv.active)
        smp = sample(pos, jbinom(1, length(pos), ustrip(v))[1])
        eco.abenv.budget.matrix[smp] = 0.0
    else
        pos = find(eco.abenv.budget.matrix .== 0.0)
        width = size(eco.abenv.budget.matrix,1)
        x, y = convert_coords(pos, width)
        neighbours = get_neighbours(eco.abenv.budget.matrix, x, y, 8)
        smp = sample(1:size(neighbours,1), jbinom(1, size(neighbours,1), ustrip(v))[1])
        i = convert_coords(x, y, width)
        eco.abenv.budget.matrix[i]=0.0
    end
end
