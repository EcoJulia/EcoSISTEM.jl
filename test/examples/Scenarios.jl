using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

import Simulation: RateType, get_neighbours

"""
    TempIncrease!(eco::Ecosystem, timestep::Unitful.Time,
        rate::typeof(1.0K/year))

A function that constantly increases temperature at a certain rate.
"""
function TempIncrease!(eco::Ecosystem, timestep::Unitful.Time, rate::typeof(1.0K/year))
    resetrate!(eco, rate)
    eco.abenv.habitat.matrix[eco.abenv.habitat.matrix .< 0K] .= 0K
end

"""
    RainDecrease!(eco::Ecosystem, timestep::Unitful.Time,
        rate::typeof(1.0mm/year))

A function that constantly decreases rainfall at a certain rate.
"""
function RainDecrease!(eco::Ecosystem, timestep::Unitful.Time, rate::typeof(1.0mm/year))
    resetrate!(eco, rate)
    eco.abenv.habitat.matrix[eco.abenv.habitat.matrix .< 0mm] .= 0mm
    v = uconvert(mm/unit(timestep), rate)
    eco.abenv.budget.b2.matrix .= v * timestep
    eco.abenv.budget.b2.matrix[eco.abenv.budget.b2.matrix .< 0mm] .= 0mm
end

"""
    RandHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that randomly removes a portion of habitat at a certain rate.
"""
function RandHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
    v = uconvert(unit(timestep)^-1, rate)
    pos = findall(eco.abenv.budget.b1.matrix .> 0.0kJ)
    howmany = jbinom(1, length(pos), ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.b1.matrix[smp] .= 0.0kJ
    eco.abenv.budget.b2.matrix[smp] .= 0.0mm
    eco.abundances.grid[:, smp] .= 0.0
    eco.abenv.active[smp] .= false
end

"""
    ClustHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that removes a portion of habitat at a certain rate, which spreads through
joined pieces of habitat.
"""
function ClustHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
    v = uconvert(unit(timestep)^-1, rate)
    # If no habitat lost then pick starting places
    if all(eco.abenv.budget.b1.matrix .> 0.0kJ)
        pos = findall(eco.abenv.active)
        howmany = jbinom(1, length(pos), ustrip(v))[1]
        smp = sample(pos, howmany)
        eco.abenv.budget.b1.matrix[smp] .= 0.0kJ
        eco.abenv.budget.b2.matrix[smp] .= 0.0mm
        eco.abundances.grid[:, smp] .= 0.0
        eco.abenv.active[smp] .= false
    # Else choose from neighbours of squares already lost
    else
        pos = findall(eco.abenv.budget.b1.matrix .== 0.0kJ)
        howmany = jbinom(1, length(findall(eco.abenv.active)), ustrip(v))[1]
        width = size(eco.abenv.budget.b1.matrix,1)
        x, y = [ x[1] for x in pos ], [ x[2] for x in pos ]
        neighbours = get_neighbours(eco.abenv.budget.b1.matrix, x, y, 8)
        smp = sample(1:size(neighbours,1), howmany)
        i = convert_coords(neighbours[smp, 1], neighbours[smp, 2], width)
        eco.abenv.budget.b1.matrix[i] .= 0.0kJ
        eco.abenv.budget.b2.matrix[i] .= 0.0mm
        eco.abundances.matrix[:, i] .= 0.0
        eco.abenv.active[i] .= false
    end
    # Add in additional start points
    howmany = jbinom(1, 1, ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.b1.matrix[smp] .= 0.0kJ
    eco.abenv.budget.b2.matrix[smp] .= 0.0mm
    eco.abundances.grid[:, smp] .= 0.0
    eco.abenv.active[smp] .= false
end

"""
    GeneralistInvasive(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that introduces an invasive species into the ecosystem, `eco`, that gains abundance at each `timestep` at a particular rate, `rate`. The invasive species are introduced in one side of the grid and have the highest variance for trait preference.
"""
function GeneralistInvasive(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
    qual = eco.spplist.native .== false
    invasive = findall(qual)
    natives = findall(.!qual)
    eco.spplist.traits.var[invasive] .= maximum(eco.spplist.traits.var)
    invasive_abun = eco.spplist.abun[invasive]
    avgain = uconvert(NoUnits, rate * timestep)
    for i in eachindex(invasive)
        gains = rand(Poisson(avgain))
        pos = 1:size(eco.abenv.habitat.matrix, 1)
        smp = sample(pos, gains)
        eco.abundances.grid[invasive[i], :, 1] .+=  map(x -> sum(smp .== x), pos)
    end
end

"""
    SpecialistInvasive(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that introduces an invasive species into the ecosystem, `eco`, that gains abundance at each `timestep` at a particular rate, `rate`. The invasive species are introduced at the hottest end of the grid and take on the trait preference of that grid cell.
"""
function SpecialistInvasive(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
    qual = eco.spplist.native .== false
    invasive = findall(qual)
    natives = findall(.!qual)
    eco.spplist.traits.mean[invasive] .= mean(eco.abenv.habitat.matrix[end, :])
    eco.spplist.traits.var[invasive] .= minimum(eco.spplist.traits.var)
    invasive_abun = eco.spplist.abun[invasive]
    avgain = uconvert(NoUnits, rate * timestep)
    for i in eachindex(invasive)
        gains = rand(Poisson(avgain))
        pos = 1:size(eco.abenv.habitat.matrix, 2)
        smp = sample(pos, gains, replace = true)
        eco.abundances.grid[invasive[i], end, :] .+= map(x -> sum(smp .== x), pos)
    end
end

function TempFluct!(eco::Ecosystem, timestep::Unitful.Time, rate::Quantity{Float64, 𝚯*𝐓^-1}, currentstep::Unitful.Time, startarray::Array{Unitful.Temperature{Float64}, 2})
    v = uconvert(K/year, rate)
    eco.abenv.habitat.matrix .= (v * year) .* sin(collect(-π:π/6:π)[mod(Int64(uconvert(NoUnits, currentstep/timestep)),12) + 1]) .+ startarray
end

function runscenario!(eco::Ecosystem, timestep::Unitful.Time, scenario::FluctScenario, currentstep::Unitful.Time)
    scenario.fun(eco, timestep, scenario.rate, currentstep, scenario.startarray)
end
function runscenario!(eco::Ecosystem, timestep::Unitful.Time, scenario::SimpleScenario, currentstep::Unitful.Time)
    scenario.fun(eco, timestep, scenario.rate)
end

function runscenario!(eco::Ecosystem, timestep::Unitful.Time, scenario::MultiScenario, currentstep::Unitful.Time)
    scenario.sc1.fun(eco, timestep, scenario.sc1.rate)
    scenario.sc2.fun(eco, timestep, scenario.sc2.rate, currentstep, scenario.sc2.startarray)
end
