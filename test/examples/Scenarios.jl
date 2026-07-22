# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units

import EcoSISTEM: RateType, get_neighbours

"""
    TempIncrease!(eco::Ecosystem, timestep::Unitful.Time,
        rate::typeof(1.0K/year))

A function that constantly increases temperature at a certain rate.
"""
function TempIncrease!(eco::Ecosystem, timestep::Unitful.Time,
                       rate::typeof(1.0K / year))
    resetrate!(eco, rate)
    return eco.habitat.regime.matrix[eco.habitat.regime.matrix .< 0K] .= 0K
end

"""
    RainDecrease!(eco::Ecosystem, timestep::Unitful.Time,
        rate::typeof(1.0mm/year))

A function that constantly decreases rainfall at a certain rate.
"""
function RainDecrease!(eco::Ecosystem, timestep::Unitful.Time,
                       rate::typeof(1.0mm / year))
    resetrate!(eco, rate)
    eco.habitat.regime.matrix[eco.habitat.regime.matrix .< 0mm] .= 0mm
    v = uconvert(mm / unit(timestep), rate)
    eco.habitat.supply.two.matrix .= v * timestep
    return eco.habitat.supply.two.matrix[eco.habitat.supply.two.matrix .< 0mm] .= 0mm
end

"""
    RandHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that randomly removes a portion of regime at a certain rate.
"""
function RandHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time,
                          rate::RateType)
    v = uconvert(unit(timestep)^-1, rate)
    pos = findall(eco.habitat.supply.one.matrix .> 0.0kJ)
    howmany = jbinom(1, length(pos), ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.habitat.supply.one.matrix[smp] .= 0.0kJ
    eco.habitat.supply.two.matrix[smp] .= 0.0mm
    eco.abundances.grid[:, smp] .= 0.0
    return eco.habitat.active[smp] .= false
end

"""
    ClustHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that removes a portion of regime at a certain rate, which spreads through
joined pieces of regime.
"""
function ClustHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time,
                           rate::RateType)
    v = uconvert(unit(timestep)^-1, rate)
    # If no regime lost then pick starting places
    if all(eco.habitat.supply.one.matrix .> 0.0kJ)
        pos = findall(eco.habitat.active)
        howmany = jbinom(1, length(pos), ustrip(v))[1]
        smp = sample(pos, howmany)
        eco.habitat.supply.one.matrix[smp] .= 0.0kJ
        eco.habitat.supply.two.matrix[smp] .= 0.0mm
        eco.abundances.grid[:, smp] .= 0.0
        eco.habitat.active[smp] .= false
        # Else choose from neighbours of squares already lost
    else
        pos = findall(eco.habitat.supply.one.matrix .== 0.0kJ)
        howmany = jbinom(1, length(findall(eco.habitat.active)), ustrip(v))[1]
        width = size(eco.habitat.supply.one.matrix, 1)
        x, y = [x[1] for x in pos], [x[2] for x in pos]
        neighbours = get_neighbours(eco.habitat.supply.one.matrix, x, y, 8)
        smp = sample(Base.axes(neighbours, 1), howmany)
        i = convert_coords(neighbours[smp, 1], neighbours[smp, 2], width)
        eco.habitat.supply.one.matrix[i] .= 0.0kJ
        eco.habitat.supply.two.matrix[i] .= 0.0mm
        eco.abundances.matrix[:, i] .= 0.0
        eco.habitat.active[i] .= false
    end
    # Add in additional start points
    howmany = jbinom(1, 1, ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.habitat.supply.one.matrix[smp] .= 0.0kJ
    eco.habitat.supply.two.matrix[smp] .= 0.0mm
    eco.abundances.grid[:, smp] .= 0.0
    return eco.habitat.active[smp] .= false
end

"""
    GeneralistInvasive(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that introduces an invasive species into the ecosystem, `eco`, that gains abundance at each `timestep` at a particular rate, `rate`. The invasive species are introduced in one side of the grid and have the highest variance for trait preference.
"""
function GeneralistInvasive(eco::Ecosystem, timestep::Unitful.Time,
                            rate::RateType)
    qual = eco.spplist.native .== false
    invasive = findall(qual)
    natives = findall(.!qual)
    eco.spplist.tolerance.var[invasive] .= maximum(eco.spplist.tolerance.var)
    invasive_abun = eco.spplist.abun[invasive]
    avgain = uconvert(NoUnits, rate * timestep)
    for i in eachindex(invasive)
        gains = rand(Poisson(avgain))
        pos = Base.axes(eco.habitat.regime.matrix, 1)
        smp = sample(pos, gains)
        eco.abundances.grid[invasive[i], :, 1] .+= map(x -> sum(smp .== x), pos)
    end
end

"""
    SpecialistInvasive(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that introduces an invasive species into the ecosystem, `eco`, that gains abundance at each `timestep` at a particular rate, `rate`. The invasive species are introduced at the hottest end of the grid and take on the trait preference of that grid cell.
"""
function SpecialistInvasive(eco::Ecosystem, timestep::Unitful.Time,
                            rate::RateType)
    qual = eco.spplist.native .== false
    invasive = findall(qual)
    natives = findall(.!qual)
    eco.spplist.tolerance.mean[invasive] .= mean(eco.habitat.regime.matrix[end,
                                                                           :])
    eco.spplist.tolerance.var[invasive] .= minimum(eco.spplist.tolerance.var)
    invasive_abun = eco.spplist.abun[invasive]
    avgain = uconvert(NoUnits, rate * timestep)
    for i in eachindex(invasive)
        gains = rand(Poisson(avgain))
        pos = Base.axes(eco.habitat.regime.matrix, 2)
        smp = sample(pos, gains, replace = true)
        eco.abundances.grid[invasive[i], end, :] .+= map(x -> sum(smp .== x),
                                                         pos)
    end
end

function TempFluct!(eco::Ecosystem,
                    timestep::Unitful.Time,
                    rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                    currentstep::Unitful.Time,
                    startarray::Matrix{Unitful.Temperature{Float64}})
    v = uconvert(K / year, rate)
    return eco.habitat.regime.matrix .= (v * year) .*
                                        sin(collect((-π):(π / 6):π)[mod(Int64(uconvert(NoUnits,
                                                                                       currentstep / timestep)),
                                                                        12) + 1]) .+
                                        startarray
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
                      scenario::SimpleScenario,
                      currentstep::Unitful.Time)
    return scenario.fun(eco, timestep, scenario.rate)
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
