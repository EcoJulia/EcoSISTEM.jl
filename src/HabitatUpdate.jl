# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols

"""
    TempChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to increase the temperature for one timestep of the ecosystem using HabitatUpdate information.
"""
function TempChange(eco::AbstractEcosystem, hab::ContinuousHab,
                    timestep::Unitful.Time)
    val = hab.change.rate
    v = uconvert(K / unit(timestep), val)
    return hab.matrix .+= v * timestep
end

"""
    RainfallChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to change the rainfall for one timestep of the ecosystem using HabitatUpdate information.
"""
function RainfallChange(eco::AbstractEcosystem, hab::ContinuousHab,
                        timestep::Unitful.Time)
    val = hab.change.rate
    v = uconvert(mm / unit(timestep), val)
    return hab.matrix .+= v * timestep
end

"""
    TempFluct(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to fluctuate the temperature for one timestep of the ecosystem using HabitatUpdate information.
"""
function TempFluct(eco::AbstractEcosystem, hab::ContinuousHab,
                   timestep::Unitful.Time)
    val = hab.change.rate
    v = uconvert(K / unit(timestep), val) * timestep
    offset = v / pi
    return hab.matrix .+= (sin.(hab.matrix ./ offset) .* v)
end

"""
    NoChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to keep the habitat the same for one timestep of the model.
"""
function NoChange(eco::AbstractEcosystem, hab::AbstractHabitat,
                  timestep::Unitful.Time)
end

ChangeLookup = Dict(K => TempChange, NoUnits => NoChange)

"""
    eraChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to step the ERA climate forward by one timestep.
"""
function eraChange(eco::AbstractEcosystem, hab::ContinuousTimeHab,
                   timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    hab.time += round(Int64, monthstep / month)
    if hab.time > size(hab.matrix, 3)
        hab.time = 1
        @warn "More timesteps than available, have repeated"
    end
end

"""
    worldclimChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to step the Worldclim climate forward by one timestep.
"""
function worldclimChange(eco::AbstractEcosystem, hab::ContinuousTimeHab,
                         timestep::Unitful.Time)
    last = size(hab.matrix, 3)
    monthstep = convert(typeof(1.0month), timestep)
    hab.time = hab.time + round(Int64, ustrip(monthstep))
    if hab.time > last
        hab.time = 1
        @warn "More timesteps than available, have repeated"
    end
end

"""
    HabitatLoss(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to destroy habitat for one timestep of the ecosystem using HabitatUpdate information.
"""
function HabitatLoss(eco::AbstractEcosystem, hab::AbstractHabitat,
                     timestep::Unitful.Time)
    val = hab.change.rate
    v = uconvert(unit(timestep)^-1, val)
    pos = find(eco.abenv.active)
    smp = sample(pos, jbinom(1, length(pos), ustrip(v))[1])
    eco.abenv.budget.matrix[smp] = 0.0
    return eco.abundances.matrix[:, smp] .= 0.0
end

"""
    habitatupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)

Function to update the habitat of an ecosystem for one timestep.
"""
function habitatupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return _habitatupdate!(eco, eco.abenv.habitat, timestep)
end
function _habitatupdate!(eco::AbstractEcosystem,
                         hab::Union{DiscreteHab, ContinuousHab,
                                    ContinuousTimeHab}, timestep::Unitful.Time)
    return hab.change.changefun(eco, hab, timestep)
end

function _habitatupdate!(eco::AbstractEcosystem, hab::HabitatCollection2,
                         timestep::Unitful.Time)
    _habitatupdate!(eco, hab.h1, timestep)
    return _habitatupdate!(eco, hab.h2, timestep)
end

"""
    budgetupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)

Function to update the budget of an ecosystem for one timestep.
"""
function budgetupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return _budgetupdate!(eco, eco.abenv.budget, timestep)
end
function _budgetupdate!(eco::AbstractEcosystem, budget::SimpleBudget,
                        timestep::Unitful.Time)
    return budget
end
function _budgetupdate!(eco::AbstractEcosystem, budget::SolarBudget,
                        timestep::Unitful.Time)
    return budget
end
function _budgetupdate!(eco::AbstractEcosystem, budget::WaterBudget,
                        timestep::Unitful.Time)
    return budget
end
function _budgetupdate!(eco::AbstractEcosystem, budget::SolarTimeBudget,
                        timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time += round(Int64, monthstep / month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        @warn "More timesteps than available, have repeated"
    end
end
function _budgetupdate!(eco::AbstractEcosystem, budget::WaterTimeBudget,
                        timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time += round(Int64, monthstep / month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        @warn "More timesteps than available, have repeated"
    end
end
function _budgetupdate!(eco::AbstractEcosystem, budget::VolWaterTimeBudget,
                        timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time += round(Int64, monthstep / month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        @warn "More timesteps than available, have repeated"
    end
end

function _budgetupdate!(eco::AbstractEcosystem,
                        budget::BudgetCollection2{B1, B2},
                        timestep::Unitful.Time) where {B1,
                                                       B2 <: AbstractTimeBudget}
    budget.b1.time = eco.abenv.habitat.h1.time
    return budget.b2.time = eco.abenv.habitat.h1.time
end

function _budgetupdate!(eco::AbstractEcosystem, budget::BudgetCollection2,
                        timestep::Unitful.Time)
    _budgetupdate!(eco, eco.abenv.budget.b1, timestep)
    return _budgetupdate!(eco, eco.abenv.budget.b2, timestep)
end
