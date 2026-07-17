# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols

"""
    TempChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Increase the temperature for one timestep of the ecosystem using
[`HabitatUpdate`](@ref) information.
"""
function TempChange(eco::AbstractEcosystem, hab::ContinuousHab,
                    timestep::Unitful.Time)
    val = hab.change.rate
    v = uconvert(K / unit(timestep), val)
    return hab.matrix .+= v * timestep
end

"""
    RainfallChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Change the rainfall for one timestep of the ecosystem using
[`HabitatUpdate`](@ref) information.
"""
function RainfallChange(eco::AbstractEcosystem, hab::ContinuousHab,
                        timestep::Unitful.Time)
    val = hab.change.rate
    v = uconvert(mm / unit(timestep), val)
    return hab.matrix .+= v * timestep
end

"""
    TempFluct(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Fluctuate the temperature for one timestep of the ecosystem using
[`HabitatUpdate`](@ref) information.
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

Keep the habitat the same for one timestep of the model.
"""
function NoChange(eco::AbstractEcosystem, hab::AbstractHabitat,
                  timestep::Unitful.Time) end

ChangeLookup = Dict(K => TempChange, NoUnits => NoChange)

"""
    cyclicChange(eco::AbstractEcosystem, hab::ContinuousTimeHab, timestep::Unitful.Time)

Advance a time-varying (monthly) climate layer by one timestep, wrapping back to the
first month once the stored time series is exhausted (with a warning). The per-source
`eraChange`/`worldclimChange` names are aliases of this — the time-step-and-wrap logic is
identical regardless of source.
"""
function cyclicChange(eco::AbstractEcosystem, hab::ContinuousTimeHab,
                      timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    hab.time += round(Int64, monthstep / month)
    if hab.time > size(hab.matrix, 3)
        hab.time = 1
        @warn "More timesteps than available, have repeated"
    end
end

const eraChange = cyclicChange
const worldclimChange = cyclicChange

"""
    HabitatLoss(eco::AbstractEcosystem, hab::AbstractHabitat, timestep::Unitful.Time)

Destroy habitat for one timestep of the ecosystem using [`HabitatUpdate`](@ref)
information. The habitat's `change.rate` is a loss rate (per unit time); over
`timestep` it gives the per-cell probability that an active cell is lost. That many
active cells are drawn at random and have their budget and abundances zeroed.
"""
function HabitatLoss(eco::AbstractEcosystem, hab::AbstractHabitat,
                     timestep::Unitful.Time)
    # Loss rate × timestep is the (dimensionless) probability of losing a cell.
    prob = uconvert(NoUnits, hab.change.rate * timestep)
    pos = findall(vec(eco.abenv.active))
    smp = sample(pos, rand(Binomial(length(pos), prob)); replace = false)
    eco.abenv.budget.matrix[smp] .= zero(eltype(eco.abenv.budget.matrix))
    eco.abundances.matrix[:, smp] .= 0
    return eco
end

"""
    habitatupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)

Update the habitat of an ecosystem for one timestep.
"""
function habitatupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return _habitatupdate!(eco, eco.abenv.habitat, timestep)
end
function _habitatupdate!(eco::AbstractEcosystem,
                         hab::Union{DiscreteHab, ContinuousHab,
                                    ContinuousTimeHab},
                         timestep::Unitful.Time)
    return hab.change.changefun(eco, hab, timestep)
end

function _habitatupdate!(eco::AbstractEcosystem,
                         hab::HabitatCollection2,
                         timestep::Unitful.Time)
    _habitatupdate!(eco, hab.h1, timestep)
    return _habitatupdate!(eco, hab.h2, timestep)
end

"""
    budgetupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)

Update the budget of an ecosystem for one timestep.
"""
function budgetupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return _budgetupdate!(eco, eco.abenv.budget, timestep)
end
function _budgetupdate!(eco::AbstractEcosystem,
                        budget::SimpleBudget,
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
function _budgetupdate!(eco::AbstractEcosystem,
                        budget::SolarTimeBudget,
                        timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time += round(Int64, monthstep / month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        @warn "More timesteps than available, have repeated"
    end
end
function _budgetupdate!(eco::AbstractEcosystem,
                        budget::WaterTimeBudget,
                        timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time += round(Int64, monthstep / month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        @warn "More timesteps than available, have repeated"
    end
end
function _budgetupdate!(eco::AbstractEcosystem,
                        budget::VolWaterTimeBudget,
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

function _budgetupdate!(eco::AbstractEcosystem,
                        budget::BudgetCollection2,
                        timestep::Unitful.Time)
    _budgetupdate!(eco, eco.abenv.budget.b1, timestep)
    return _budgetupdate!(eco, eco.abenv.budget.b2, timestep)
end
