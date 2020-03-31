using Unitful
using Unitful.DefaultSymbols

"""
    TempChange(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to increase the temperature for one timestep of the ecosystem using HabitatUpdate information.
"""
function TempChange(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)
  val = hab.change.rate
  v = uconvert(K/unit(timestep), val)
  hab.matrix .+= v * timestep
end

"""
    RainfallChange(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to change the rainfall for one timestep of the ecosystem using HabitatUpdate information.
"""
function RainfallChange(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)
  val = hab.change.rate
  v = uconvert(mm/unit(timestep), val)
  hab.matrix .+= v * timestep
end

"""
    TempFluct(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to fluctuate the temperature for one timestep of the ecosystem using HabitatUpdate information.
"""
function TempFluct(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)
  val = hab.change.rate
  v = uconvert(K/unit(timestep), val) * timestep
  offset = v/pi
  hab.matrix .+= (sin.(hab.matrix ./ offset) .* v)
end

"""
    NoChange(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to keep the habitat the same for one timestep of the model.
"""
function NoChange(eco::Ecosystem, hab::AbstractHabitat, timestep::Unitful.Time)
end

ChangeLookup = Dict(K => TempChange, NoUnits => NoChange)

"""
    eraChange(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to step the ERA climate forward by one timestep.
"""
function eraChange(eco::Ecosystem, hab::ContinuousTimeHab, timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    hab.time += round(Int64, monthstep/month)
    if hab.time > size(hab.matrix, 3)
        hab.time = 1
        @warn "More timesteps than available, have repeated"
    end
end

"""
    worldclimChange(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to step the Worldclim climate forward by one timestep.
"""
function worldclimChange(eco::Ecosystem, hab::ContinuousTimeHab, timestep::Unitful.Time)
    last = size(hab.matrix, 3)
    monthstep = convert(typeof(1.0month), timestep)
    hab.time = hab.time + round(Int64,ustrip(monthstep))
    if hab.time > last
        hab.time = 1
        @warn "More timesteps than available, have repeated"
    end
end

"""
    HabitatLoss(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to destroy habitat for one timestep of the ecosystem using HabitatUpdate information.
"""
function HabitatLoss(eco::Ecosystem, hab::AbstractHabitat, timestep::Unitful.Time)
    val = hab.change.rate
    v = uconvert(unit(timestep)^-1, val)
    pos = find(eco.abenv.active)
    smp = sample(pos, jbinom(1, length(pos), ustrip(v))[1])
    eco.abenv.budget.matrix[smp] = 0.0
    eco.abundances.matrix[:, smp] .= 0.0
end

"""
    habitatupdate!(eco::Ecosystem, timestep::Unitful.Time)

Function to update the habitat of an ecosystem for one timestep.
"""
function habitatupdate!(eco::Ecosystem, timestep::Unitful.Time)
  _habitatupdate!(eco, eco.abenv.habitat, timestep)
end
function _habitatupdate!(eco::Ecosystem, hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab}, timestep::Unitful.Time)
    hab.change.changefun(eco, hab, timestep)
end

function _habitatupdate!(eco::Ecosystem, hab::HabitatCollection2, timestep::Unitful.Time)
    _habitatupdate!(eco, hab.h1, timestep)
    _habitatupdate!(eco, hab.h2, timestep)
end

"""
    budgetupdate!(eco::Ecosystem, timestep::Unitful.Time)

Function to update the budget of an ecosystem for one timestep.
"""
function budgetupdate!(eco::Ecosystem, timestep::Unitful.Time)
    _budgetupdate!(eco, eco.abenv.budget, timestep)
end
function _budgetupdate!(eco::Ecosystem, budget::SimpleBudget, timestep::Unitful.Time)
    return budget
end
function _budgetupdate!(eco::Ecosystem, budget::SolarBudget, timestep::Unitful.Time)
    return budget
end
function _budgetupdate!(eco::Ecosystem, budget::WaterBudget, timestep::Unitful.Time)
    return budget
end
function _budgetupdate!(eco::Ecosystem, budget::SolarTimeBudget, timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time +=
    round(Int64, monthstep/month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        @warn "More timesteps than available, have repeated"
    end
end
function _budgetupdate!(eco::Ecosystem, budget::WaterTimeBudget, timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time +=
    round(Int64, monthstep/month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        @warn "More timesteps than available, have repeated"
    end
end
function _budgetupdate!(eco::Ecosystem, budget::VolWaterTimeBudget, timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time +=
    round(Int64, monthstep/month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        @warn "More timesteps than available, have repeated"
    end
end

function _budgetupdate!(eco::Ecosystem, budget::BudgetCollection2{B1, B2}, timestep::Unitful.Time) where {B1, B2 <: AbstractTimeBudget}
        budget.b1.time = eco.abenv.habitat.h1.time
        budget.b2.time = eco.abenv.habitat.h1.time
end

function _budgetupdate!(eco::Ecosystem, budget::BudgetCollection2, timestep::Unitful.Time)
    _budgetupdate!(eco, eco.abenv.budget.b1, timestep)
    _budgetupdate!(eco, eco.abenv.budget.b2, timestep)
end
