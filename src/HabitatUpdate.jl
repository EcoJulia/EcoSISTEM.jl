using Compat
using Unitful
using Unitful.DefaultSymbols

function TempChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)
  val = hab.change.rate
  v = uconvert(K/unit(timestep), val)
  hab.matrix .+= v * timestep
end

function RainfallChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)
  val = hab.change.rate
  v = uconvert(mm/unit(timestep), val)
  hab.matrix .+= v * timestep
end

function TempFluct(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)
  val = hab.change.rate
  v = uconvert(K/unit(timestep), val) * timestep
  offset = v/pi
  hab.matrix .+= (sin.(hab.matrix ./ offset) .* v)
end

GLOBAL_funcdict["TempChange"] = TempChange
function NoChange(eco::AbstractEcosystem, hab::AbstractHabitat, timestep::Unitful.Time)
end
GLOBAL_funcdict["NoChange"] = NoChange
ChangeLookup = Dict(K => TempChange, NoUnits => NoChange)
function eraChange(eco::AbstractEcosystem, hab::ContinuousTimeHab, timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    hab.time += round(Int64, monthstep/month)
    if hab.time > size(hab.matrix, 3)
        hab.time = 1
        Compat.@warn "More timesteps than available, have repeated"
    end
end
GLOBAL_funcdict["eraChange"] = eraChange
function worldclimChange(eco::AbstractEcosystem, hab::ContinuousTimeHab, timestep::Unitful.Time)
    last = size(hab.matrix, 3)
    monthstep = convert(typeof(1.0month), timestep)
    hab.time = hab.time + round(Int64,ustrip(monthstep))
    if hab.time > last
        hab.time = 1
        Compat.@warn "More timesteps than available, have repeated"
    end
end
GLOBAL_funcdict["worldclimChange"] = worldclimChange
function HabitatLoss(eco::AbstractEcosystem, hab::AbstractHabitat, timestep::Unitful.Time)
    val = hab.change.rate
    v = uconvert(unit(timestep)^-1, val)
    pos = find(eco.abenv.active)
    smp = sample(pos, jbinom(1, length(pos), ustrip(v))[1])
    eco.abenv.budget.matrix[smp] = 0.0
    eco.abundances.matrix[:, smp] .= 0.0
end
GLOBAL_funcdict["HabitatLoss"] = HabitatLoss
function habitatupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
  _habitatupdate!(eco, eco.abenv.habitat, timestep)
end
function _habitatupdate!(eco::AbstractEcosystem, hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab}, timestep::Unitful.Time)
    hab.change.changefun(eco, hab, timestep)
end

function _habitatupdate!(eco::AbstractEcosystem, hab::HabitatCollection2, timestep::Unitful.Time)
    _habitatupdate!(eco, hab.h1, timestep)
    _habitatupdate!(eco, hab.h2, timestep)
end

function budgetupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    _budgetupdate!(eco, eco.abenv.budget, timestep)
end
function _budgetupdate!(eco::AbstractEcosystem, budget::SimpleBudget, timestep::Unitful.Time)
    return budget
end
function _budgetupdate!(eco::AbstractEcosystem, budget::SolarBudget, timestep::Unitful.Time)
    return budget
end
function _budgetupdate!(eco::AbstractEcosystem, budget::WaterBudget, timestep::Unitful.Time)
    return budget
end
function _budgetupdate!(eco::AbstractEcosystem, budget::SolarTimeBudget, timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time +=
    round(Int64, monthstep/month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        Compat.@warn "More timesteps than available, have repeated"
    end
end
function _budgetupdate!(eco::AbstractEcosystem, budget::WaterTimeBudget, timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time +=
    round(Int64, monthstep/month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        Compat.@warn "More timesteps than available, have repeated"
    end
end
function _budgetupdate!(eco::AbstractEcosystem, budget::VolWaterTimeBudget, timestep::Unitful.Time)
    monthstep = uconvert(month, timestep)
    budget.time +=
    round(Int64, monthstep/month)
    if budget.time > size(budget.matrix, 3)
        budget.time = 1
        Compat.@warn "More timesteps than available, have repeated"
    end
end

function _budgetupdate!(eco::AbstractEcosystem, budget::BudgetCollection2{B1, B2}, timestep::Unitful.Time) where {B1, B2 <: AbstractTimeBudget}
        budget.b1.time = eco.abenv.habitat.h1.time
        budget.b2.time = eco.abenv.habitat.h1.time
end

function _budgetupdate!(eco::AbstractEcosystem, budget::BudgetCollection2, timestep::Unitful.Time)
    _budgetupdate!(eco, eco.abenv.budget.b1, timestep)
    _budgetupdate!(eco, eco.abenv.budget.b2, timestep)
end
