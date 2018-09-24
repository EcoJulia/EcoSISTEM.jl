function TempChange(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)
  val = hab.change.rate
  v = uconvert(°C/unit(timestep), val)
  hab.matrix .+= v * timestep
end

function NoChange(eco::Ecosystem, hab::AbstractHabitat, timestep::Unitful.Time)
end
function eraChange(eco::Ecosystem, hab::ContinuousTimeHab, timestep::Unitful.Time)
    last = size(hab.matrix, 3)
    monthstep = convert(typeof(1.0month), timestep)
    hab.time = hab.time + round(Int64,ustrip(monthstep))
    if hab.time > last
        hab.time = 1
        warn("More timesteps than available, have repeated")
    end
end
function worldclimChange(eco::Ecosystem, hab::ContinuousTimeHab, timestep::Unitful.Time)
    last = size(hab.matrix, 3)
    monthstep = convert(typeof(1.0month), timestep)
    hab.time = hab.time + round(Int64,ustrip(monthstep))
    if hab.time > last
        hab.time = 1
        warn("More timesteps than available, have repeated")
    end
end

function HabitatLoss(eco::Ecosystem, hab::AbstractHabitat, timestep::Unitful.Time)
    val = hab.change.rate
    v = uconvert(unit(timestep)^-1, val)
    pos = find(eco.abenv.active)
    smp = sample(pos, jbinom(1, length(pos), ustrip(v))[1])
    eco.abenv.budget.matrix[smp] = 0.0
    eco.abundances.matrix[:, smp] .= 0.0
end

function habitatupdate!(eco::Ecosystem, timestep::Unitful.Time)
  _habitatupdate!(eco, eco.abenv.habitat, timestep)
end
function _habitatupdate!(eco::Ecosystem, hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab}, timestep::Unitful.Time)
    hab.change.changefun(eco, hab, timestep)
end
function _habitatupdate!(eco::Ecosystem, hab::Union{HabitatCollection2, HabitatCollection3}, timestep::Unitful.Time)
    habnames = fieldnames(hab)
    results = map(length(habnames)) do h
        thishab = gethabitat(hab, habnames[h])
        _habitatupdate!(eco, thishab, timestep::Unitful.Time)
    end
end

function budgetupdate!(eco::Ecosystem, timestep::Unitful.Time)
    _budgetupdate!(eco, eco.abenv.budget, timestep)
end
function _budgetupdate!(eco::Ecosystem, budget::SimpleBudget, timestep::Unitful.Time)
    return budget
end
function _budgetupdate!(eco::Ecosystem, budget::SolarBudget, timestep::Unitful.Time)
    lastE = size(budget.matrix, 3)
    monthstep = convert(typeof(1.0month), timestep)
    eco.abenv.budget.time = eco.abenv.budget.time +
    round(Int64,ustrip(monthstep))
    if eco.abenv.budget.time > lastE
        eco.abenv.budget.time = 1
        warn("More timesteps than available, have repeated")
    end
end
function _budgetupdate!(eco::Ecosystem, budget::WaterBudget, timestep::Unitful.Time)
    lastE = size(budget.matrix, 3)
    monthstep = convert(typeof(1.0month), timestep)
    eco.abenv.budget.time = eco.abenv.budget.time +
    round(Int64,ustrip(monthstep))
    if eco.abenv.budget.time > lastE
        eco.abenv.budget.time = 1
        warn("More timesteps than available, have repeated")
    end
end
function _budgetupdate!(eco::Ecosystem, budget::BudgetCollection2, timestep::Unitful.Time)
    for i in fieldnames(budget)
        bud = getfield(budget, i)
        lastE = size(bud.matrix, 3)
        monthstep = convert(typeof(1.0month), timestep)
        getfield(eco.abenv.budget, i).time = getfield(eco.abenv.budget, i).time +
        round(Int64,ustrip(monthstep))
        if getfield(eco.abenv.budget, i).time > lastE
            getfield(eco.abenv.budget, i).time = 1
            warn("More timesteps than available, have repeated")
        end
    end
end
