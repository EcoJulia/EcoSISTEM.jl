function TempChange(eco::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)
  val = hab.change.rate
  v = uconvert(°C/unit(timestep), val)
  hab.matrix .+= v * timestep
end

function NoChange(eco::Ecosystem, hab::AbstractHabitat, timestep::Unitful.Time)
end

function HabitatLoss(eco::Ecosystem, hab::AbstractHabitat, timestep::Unitful.Time)
    val = hab.change.rate
    v = uconvert(unit(timestep)^-1, val)
    pos = find(eco.abenv.active)
    smp = sample(pos, jbinom(1, length(pos), ustrip(v))[1])
    eco.abenv.budget.matrix[smp] = 0.0
end

function habitatupdate!(eco::Ecosystem, timestep::Unitful.Time)
  _habitatupdate!(eco, eco.abenv.habitat, timestep)
end
function _habitatupdate!(eco::Ecosystem, hab::Union{DiscreteHab, ContinuousHab}, timestep::Unitful.Time)
    hab.change.changefun(eco, hab, timestep)
end
function _habitatupdate!(eco::Ecosystem, hab::Union{HabitatCollection2, HabitatCollection3}, timestep::Unitful.Time)
    habnames = fieldnames(hab)
    results = map(length(habnames)) do h
        thishab = gethabitat(hab, habnames[h])
        _habitatupdate!(eco, thishab, timestep::Unitful.Time)
    end
end
