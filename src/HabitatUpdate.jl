function TempChange(eco::Ecosystem, timestep::Unitful.Time)
  val = gethabitat(eco).change.rate
  v = uconvert(°C/unit(timestep), val)
  gethabitat(eco).matrix .+= v
end

function NoChange(eco::Ecosystem, timestep::Unitful.Time)
end

function HabitatLoss(eco::Ecosystem, timestep::Unitful.Time)
    val = gethabitat(eco).change.rate
    v = uconvert(unit(timestep)^-1, val)
    pos = find(eco.abenv.active)
    smp = sample(pos, jbinom(1, length(pos), ustrip(v))[1])
    eco.abenv.budget.matrix[smp] = 0.0
end

function getchangefun(eco::Ecosystem)
  return gethabitat(eco).change.changefun
end
