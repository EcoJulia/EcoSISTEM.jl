function TempChange(eco::Ecosystem, timestep::Unitful.Time)
  val = gethabitat(eco).change.rate
  v = uconvert(°C/unit(timestep), val)
  gethabitat(eco).matrix .+= v
end

function NoChange(eco::Ecosystem, timestep::Unitful.Time)
end

function getchangefun(eco::Ecosystem)
  return gethabitat(eco).change.changefun
end
