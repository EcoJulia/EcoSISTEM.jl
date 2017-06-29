function TempChange(eco::Ecosystem)
  val = gethabitat(eco).change.rate
  gethabitat(eco).matrix .+= val
end

function getchangefun(eco::Ecosystem)
  return gethabitat(eco).change.changefun
end
