

function setdata!(x::TreeComponent, data)
  x.data = Nullable(data)
  return x
end
function hasdata(x::TreeComponent)
  return !isnull(x.data)
end
