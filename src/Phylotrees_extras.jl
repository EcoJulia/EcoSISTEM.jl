abstract TreeComponent

"""
A parametric node of phylogenetic tree
"""
type Node{N} <: TreeComponent
  label::Nullable{String}
  in::Vector{Int64}
  out::Vector{Int64}
  data::Nullable{N}

  function Node()
    new(Nullable{String}(), Int64[], Int64[], Nullable())
  end

  function Node(label::String)
    new(Nullable(label), Int64[], Int64[], Nullable())
  end
end


function setdata!(x::TreeComponent, data)
  x.data = Nullable(data)
  return x
end
function hasdata(x::TreeComponent)
  return !isnull(x.data)
end
