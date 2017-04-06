using RCall
function plot_move(eco::Ecosystem, x::Int64, y::Int64, spp::Int64)
  table = eco.lookup[spp] .+ [x y 0]
  maxGrid = maximum(size(eco.abenv.habitat.matrix))
  # Can't go over maximum dimension
  lower  = find(mapslices(x->all(x.>0), table, 2))
  upper = find(mapslices(x->all(x.<= maxGrid), table, 2))
  valid = intersect(lower, upper)
  table = table[valid, :]
  table[:, 3] = table[:, 3]/sum(table[:, 3])
  table
  A = zeros(size(eco.abenv.habitat.matrix))
  for i in eachindex(table[:, 1])
    A[table[i, 1], table[i, 2]] = table[i, 3]
  end
  @rput A
  R"library(fields);
  A[A==0]=NA
  image.plot(A)"
end
