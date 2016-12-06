tree=jcoal(3,5)
function treeplot(tree::Tree)
  nodequeue = findroots(tree)
  treesize = descendantcount(tree, nodequeue) + 1
  distances = distance(tree, nodequeue)
  countpct = treesize / sum(treesize)
  height = cumsum(countpct) .- (0.5 * countpct)
  queueposition = 1
  while(queueposition <= length(nodequeue))
    children = childnodes(tree, nodequeue[queueposition])
    if length(children) > 0
      append!(nodequeue, children)
      subtreesize = descendantcount(tree, children) + 1
      append!(distances, distance(tree, children))
      append!(countpct, subtreesize / sum(treesize))
      append!(height, (height[queueposition] - (countpct[queueposition] / 2)) + (cumsum(countpct[end-length(subtreesize)+1:end]) .- (0.5 * countpct[end-length(subtreesize)+1:end])))
    end
    queueposition += 1
  end
  processorder = fill(Nullable{Int64}(), length(tree.nodes))
  for i = 1:length(nodequeue)
    processorder[nodequeue[i]] = i
  end
  tree_x = Vector{Float64}[]
  tree_y = Vector{Float64}[]
  xmax = Float64[]
  for i in nodequeue
    if !isroot(tree, i)
      push!(tree_x, distances[[get(processorder[i]), get(processorder[parentnode(tree, i)]), get(processorder[parentnode(tree, i)])]])
      push!(tree_y, height[[get(processorder[i]), get(processorder[i]), get(processorder[parentnode(tree, i)])]])
      push!(xmax, distances[get(processorder[i])])
    end
  end
  return tree_x, tree_y, maximum(xmax)
end

@recipe function plot(tree::Tree)
  tree_x, tree_y, xmax = treeplot(tree)
  seriestype := :path
  linecolor --> :black
  legend := false
  yticks --> nothing
  xlims --> (-1., xmax+1.)
  ylims --> (0., 1.)
  tree_x, tree_y
end


plot(tree)
trt=get_traits(trait_tree)
tree=trait_tree

tree_x, tree_y=treeplot(tree)[1],treeplot(tree)[2]

num_traits=replace(trt, "A", "red")
trt=map(i->replace(trt[i], "A", "red"),1:5)
trt=map(i->replace(trt[i], "B", "blue"),1:5)
for i in 1:length(tree_x)
  
plot!(tree_x[i], tree_y[i],
markershape=:circle,
markercolor=[:blue,false, :red],
markerstrokecolor=[:black, false, :black],
linecolor=:black)

plot!(tree_x[2], tree_y[2],
markershape=:circle,
markercolor=[:red, false, :red],
markerstrokecolor=[:black, false, :black],
linecolor=:black)

plot!(tree_x[3], tree_y[3],
markershape=:circle,
markercolor=[:red, false, :red],
markerstrokecolor=[:black, false, :black],
linecolor=:black)
