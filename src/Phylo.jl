using PhyloTrees
using Distributions
using Plots
#using Simulation
#Function to create random tree with set number of tips
function jtree(n::Int64, T::Type = String, dist::Distribution = Uniform(0,n))
  #Create tree
  tree1 = NodeTree(n, nodetype = Vector{T})
  # If only one branch
  leaf_names = getleafnames(tree1)
  inner_names = addnodes!(tree1, n-1)
  if n == 1
    error("More species needed")
  elseif n==2
    map(x -> addbranch!(tree1, inner_names[1], x, rand(dist)), leaf_names)
  else
    inner_names_copy = copy(inner_names)
    choose_node = inner_names[1]
    avail_nodes = inner_names
    # Set up internal structure
    # Connect all internal nodes together
    while !isempty(avail_nodes)
      deleteat!(inner_names, find(inner_names .== choose_node))

      if length(avail_nodes) > 1
        howmany = rand(1:2)
      else
        howmany = 1
      end

      targets = sample(avail_nodes, howmany, replace = false)
      map(x -> addbranch!(tree1, choose_node, x, rand(dist)), targets)

      avail_nodes = inner_names[map(x -> indegree(tree1, x) < 1, inner_names)]
      choose_node = rand(targets, 1)[1]
    end

    # Add on leaf nodes
    count = 0
    while !isempty(leaf_names) && count <= length(inner_names_copy)
      count = count + 1
      if outdegree(tree1, inner_names_copy[count]) < 2
        targets = sample(leaf_names, 2 - outdegree(tree1, inner_names_copy[count]),
                            replace = false)
        map(x -> addbranch!(tree1, inner_names_copy[count], x, rand(dist)), targets)
        deleteat!(leaf_names, findin(leaf_names, targets))
      end
    end
  end
  tree1
end
tree = jtree(5)
plot(tree)

# Function to produce a matrix of all source and target nodes within a tree
function sou_tar(tree::AbstractTree, len::Bool)
  # Create an empty array
  path=Array{String}(length(tree.branches),2)
  # Run through each branch
  map(x -> path[x, :] = [tree.branches[x].source, tree.branches[x].target], 1:size(path, 1))
  path = hcat(path, repmat([0.0], size(path, 1)))
    # Set length
    if len
      # Get length from tree
      map(x -> path[x, 3] = tree.branches[x].length, 1:size(path, 1))
  end
  path
end
# Function to split vector into pairs of numbers
function pair(vec)
  # Calc number of pairs
  npairs=length(vec)-1
  # Create empty array
  newvec=Array{String}(npairs,2)
  # Split into pairs
  for i in collect(1:npairs)
    newvec[i,:]=vec[i:(i+1)]
  end
  newvec
end

# Function to find if a node is a leaf
#function isleaf(tree, node)
#  length(tree.nodes[node].out)==0
#end
# Function to find if a node is a root
#function isroot(tree, node)
#  length(tree.nodes[node].in)==0
#end

# Function to find which row source and target node are in
function find_row(mat, sou, tar)
  row=map(==, mat[:,1:2], repmat([sou tar],size(mat,1)))
  findin(row[:,1]&row[:,2], true)
end
# Function to find which row source and target nodes are in
function find_rows(mat, sour, targ)
  hcat(map(find_row, repmat([mat], length(sour)), sour, targ)...)
end
# Function to find tip nodes of tree
#function findtips(tree)
#  findin(map(isleaf,repmat([tree],length(tree.nodes)),collect(1:length(tree.nodes))), true)
#end
#import PhyloTrees.findroot
#function find_root(tree)
#  findin(map(isroot,repmat([tree],length(tree.nodes)),collect(1:length(tree.nodes))), true)
#end
function root_to_tips(tree)
  tips=findleaves(tree)
  root=findroots(tree)[1]
  map(nodepath,repmat([tree],endof(tips)), repmat([root],endof(tips)), tips)
end
# Function to create a random ultrametric tree
function jcoal(n::Int64, len::Real, T::Type=String)
  # Create random tree
  tree2=jtree(n, T)

  if n == 1
    tree2.branches[1].length = len
  else
  # Find tips of tree
  tips=findleaves(tree2)
  # Find paths of each tip from root
  paths=root_to_tips(tree2)
  # Create array of source and target nodes
  paths_mat=sou_tar(tree2, false)

  # Loop through all the paths and calculate branch lengths
  while length(paths)>0
    # How many branches does each path have?
    howmany=map(length, paths)
    # Select maximum
    maxpath=findmax(howmany)
    choosepath=paths[maxpath[2]]
    # Split path into pairs of nodes
    pairs=pair(choosepath)
    # Find the row each pair of nodes corresponds to
    rows=find_rows(paths_mat,pairs[:,1],pairs[:,2])
    # Calculate how long the path is already
    length_left=sum(paths_mat[rows,3])
    # Remove rows that have already been assigned a length
    rows=rows[map(==,vcat(paths_mat[rows, 3]...),repmat([0],length(rows)))]
    # Calculate branch lengths from Dirichlet distribution
    pathlength=jdir(length(rows),1)*(len-length_left)
    # Assign branch lengths to matrix
    paths_mat[rows,3]=pathlength
    # Remove completed path
    deleteat!(paths,maxpath[2])
  end
  # Loop through and assign new branch lengths to tree
  for j in collect(1:size(paths_mat,1))
      tree2.branches[j].length=paths_mat[j,3]
  end
end
  # Return ultrametric tree
  tree2
end
function hasnoderecord(tree::AbstractTree, node::String)
  return !isempty(getnoderecord(tree, node))
end
function hasnoderecords(tree::AbstractTree, nodes::Vector{String})
  map(x -> hasnoderecord(tree, x), nodes)
end
function findall(tree::AbstractTree)
  return vcat(keys(tree.nodes)...)
end


function assign_traits!(tree::NodeTree, switch_rate::Real, traits::Vector)
  # Check if tree already assigned
  check = map(a->hasnoderecord(tree, a), keys(tree.nodes))
  !all(check) || error("Some nodes already assigned traits")
  # Calculate all branch paths from root to tips
  tips=findleaves(tree)
  paths = map(i-> reverse(nodepath(tree, i)), tips)
  # Assign first node a trait randomly
  setnoderecord!(tree, "Node 1", [sample(traits)])
  # Loop through all paths
  for i in (1 : length(paths))
    # Choose first path
    choosepath = paths[i]
    # Split path into pairs of nodes
    pairs = pair(choosepath)
    assigned=false
    # Test if any branches have been assigned already
    if size(pairs, 1) > 1
      test_assigned = mapslices(a -> hasnoderecords(tree, a), pairs, 1)
      assigned = mapslices(all, test_assigned, 2)
      assigned = vcat(assigned...)
      # If they have been assigned, remove from node pairs
      pairs = pairs[!assigned, :]
    end

  # Calculate how long the path is already
  path = mapslices(a -> branchpath(tree, a[1], a[2]), pairs, 2)
  len = sum(map(a -> tree.branches[a].length, path))


  # Calculate time to next switch
  sumtimes = Array{Vector{Float64}}(length(paths))
  times = Array{Float64}(0)
  # Draw switch times from exponential distribution
  # Stop when they are larger than the length of the path
  while(sum(times) < len)
  time_switch = jexp(switch_rate*len)
  append!(times, time_switch)
  end
  # Sum up the event times cumulatively
  cum_times = cumsum(times)

  # Run through the branches for the path, assigning a trait
    for j in 1 : size(pairs, 1)

      sel_pair = pairs[j, :]
      #Get the branch the node pairs are found on
      branch_path = branchpath(tree, sel_pair[1], sel_pair[2])
      # Calculate the length of the branch
      branch_len = tree.branches[branch_path[1]].length
      # If previous branches have been assigned then need to add previous branch lengths to switch times
      prev_branch_len = 0.0
      if any(assigned)
        prev_path=branchpath(tree, sel_pair[1])
        prev_branch_len=sum(map(i-> tree.branches[i].length, prev_path))
      end
      # Assign switch times as branch data
      #tree.branches[branch_path[1]].data = cum_times[cum_times .< branch_len] + prev_branch_len
      # Calculate switch times
      num_switches = sum(cum_times .< branch_len)

      # Find trait of last node
      labels = map(a -> hasnoderecord(tree, a), sel_pair)
      last_node = maximum(sel_pair[labels])
      last_label = getnoderecord(tree, last_node)

      # If there are no switches, give same trait as previous node
      if num_switches == 0
        set_node = minimum(sel_pair[!labels])
        setnoderecord!(tree, set_node, last_label)
      else
      # Else loop through for the required number of switches, sampling from list of traits
        while num_switches > 0
          set_node = minimum(sel_pair[!labels])
          setnoderecord!(tree, set_node, [sample(traits[traits .!= last_label])])
          last_label = getnoderecord(tree, set_node)
          num_switches = num_switches - 1
        end
      end
    end
  end
end


function get_traits(tree::NodeTree, tips::Bool=true)
  if tips
    tipnodes= findleaves(tree)
    return hcat(map(a->getnoderecord(tree, a)[1], tipnodes), tipnodes)
  else
    allnodes= findall(tree)
    return hcat(map(a->getnoderecord(tree, a)[1], allnodes), allnodes)
  end
end

#function get_times(tree::Tree)
#  map(a->get(tree.branches[a]), 1:length(tree.branches))
#end



function BM(T::Real,σ²::Float64, start::Float64, lab::String="")
  t = 0:T  # time
## first, simulate a set of random deviates
x = jnorm(0, sqrt(σ²),length(t) - 1)
## now compute their cumulative sum
x = cumsum(append!([start], x))
#Plots.plot(t, x, ylims=collect(extrema(x)).*[0.9,1.1], label=lab)
end



tree = jcoal(5, 10, Float64)


function assign_traits!(tree::NodeTree, start::Float64, σ²:: Float64)
  check = map(a->hasnoderecord(tree, a), findall(tree))
  !all(check) || error("Some nodes already assigned traits")
  # Start out with branches 1 & 2
  len = map(x-> tree.branches[x].length, 1:2)
  traits = map(x->BM(x, σ², start)[end], len)
  # Find all names of nodes
  names =  findall(tree)
  # Sort by distance from root
  dist = map(x -> distance(tree, "Node 1", x), names)
  names = names[sortperm(dist)]
  # Loop through nodes in order of appearance
  for i in names
    if isroot(tree, i)
      setnoderecord!(tree, i, [start])
    else
      pnt = parentnode(tree, i)
      srt = getnoderecord(tree, pnt)[1]
      path = branchpath(tree, pnt, i)[1]
      ln = tree.branches[path].length
      setnoderecord!(tree, i, [BM(ln, σ², srt[1])[end]])

    end
  end
  tree
end

assign_traits!(tree, 2.0, 0.1)
