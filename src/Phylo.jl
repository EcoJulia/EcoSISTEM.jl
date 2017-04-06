using PhyloTrees
#using Simulation
#Function to create random tree with set number of tips
function jtree(n::Int64, T::Type = String, dist::Distribution = Uniform(0,n))
  #Create tree and set initial two branches
  tree1 = Tree{T, Vector}(); addnodes!(tree1, 2*n-1)
  # If only one branch
  if n==1
    addnodes!(tree1, 1)
    addbranch!(tree1, 1, 2, rand(dist))
  else
  map(target->addbranch!(tree1, 1, target, rand(dist)),2:3)
  if n>2
  nodes=3
  #Run through randomly choosing the next node to split
  i=collect(2:3)
  # Run until reaches max number of nodes
    while nodes<(2*n-1)
      #Choose random node
      n1=rand(i)
      # Split into two
      nodes=nodes+2
      # Make new branches random length
      map(var-> addbranch!(tree1,n1, var, rand(dist)),(nodes-1):(nodes))
      # Remove nodes already split
      deleteat!(i, findin(i,n1))
      # Add newly created tips into the mix
      append!(i,collect((nodes-1):nodes))
    end
  end
end
# Return tree
tree1
end

# Function to produce a matrix of all source and target nodes within a tree
function sou_tar(tree, len::Bool)
  # Create an empty array
  path=Array{Real}(length(tree.branches),3)
  # Run through each branch
  for x in 1:length(tree.branches)
    # Get source and target nodes of branch
    path[x,[1,2]]=[tree.branches[x].source,tree.branches[x].target]
    # Set length
    if len
      # Get length from tree
      path[x,3]=tree.branches[x].length.value
    else
      # Create blank
    path[x,3]=0
  end
  end
  path
end
# Function to split vector into pairs of numbers
function pair(vec)
  # Calc number of pairs
  npairs=length(vec)-1
  # Create empty array
  newvec=Array{Int64}(npairs,2)
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
function findtips(tree)
  findin(map(isleaf,repmat([tree],length(tree.nodes)),collect(1:length(tree.nodes))), true)
end
import PhyloTrees.findroot
function find_root(tree)
  findin(map(isroot,repmat([tree],length(tree.nodes)),collect(1:length(tree.nodes))), true)
end
function root_to_tips(tree)
  tips=findtips(tree)
  root=find_root(tree)
  map(nodepath,repmat([tree],endof(tips)), repmat(root,endof(tips)), tips)
end
# Function to create a random ultrametric tree
function jcoal(n::Int64, len::Real, T::Type=String)
  # Create random tree
  tree2=jtree(n, T)

  if n == 1
    tree2.branches[1].length = len
  else
  # Find tips of tree
  tips=findtips(tree2)
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


function assign_traits!(tree::Tree, switch_rate::Real, traits::Vector)
  # Check if tree already assigned
  check = map(a->haslabel(tree.nodes[a]), 1:length(tree.nodes))
  !all(check) || error("Some nodes already assigned traits")
  # Calculate all branch paths from root to tips
  tips=findleaves(tree)
  paths = map(i-> reverse(nodepath(tree, i)), tips)
  # Assign first node a trait randomly
  setlabel!(tree.nodes[1], sample(traits))
  # Loop through all paths
  for i in (1 : length(paths))
    # Choose first path
    choosepath = paths[i]
    # Split path into pairs of nodes
    pairs = pair(choosepath)
    assigned=false
    # Test if any branches have been assigned already
    if size(pairs, 1) > 1
      test_assigned = map(a -> haslabel(tree.nodes[a]), pairs)
      assigned = mapslices(sum, test_assigned, 2) .== 2
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
      tree.branches[branch_path[1]].data = cum_times[cum_times .< branch_len] + prev_branch_len
      # Calculate switch times
      num_switches = sum(cum_times .< branch_len)

      # Find trait of last node
      labels = map(a -> haslabel(tree.nodes[a]), sel_pair)
      last_node = maximum(sel_pair[labels])
      last_label = getlabel(tree.nodes[last_node])

      # If there are no switches, give same trait as previous node
      if num_switches == 0
        set_node = minimum(sel_pair[!labels])
        setlabel!(tree.nodes[set_node], last_label)
      else
      # Else loop through for the required number of switches, sampling from list of traits
        while num_switches > 0
          set_node = minimum(sel_pair[!labels])
          setlabel!(tree.nodes[set_node], sample(traits[traits .!= last_label]))
          last_label = getlabel(tree.nodes[set_node])
          num_switches = num_switches - 1
        end
      end
    end
  end
# Return tree
#tree
end


function get_traits(tree::Tree, tips::Bool=false)
  if tips
    tipnodes= findleaves(tree)
    map(a->getlabel(tree.nodes[a]), 1:length(tree.nodes))[tipnodes]
  else
    map(a->getlabel(tree.nodes[a]), 1:length(tree.nodes))
  end
end

function get_times(tree::Tree)
  map(a->get(tree.branches[a]), 1:length(tree.branches))
end



function BM(T::Real,σ²::Float64, start::Float64, lab::String="")
  t = 0:T  # time
## first, simulate a set of random deviates
x = jnorm(0, sqrt(σ²),length(t) - 1)
## now compute their cumulative sum
x = cumsum(append!([start], x))
#Plots.plot(t, x, ylims=collect(extrema(x)).*[0.9,1.1], label=lab)
end

tree=jcoal(5, 10, Real)

function assign_traits!(tree::Tree, start::Float64, σ²:: Float64)
  check = map(a->!isnull(tree.nodes[a].data), 1:length(tree.nodes))
  !all(check) || error("Some nodes already assigned traits")
  # Start out with branches 1 & 2
  len = map(x-> tree.branches[x].length, 1:2)
  traits = map(x->BM(x, σ², start)[end], len)
  #loop through nodes
  for i in eachindex(tree.nodes)
    if isroot(tree.nodes[i])
      tree.nodes[i].data = start
    else
      pnt = parentnode(tree,i)
      srt = get(tree.nodes[pnt].data)
      path = branchpath(tree, pnt, i)[1]
      ln = tree.branches[path].length
      tree.nodes[i].data = BM(ln, σ², srt)[end]
    end
  end
  tree
end
function get_data(tree::Tree, tips::Bool=false)
  if tips
    tipnodes= findleaves(tree)
    map(a->get(tree.nodes[a].data), 1:length(tree.nodes))[tipnodes]
  else
    map(a->get(tree.nodes[a].data), 1:length(tree.nodes))
  end
end
assign_traits!(tree, 2.0, 0.5)
get_data(tree)
