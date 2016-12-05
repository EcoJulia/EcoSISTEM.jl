
%Load packages
using PhyloTrees
using Plots
using Distributions
#Pkg.update()


#Function to create random tree with set number of tips
function jtree(n::Int64, dist::Distribution = Uniform(0,n))
  #Function to add branch of random length
  function f(sou,var, dist)
    addbranch!(tree1, sou, var, rand(dist))
  end
  #Create tree and set initial two branches
  tree1 = Tree{String, Void}(); addnodes!(tree1, 2*n-1)
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
      map(f,[n1,n1],(nodes-1):(nodes), [dist,dist])
      # Remove nodes already split
      deleteat!(i, findin(i,n1))
      # Add newly created tips into the mix
      append!(i,collect((nodes-1):nodes))
    end
  end
# Return tree
tree1
end
randtree=jtree(17, Exponential(0.1))
Plots.plot(randtree)

# Function to sample randomly from the Uniform distribution
function junif(a,b)
  rand(Uniform(a,b))
end

# Function to sample randomly from the Dirichlet distribution
function jdir(k,a)
  rand(Dirichlet(k,a))
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
  newvec=Array{Any}(npairs,2)
  # Split into pairs
  for i in collect(1:npairs)
    newvec[i,:]=vec[i:(i+1)]
  end
  newvec
end

# Function to find if a node is a leaf
function isleaf(tree, node)
  length(tree.nodes[node].out)==0
end
# Function to find if a node is a root
function isroot(tree, node)
  length(tree.nodes[node].in)==0
end

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
function findroot(tree)
  findin(map(isroot,repmat([tree],length(tree.nodes)),collect(1:length(tree.nodes))), true)
end
function root_to_tips(tree)
  tips=findtips(tree)
  root=findroot(tree)
  map(nodepath,repmat([tree],endof(tips)), repmat(root,endof(tips)), tips)
end
# Function to create a random ultrametric tree
function jcoal(n::Int64, len::Real)
  # Create random tree
  tree2=jtree(n)
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
  # Return ultrametric tree
  tree2
end
coaltree=jcoal(14, 5)
Plots.plot(coaltree)

function jexp(theta, n::Int64=1)
  rand(Exponential(theta), n)
end
function jpois(gamma, n::Int64=1)
  rand(Poisson(gamma), n)
end

function node_find(tree)
  # Create an empty array
  path=Array{Any}(length(tree.nodes),3)
  # Run through each branch
  for x in 1:length(tree.nodes)
    # Get in and out nodes of branch
    path[x,1]=x
    path[x,2]=map(a->a.source, tree.branches[tree.nodes[x].in])
    path[x,3]= map(a->a.target, tree.branches[tree.nodes[x].out])
  end
  path
end
function jbinom(n::Int64, p::Real)
  rand(Binomial(n,p), n)
end
tree=jcoal(3, 10)
Plots.plot(tree)
switch_rate=0.5
function assign_trait(tree, switch_rate::Real, traits)
  # Calculate all branch paths
paths=root_to_tips(tree)
paths_mat=sou_tar(tree, true)
# Calculate all nodes
nodes_mat=node_find(tree)
nodes_mat=hcat(nodes_mat, zeros(size(nodes_mat,1)))
# Assign first node a trait randomly
nodes_mat[1,4]=sample(traits)
for i in (1:length(paths))
choosepath=paths[i]
# Split path into pairs of nodes
pairs=pair(choosepath)
# Find the row each pair of nodes corresponds to
rows=find_rows(paths_mat,pairs[:,1],pairs[:,2])
# Calculate how long the path is already
len=sum(paths_mat[rows,3])
# Check if any of the nodes have already been assigned?

# Calculate time to next switch
sumtimes=Array{Vector{Float64}}(length(paths))
times=Array{Float64}(0)
while(sum(times)<len)
time_switch= jexp(switch_rate*len)
append!(times,time_switch)
end
cum_times=cumsum(times)
for j in eachindex(pairs)
num_switches= sum(cum_times.<len)
# Find trait of last node

sample(traits, num_switches)

end


# Return trait matrix
nodes_mat

end

tree=jcoal(3,10)
traits=assign_trait(tree, 0.9, [1,2])
Plots.plot(tree,markershape=:circle, markercolor=hcat(traits[:,4]...))

function discrete_trait(tree, switch_rate::Real, traits)
  # Calculate all branch paths
paths=root_to_tips(tree)
paths_mat=sou_tar(tree, true)
# Calculate all nodes
nodes_mat=node_find(tree)
nodes_mat=hcat(nodes_mat, zeros(size(nodes_mat,1)))
# Assign first node a trait randomly
nodes_mat[1,4]=sample(traits)
# Then loop through the rest of the paths until all filled in
for i in collect(1:length(paths))
  # Calculate rate based on branch length
  rate= switch_rate*paths_mat[i,3]
  # Calculate time to next switch
  time_switch= jexp(rate)
  # Find parent node trait
  ind=indexin([nodes_mat[i-1,4]], traits)
# If time to switch larger than distance between nodes, switch to another trait

  if time_switch[1]>paths_mat[i-1,3]
    nodes_mat[i,4]=sample(deleteat!(traits,ind))
    append!(traits,[nodes_mat[nodes_mat[i,2],4][1]])
  else
    # else stay the same as parent node
    nodes_mat[i,4]=nodes_mat[nodes_mat[i,2],4][1]
  end

end

# Return trait matrix
nodes_mat

end

function jnorm(μ,σ,n::Int=1)
  rand(Normal(μ,σ),n)
end

function BM(T::Real,σ²::Float64, start::Float64, lab::String="")
  t = 0:T  # time
## first, simulate a set of random deviates
x = jnorm(0, sqrt(σ²),length(t) - 1)
## now compute their cumulative sum
x = cumsum(append!([start], x))
Plots.plot(t, x, ylims=collect(extrema(x)).*[0.9,1.1], label=lab)
end
BM(100.0,0.01,10.0)

function continuous_trait(tree::Tree, σ²::Float64, start::Float64)


Plots.plot(jcoal(10, 10000))
cont_trait(tree)

landscape(n.traits::Int64, dim::Real,size::Real)
