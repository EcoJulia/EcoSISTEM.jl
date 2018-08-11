using Phylo
using Distributions
using DataFrames
import Phylo.NodeNameIterator

"""
    jtree(SpeciesNames::Vector{String}, T::Type = String,
    dist::Distribution = Uniform(0, length(SpeciesNames)))

Function to create a random non-ultrametric phylogenetic tree, with branch
lengths drawn from a distribution, `dist`, and node information of type, `T`.
"""
function jtree(SpeciesNames::Vector{String}, T::Type = Int64,
  dist::Distribution = Uniform(0, length(SpeciesNames)))
  #Create tree
  tree1 = BinaryTree{LeafInfo, Vector{T}}(SpeciesNames)
  n = length(SpeciesNames)
  # If only one branch
  leaf_names = getleafnames(tree1)
  inner_names = addnodes!(tree1, n-1)
  if n == 1
    error("More species needed")
  elseif n == 2
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
#tree = jtree(5)
#plot(tree)

# Function to produce a matrix of all source and target nodes within a tree
function sou_tar(tree::AbstractTree, len::Bool)
  # Create an empty array
  branches = sort(getbranchnames(tree))
  # Run through each branch
  path =  map(branches) do brn
  getsource(tree, brn), gettarget(tree, brn)
end
  path = hcat(path, repmat([0.0], size(path, 1)))
    # Set length
    if len
      # Get length from tree
      path[:, 2] = map(branches) do brn
        getlength(tree, brn)
      end
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

# Function to find which row source and target node are in
function find_row(mat, sou, tar)
  find(map(mat[:,1]) do check
    (sou, tar) == check
  end)
end
# Function to find which row source and target nodes are in
function find_rows(mat, sour, targ)
  vcat(map(sour, targ) do sour, targ
    find_row(mat, sour, targ)
  end...)
end

function root_to_tips(tree)
  tips = collect(NodeNameIterator(tree, isleaf))
  paths = map(tips) do tps
    reverse(nodehistory(tree, tps))
  end
  paths
end
"""
    jcoal(SpeciesNames::Vector{String}, len::Float64, T::Type = String)

Function to create a random ultrametric phylogenetic tree, with branches
of length `len`, and node information of type, `T`.
"""
function jcoal(SpeciesNames::Vector{String}, len::Float64, T::Type=Int64)
  # Create random tree
  tree2=jtree(SpeciesNames, T)
  n = length(SpeciesNames)
  if n == 1
    tree2.branches[1].length = len
  else
  # Find tips of tree
  tips=NodeIterator(tree2, isleaf)
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
    rows=find_rows(paths_mat, pairs[:,1], pairs[:,2])
    # Calculate how long the path is already
    length_left=sum(paths_mat[rows, 2])
    # Remove rows that have already been assigned a length
    rows=rows[paths_mat[rows, 2] .== 0]
    # Calculate branch lengths from Dirichlet distribution
    pathlength=jdir(length(rows),1)*(len-length_left)
    # Assign branch lengths to matrix
    paths_mat[rows,2]=pathlength
    # Remove completed path
    deleteat!(paths, maxpath[2])
  end
  # Loop through and assign new branch lengths to tree
  for j in collect(1:size(paths_mat,1))
      tree2.branches[j].length=paths_mat[j, 2]
  end
end
  # Return ultrametric tree
  tree2
end



function arenoderecordsempty(tree::AbstractTree, nodes::Vector{String})
  map(nodes) do nod
  isempty(getnoderecord(tree, nod))
  end
end

"""
    assign_traits!(tree::BinaryTree, switch_rate::Vector{Float64},
    traits::Vector{Vector{String}})

Function to evolve categorical functional traits through a phylogenetic tree
with a specific switching rate.

"""
function assign_traits!(tree::BinaryTree, switch_rate::Vector{Float64},
          traits::DataFrame)

  # Check if tree already assigned
  check = arenoderecordsempty(tree, getnodenames(tree))
  all(check) || error("Some nodes already assigned traits")

  # Calculate all branch paths from root to tips
  tips = NodeNameIterator(tree, isleaf)
  root = first(NodeNameIterator(tree, isroot))

  paths = root_to_tips(tree)
  samp = DataFrame(colwise(rand,traits))
  names!(samp, names(traits))
  # Assign first node a trait randomly
  setnoderecord!(tree, root, samp)

  # Loop through all paths
  for i in paths
    # Split path into pairs of nodes
    pairs = pair(i)
    # Test if any branches have been assigned already
    unassigned = arenoderecordsempty(tree, i[i .!= root])
    pairs = pairs[unassigned, :]

    # Calculate how long the path is already
    len = distance(tree, first(pairs), last(pairs))
    # Calculate time to next switch
    # Draw switch times from exponential distribution
    # Stop when they are larger than the length of the path
    alltimes = map(switch_rate) do swt
      times = Array{Float64}(0)
        while(sum(times) < len)
          time_switch = jexp(swt*len)
          append!(times, time_switch)
        end
        times
      end

    # Sum up the event times cumulatively
    cum_times = cumsum.(alltimes)

    # Run through the branches for the path, assigning a trait
    for j in 1 : size(pairs, 1)
      sel_pair = pairs[j, :]
      # Calculate the length of the branch
      branch_len = distance(tree, first(sel_pair), last(sel_pair))

      # Calculate switch times
      num_switches =
        map(cum_times) do cum
          sum(cum .< branch_len)
        end

      # Find trait of last node
      last_label = getnoderecord(tree, first(sel_pair))

      map(num_switches) do number
        # If there are no switches, give same trait as previous node
        if number == 0
          set_node = last(sel_pair)
          setnoderecord!(tree, set_node, last_label)
        else
          # Else loop through for the required number of switches, sampling from list of traits
          while number > 0
            set_node = last(sel_pair)
            newtrait = map(names(traits)) do trt
                        col = traits[trt]
                        sample(col[col .!= last_label[trt]])
                       end
            newtrait = DataFrame(newtrait)
            names!(newtrait, names(traits))
            setnoderecord!(tree, set_node, newtrait)
            last_label = getnoderecord(tree, set_node)
            number = number - 1
          end
        end
      end
    end
  end
end

function assign_traits!(tree::BinaryTree, switch_rate::Float64,
  traits::Vector{Int64})
  if length(traits) == 1
      for n in getnodenames(tree)
          setnoderecord!(tree, n, traits)
      end
  else
      assign_traits!(tree, [switch_rate], [traits])
  end
end


function BM(T::Real, σ²::Float64, start::Float64, lab::String="")
  t = 0:T  # time
## first, simulate a set of random deviates
x = jnorm(0, sqrt(σ²),length(t) - 1)
## now compute their cumulative sum
x = cumsum(append!([start], x))
#Plots.plot(t, x, ylims=collect(extrema(x)).*[0.9,1.1], label=lab)
end

"""
    assign_traits!(tree::BinaryTree, start::Vector{Float64},
      σ²::Vector{Float64})

Function to evolve continuous functional traits through a phylogenetic tree
through Brownian motion, with a starting value, `start`, and rate, `σ²`.

"""
function assign_traits!(tree::BinaryTree, start::Vector{Float64},
  σ²::Vector{Float64})

  check = arenoderecordsempty(tree, getnodenames(tree))
  all(check) || warn("Some nodes already assigned traits")

  length(start) == length(σ²) || error("Start values and variance must have
  same number of traits")
  # Find all names of nodes
  names =  getnodenames(tree)
  # Sort by distance from root
  root = first(collect(NodeNameIterator(tree, isroot)))
  dist = map(names) do node
          distance(tree, root, node)
         end
  names = names[sortperm(dist)]
  # Loop through nodes in order of appearance
  for i in names
    if isroot(tree, i)
      setnoderecord!(tree, i, start)
    else
      pnt = getparent(tree, i)
      srt = getnoderecord(tree, pnt)
      path = first(get(branchroute(tree, pnt, i)))
      ln = getlength(getbranch(tree, path))

      newtrait = map(srt, σ²) do start, sig
                  last(BM(ln, sig, start))
                 end
      setnoderecord!(tree, i, newtrait)

    end
  end
end

function assign_traits!(tree::BinaryTree, start::Float64,
  σ²::Float64)
  assign_traits!(tree, [start], [σ²])
end

"""
    get_traits(tree::BinaryTree, tips::Bool=true)

Function to retrieve functional traits assigned to a phylogenetic tree, either
just tips or all nodes.

"""
function get_traits(tree::BinaryTree, tips::Bool=true)
   check = .!arenoderecordsempty(tree, getnodenames(tree))
   all(check) || error("All node records empty")
  if tips
    nodes = collect(NodeNameIterator(tree, isleaf))
  else
    nodes = getnodenames(tree)
  end
  df = vcat(getnoderecord.(tree, nodes)...)
  df[:species] = nodes
  return df
end
