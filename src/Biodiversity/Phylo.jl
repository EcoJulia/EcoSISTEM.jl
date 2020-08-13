using Phylo
using Distributions
using DataFrames
using Phylo: NodeNameIterator
using Compat

function reroot!(tree::AbstractTree, node::String)
    root = collect(nodenamefilter(isroot, tree))[1]
    deletenode!(tree, node)
    addnode!(tree, node)
    addnode!(tree, "NewRoot")
    addbranch!(tree, "NewRoot", root)
    addbranch!(tree, "NewRoot", node)
end

function resettraits!(tree::AbstractTree)
    nodes = getnodenames(tree)
    for i in nodes
        setnodedata!(tree, i, DataFrame())
    end
end

# Function to split vector into pairs of numbers
function pair(vec)
  # Calc number of pairs
  npairs=length(vec)-1
  # Create empty array
  newvec=Array{String}(Compat.undef, npairs, 2)
  # Split into pairs
  for i in collect(1:npairs)
    newvec[i,:]=vec[i:(i+1)]
  end
  newvec
end

function root_to_tips(tree)
  tips = collect(nodenamefilter(isleaf, tree))
  paths = map(tips) do tps
    reverse(nodehistory(tree, tps))
  end
  paths
end

function arenoderecordsempty(tree::AbstractTree, nodes::Vector{String})
  map(nodes) do nod
  isempty(getnodedata(tree, nod))
  end
end

"""
    assign_traits!(tree::AbstractTree, switch_rate::Vector{Float64},
    traits::Vector{Vector{String}})

Function to evolve categorical functional traits through a phylogenetic tree
with a specific switching rate.

"""
function assign_traits!(tree::AbstractTree, switch_rate::Vector{Float64},
          traits::DataFrame)
  # Check if tree already assigned
  check = arenoderecordsempty(tree, collect(getnodenames(tree)))
  all(check) || error("Some nodes already assigned traits")

  # Calculate all branch paths from root to tips
  tips = collect(nodenamefilter(isleaf, tree))
  root = first(collect(nodenamefilter(isroot, tree)))

  paths = root_to_tips(tree)
  samp = DataFrame(hcat([rand(col) for col = eachcol(traits)]), names(traits))
  # Assign first node a trait randomly
  setnodedata!(tree, root, samp)

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
      times = Array{Float64}(Compat.undef, 0)
        while(sum(times) < len)
          time_switch = rand(Exponential(swt*len))
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
      last_label = getnodedata(tree, first(sel_pair))

      map(num_switches) do number
        # If there are no switches, give same trait as previous node
        if number == 0
          set_node = last(sel_pair)
          setnodedata!(tree, set_node, last_label)
        else
          # Else loop through for the required number of switches, sampling from list of traits
          while number > 0
            set_node = last(sel_pair)
            newtrait = map(names(traits)) do trt
                        col = traits[!, trt]
                        sample(col[col .!= last_label[:, trt]])
                       end
            newtrait = DataFrame(hcat(newtrait), names(traits))
            setnodedata!(tree, set_node, newtrait)
            last_label = getnodedata(tree, set_node)
            number = number - 1
          end
        end
      end
    end
  end
end

function assign_traits!(tree::AbstractTree, switch_rate::Float64,
  traits::DataFrame)
  return assign_traits!(tree, [switch_rate], traits)
end

"""
    BM(T::Real, σ²::Float64, start::Float64, lab::String="")

Function to evolve a Real value through Brownian motion, with a starting value,
 `start`, and rate, `σ²`.

"""
function BM(T::Real, σ²::Float64, start::Float64, lab::String="")
  t = 0:T  # time
    # first, simulate a set of random deviates
    x = rand(Normal(0, sqrt(σ²)),length(t))
    # now compute their cumulative sum
    x = cumsum(append!([start], x))
    return x
end

"""
    assign_traits!(tree::AbstractTree, start::Vector{Float64},
      σ²::Vector{Float64})

Function to evolve continuous functional traits through a phylogenetic tree
through Brownian motion, with a starting value, `start`, and rate, `σ²`.

"""
function assign_traits!(tree::AbstractTree, traits::DataFrame)

  # Warning for nodes that have already been assigned traits
  check = arenoderecordsempty(tree, collect(getnodenames(tree)))
  all(check) || Compat.@warn "Some nodes already assigned traits"

  # Check that the length of the starting values and variances are the same
  length(traits[:start]) == length(traits[:σ²]) || error("Start values and variance must have
  same number of traits")

  # Find all names of nodes
  names =  collect(getnodenames(tree))

  # Sort by distance from root
  root = first(collect(NodeNameIterator(tree, isroot)))
  dist = map(names) do node
          distance(tree, root, node)
         end
  names = names[sortperm(dist)]

  # Loop through nodes in order of appearance
  for i in names
      # Set start value of BM to root
    if isroot(tree, i)
      setnodedata!(tree, i, traits)
    else
      # Get value information for parent node
      pnt = getparent(tree, i)
      srt = getnodedata(tree, pnt)[:start]
      # Find length of path between parent and child node
      path = first(branchroute(tree, pnt, i))
      ln = getlength(tree, getbranch(tree, path))

      # Run BM model on each trait and set record
      newtrait = map(srt, traits[:σ²]) do start, sig
                  last(BM(ln, sig, start))
                 end
      newdat = DataFrame(start = newtrait, σ² = traits[:σ²])
      setnodedata!(tree, i, newdat)

    end
  end
end

"""
    get_traits(tree::AbstractTree, tips::Bool=true)

Function to retrieve functional traits assigned to a phylogenetic tree, either
just tips or all nodes.

"""
function get_traits(tree::AbstractTree, tips::Bool=true)
   check = .!arenoderecordsempty(tree, collect(getnodenames(tree)))
   all(check) || error("All node records empty")
  if tips
    nodes = collect(nodenamefilter(isleaf, tree))
  else
    nodes = getnodenames(tree)
  end
  df = vcat(map(n->getnodedata(tree, n), nodes)...)
  df[!, :species] = nodes
  return df
end
