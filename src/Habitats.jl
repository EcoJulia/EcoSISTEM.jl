using StatsBase
import Diversity.API._countsubcommunities
import Diversity.countsubcommunities
using Unitful
using Unitful.DefaultSymbols

"""
    AbstractHabitat

Abstract supertype for all habitat types
"""
abstract type AbstractHabitat{H} end

function countsubcommunities(ah::AbstractHabitat)
  return _countsubcommunities(ah)
end

"""
    ContinuousHab <: AbstractHabitat{Float64}

This habitat subtype has a matrix of floats and a float grid square size
"""
mutable struct HabitatUpdate{D <: Unitful.Dimension}
  changefun::Function
  rate::Quantity{Float64,Unitful.Dimensions{Tuple{D,
   Unitful.Dimension{:Time}(-1//1)}}}
end

mutable struct ContinuousHab{E <: Envtype} <: AbstractHabitat{Float64}
  matrix::Matrix{Unitful.Quantity{Float64}}
  size::Unitful.Length
  change::HabitatUpdate
end

function _countsubcommunities(hab::ContinuousHab)
  return length(hab.matrix)
end

"""
    DiscreteHab <: AbstractHabitat{String}

This habitat subtype has a matrix of strings and a float grid square size
"""
mutable struct DiscreteHab{E <: Envtype} <: AbstractHabitat{String}
  matrix::Matrix{String}
  size::Unitful.Length
  change::HabitatUpdate
end

function _countsubcommunities(niches::DiscreteHab)
  return length(niches.matrix)
end

# Function to create a habitat from a discrete set of types according to the
# Saura-Martinez-Millan algorithm (2000)
function _percolate!(M::AbstractMatrix, clumpiness::Real)
  for i in 1:(length(M))
    if junif(0, 1) < clumpiness
      M[i]=1
    end
  end
end

# Function to create clusters from percolated grid
function _identify_clusters!(M::AbstractMatrix)
  dimension=size(M)
  # Begin cluster count
  count=1
  # Loop through each grid square in M
  for x in 1:dimension[1]
    for y in 1:dimension[2]

      # If square is marked as 1, then apply cluster finding algorithm
      if M[x,y]==1.0
        # Find neighbours of M at this location
        neighbours=get_neighbours(M, y, x)
        # Find out if any of the neighbours also have a value of 1, thus, have
        # not been assigned a cluster yet
        cluster = vcat(mapslices(x->M[x[1],x[2]].==1, neighbours, 2)...)
        # Find out if any of the neighbours have a value > 1, thus, have already
        # been assigned a cluster
        already=vcat(mapslices(x->M[x[1],x[2]].>1, neighbours, 2)...)
        # If any already assigned neighbours, then assign the grid square to this
        # same type
          if any(already)
            neighbours=neighbours[already,:]
            M[x,y]=M[neighbours[1,1],neighbours[1,2]]
          # If none are assigned yet, then create a new cluster
          else
            count=count+1
            neighbours=neighbours[cluster,:]
            M[x,y]=count
            map(i->M[neighbours[i,1],neighbours[i,2]]=count, 1:size(neighbours,1))
        end
      end
    end
  end
end

function _fill_in!(T, M, types, wv)
  dimension = size(M)
  # Loop through grid of clusters
  for x in 1:dimension[1]
    for y in 1:dimension[2]
      # If square is zero then it is yet to be assigned
      if M[x,y]==0
        # Find neighbours of square on string grid
        neighbours=get_neighbours(T, y, x, 8)
        # Check if they have already been assigned
        already=vcat(mapslices(x->isassigned(T,x[1],x[2]), neighbours, 2)...)
        # If any already assigned then sample from most frequent neighbour traits
          if any(already)
            neighbours=neighbours[already,:]
            # Find all neighbour traits
            neighbour_traits=map(i->T[neighbours[i,1],neighbours[i,2]],
             1:size(neighbours,1))
             # Find which one is counted most often
            ind=indmax(map(x->sum(neighbour_traits.==x), types))
            # Assign this type to the grid square in T
            T[x,y]= types[ind]
          # If none are assigned in entire grid already,
          # sample randomly from traits
        elseif all(M.<=1)
            T[x,y]=sample(types, wv)
          # If some are assigned in grid, sample from these
          else
            T[x,y]=sample(T[M.>1])
        end
      end
    end
  end
end
"""
    randomniches(dimension::Tuple, types::Vector{String}, clumpiness::Float64, weights::Vector)

Function to create a `DiscreteHab` habitat of dimension `dimension`, made up of sampled
string types, `types`, that have a weighting, `weights` and clumpiness parameter,
`clumpiness`.
"""
function randomniches(dimension::Tuple, types::Vector{String}, clumpiness::Float64,
  weights::Vector, gridsquaresize::Unitful.Length)
  # Check that the proportion of coverage for each type matches the number
  # of types and that they add up to 1
  length(weights)==length(types) || error("There must be an area proportion for each type")
  sum(weights)==1 || error("Proportion of habitats must add up to 1")
  # Create weighting from proportion habitats
  wv = Weights(weights)

  # Create an empty grid of the right dimension
  M = zeros(dimension)

  # If the dimensions are too small for the algorithm, just use a weighted sample
  if dimension[1] <= 2 || dimension[2] <= 2
    T = sample(types, Weights(weights), dimension)
  else
    # Percolation step
    _percolate!(M, clumpiness)
    # Select clusters and assign types
    _identify_clusters!(M)
    # Create a string grid of the same dimensions
    T = Array{String}(dimension)
    # Fill in T with clusters already created
    map(x -> T[M.==x]=sample(types, wv), 1:maximum(M))
    # Fill in undefined squares with most frequent neighbour
    _fill_in!(T, M, types, wv)
  end

  return DiscreteHab{Niches}(T, gridsquaresize, HabitatUpdate{Unitful.Dimension{()}}(NoChange, 0.0/s))
end

function simplehabitat(val::Unitful.Quantity, size::Unitful.Length,
  dim::Tuple{Int64, Int64})
  M = Array{Unitful.Quantity}(dim)
  fill!(M, val)

  ContinuousHab{None}(M, size, HabitatUpdate{Unitful.Dimension{()}}(NoChange, 0.0/s))
end
function simplehabitat(val::Float64, size::Unitful.Length,
  dim::Tuple{Int64, Int64})
  M = Array{Float64}(dim)
  fill!(M, val)

  ContinuousHab{None}(M, size, HabitatUpdate{Unitful.Dimension{()}}(NoChange, 0.0/s))
end
function degradedhab(val::Union{Float64, Unitful.Quantity}, size::Unitful.Length,
  dim::Tuple{Int64, Int64}, rate::Quantity{Float64, typeof(𝐓^-1)})
  M = Array{Float64}(dim)
  fill!(M, val)

  ContinuousHab{None}(M, size, HabitatUpdate{Unitful.Dimension{()}}(HabitatLoss, rate))
end
function tempgrad(min::Unitful.Temperature{Float64}, max::Unitful.Temperature{Float64},
  size::Unitful.Length{Float64},
  dim::Tuple{Int64, Int64}, rate::Quantity{Float64, typeof(𝚯*𝐓^1)})
  M = Array{Unitful.Temperature}(dim)
  total = dim[1]
  temp_range = collect(linspace(min, max, total))
  map(1:total) do seq
    M[seq, :] = temp_range[seq]
  end
  ContinuousHab{Temp}(M, size, HabitatUpdate{Unitful.Dimension{:Temperature}}(TempChange, rate))
end
