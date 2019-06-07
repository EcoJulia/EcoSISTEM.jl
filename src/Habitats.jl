using StatsBase
import Diversity.API._countsubcommunities
import Diversity.countsubcommunities
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using Compat
using RecipesBase
using EcoBase
import EcoBase: xmin, ymin, xcellsize, ycellsize, xcells, ycells, cellsize,
cells, xrange, yrange, xmax, ymax, indices, coordinates
import Plots: px

"""
    AbstractHabitat

Abstract supertype for all habitat types
"""
abstract type AbstractHabitat{H} <: EcoBase.AbstractGrid end

function countsubcommunities(ah::AbstractHabitat)
  return _countsubcommunities(ah)
end
xmin(ah::AbstractHabitat) = 0
ymin(ah::AbstractHabitat) = 0
xcellsize(ah::AbstractHabitat) = Float64(ah.size/km)
ycellsize(ah::AbstractHabitat) = Float64(ah.size/km)
xcells(ah::AbstractHabitat) = size(ah.matrix, 1)
ycells(ah::AbstractHabitat) = size(ah.matrix, 2)
indices(ah::AbstractHabitat) =
    hcat(collect.(convert_coords.(1:length(ah.matrix), xcells(ah)))...)'
indices(ah::AbstractHabitat, idx) = indices(ah)[:, idx]
coordinates(ah::AbstractHabitat) = indices(ah)
"""
    ContinuousHab <: AbstractHabitat{Float64}

This habitat subtype has a matrix of floats and a float grid square size
"""
mutable struct HabitatUpdate{D <: Unitful.Dimensions,
                             DT}
  changefun::Function
  rate::DT
end

function HabitatUpdate{D}(changefun, rate::DT) where {D <: Unitful.Dimensions, DT}
    typeof(dimension(rate * 1month)) == D || error("Failed to match types")
    return HabitatUpdate{D, DT}(changefun, rate)
end
GLOBAL_typedict["HabitatUpdate"] = HabitatUpdate

mutable struct ContinuousHab{C <: Number} <: AbstractHabitat{C}
  matrix::Array{C, 2}
  size::Unitful.Length
  change::HabitatUpdate
end

iscontinuous(hab::ContinuousHab{C}) where C = true
function eltype(hab::ContinuousHab{C}) where C
    return C
end
@recipe function f(H::ContinuousHab{C}) where C
    unitdict= Dict(K => "Temperature (K)", °C => "Temperature (°C)", mm => "Rainfall (mm)", kJ => "Solar Radiation (kJ)")
    h = ustrip.(H.matrix)
    seriestype  :=  :heatmap
    grid --> false
    right_margin --> 0.1px
    margin --> 10px
    aspect_ratio --> 1
    title --> unitdict[unit(C)]
    clims --> (minimum(h) * 0.99, maximum(h) * 1.01)
    xrange(H), yrange(H), h
end
@recipe function f(H::ContinuousHab{C}) where C <: Unitful.Temperature
    unitdict= Dict(K => "Temperature (K)", °C => "Temperature (°C)", mm => "Rainfall (mm)", kJ => "Solar Radiation (kJ)")
    h = ustrip.(uconvert.(°C, H.matrix))
    seriestype  :=  :heatmap
    grid --> false
    right_margin --> 0.1px
    margin --> 10px
    aspect_ratio --> 1
    title --> unitdict[unit(C)]
    clims --> (minimum(h) * 0.99, maximum(h) * 1.01)
    xrange(H), yrange(H), h
end
GLOBAL_typedict["ContinuousHab"] = ContinuousHab

mutable struct ContinuousTimeHab{C <: Number} <: AbstractHabitat{C}
  matrix::Array{C, 3}
  time::Int64
  size::Unitful.Length
  change::HabitatUpdate
end
@recipe function f(H::ContinuousTimeHab{C}, time::Int64) where C
    h = ustrip.(H.matrix)
    seriestype  :=  :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(C)]
    clims --> (minimum(h[:,:,time]) * 0.99, maximum(h[:,:,time]) * 1.01)
    xrange(H), yrange(H), h[:,:,time]
end

iscontinuous(hab::ContinuousTimeHab{C}) where C = true
function eltype(hab::ContinuousTimeHab{C}) where C
    return C
end
function _resettime!(hab::ContinuousTimeHab)
    hab.time = 1
end

function _countsubcommunities(hab::ContinuousHab)
  return length(hab.matrix)
end

function _countsubcommunities(hab::ContinuousTimeHab)
  return length(hab.matrix[:, :, 1])
end
GLOBAL_typedict["ContinuousTimeHab"] = ContinuousTimeHab

"""
    DiscreteHab <: AbstractHabitat{String}

This habitat subtype has a matrix of strings and a float grid square size
"""
mutable struct DiscreteHab{D} <: AbstractHabitat{D}
  matrix::Array{D, 2}
  size::Unitful.Length
  change::HabitatUpdate
end
@recipe function f(H::DiscreteHab{D}) where D
    h = ustrip.(H.matrix)
    seriestype  :=  :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(D)]
    clims --> (minimum(h) * 0.99, maximum(h) * 1.01)
    xrange(H), yrange(H), h
end


iscontinuous(hab::DiscreteHab) = false
function eltype(hab::DiscreteHab{D}) where D
    return D
end
function _countsubcommunities(hab::DiscreteHab)
  return length(hab.matrix)
end
GLOBAL_typedict["DiscreteHab"] = DiscreteHab

function _getdimension(hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab})
    return (size(hab.matrix, 1), size(hab.matrix, 2))
end
function _getsize(hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab})
  x = hab.size * size(hab.matrix, 1)
  y = hab.size * size(hab.matrix, 2)
  return x * y
end
import Base.size
function size(hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab}, d)
    return size(hab.matrix, d)
end

function _getgridsize(hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab})
  return hab.size
end

mutable struct HabitatCollection2{H1, H2} <: AbstractHabitat{Tuple{H1, H2}}
    h1::H1
    h2::H2
end
iscontinuous(hab::HabitatCollection2{H1, H2}) where {H1, H2} = [iscontinuous(hab.h1), iscontinuous(hab.h2)]
function eltype(hab::HabitatCollection2)
    return [eltype(hab.h1), eltype(hab.h2)]
end
@recipe function f(H::HabitatCollection2{H1, H2}) where {H1, H2}
    x, y = H.h1, H.h2
    layout := 2
    @series begin
        subplot := 1
        x
    end
    @series begin
        subplot := 2
        y
    end
end
GLOBAL_typedict["HabitatCollection2"] = HabitatCollection2
function _resettime!(hab::HabitatCollection2)
    _resettime!(hab.h1)
    _resettime!(hab.h2)
end


mutable struct HabitatCollection3{H1, H2, H3} <: AbstractHabitat{Tuple{H1, H2, H3}}
    h1::H1
    h2::H2
    h3::H3
end
iscontinuous(hab::HabitatCollection3) = [iscontinuous(hab.h1),
    iscontinuous(hab.h2), iscontinuous(hab.h3)]
function eltype(hab::HabitatCollection3)
    return [eltype(hab.h1), eltype(hab.h2), eltype(hab.h3)]
end

@recipe function f(H::HabitatCollection3{H1, H2, H3}) where {H1, H2, H3}
    x, y, z = H.h1, H.h2, H.h3
    layout := 3
    @series begin
        subplot := 1
        x
    end
    @series begin
        subplot := 2
        y
    end
    @series begin
        subplot := 3
        z
    end
end
GLOBAL_typedict["HabitatCollection3"] = HabitatCollection3
function _resettime!(hab::HabitatCollection3)
    _resettime!(hab.h1)
    _resettime!(hab.h2)
    _resettime!(hab.h3)
end

function _getdimension(hab::Union{HabitatCollection2, HabitatCollection3})
    return size(hab.h1.matrix, 1, 2)
end
function _getsize(hab::Union{HabitatCollection2, HabitatCollection3})
  return _getsize(hab.h1)
end
function size(hab::Union{HabitatCollection2, HabitatCollection3}, d)
    return size(hab.h1, d)
end

function _getgridsize(hab::Union{HabitatCollection2, HabitatCollection3})
  return _getgridsize(hab.h1)
end

function gethabitat(hab::AbstractHabitat, pos::Int64)
    x, y = convert_coords(pos, size(hab.matrix, 1))
    return hab.matrix[x, y]
end
function gethabitat(hab::AbstractHabitat, field::Symbol)
    return getfield(hab, field)
end
function gethabitat(hab::ContinuousTimeHab, pos::Int64)
    x, y = convert_coords(pos, size(hab.matrix, 1))
    return hab.matrix[x, y, hab.time]
end

function _countsubcommunities(hab::HabitatCollection2)
  return _countsubcommunities(hab.h1)
end
function _countsubcommunities(hab::HabitatCollection3)
  return _countsubcommunities(hab.h1)
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
        cluster = vcat(mapslices(x->M[x[1],x[2]].==1, neighbours, dims=2)...)
        # Find out if any of the neighbours have a value > 1, thus, have already
        # been assigned a cluster
        already=vcat(mapslices(x->M[x[1],x[2]].>1, neighbours, dims=2)...)
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
        already=vcat(mapslices(x->isassigned(T,x[1],x[2]), neighbours, dims=2)...)
        # If any already assigned then sample from most frequent neighbour traits
          if any(already)
            neighbours=neighbours[already,:]
            # Find all neighbour traits
            neighbour_traits=map(i->T[neighbours[i,1],neighbours[i,2]],
             1:size(neighbours,1))
             # Find which one is counted most often
            ind=argmax(map(x->sum(neighbour_traits.==x), types))
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
function randomniches(dimension::Tuple, types::Vector{Int64}, clumpiness::Float64,
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
    T = Array{Int64}(Compat.undef, dimension)
    # Fill in T with clusters already created
    map(x -> T[M.==x] .= sample(types, wv), 1:maximum(M))
    # Fill in undefined squares with most frequent neighbour
    _fill_in!(T, M, types, wv)
  end

  return DiscreteHab(T, gridsquaresize, HabitatUpdate{Unitful.Dimensions{()}}(NoChange, 0.0/s))
end

function simplehabitat(val::Unitful.Quantity, size::Unitful.Length,
  dim::Tuple{Int64, Int64})
  M = fill(val, dim)
  func = ChangeLookup[unit(val)]
  rate = 0.0 * unit(val)/s
  ContinuousHab(M, size, HabitatUpdate{typeof(dimension(val))}(func, rate))
end

function simplehabitat(val::Float64, size::Unitful.Length,
  dim::Tuple{Int64, Int64})
  M = fill(val, dim)

  ContinuousHab(M, size, HabitatUpdate{Unitful.Dimensions{()}}(NoChange, 0.0/s))
end

function tempgrad(minT::Unitful.Temperature{Float64}, maxT::Unitful.Temperature{Float64},
  size::Unitful.Length{Float64},
  dim::Tuple{Int64, Int64}, rate::Quantity{Float64, 𝚯*𝐓^-1})
  dim[1] > 1 ||
  error("First dimension should be greater than 1 for temperature gradient")
  M = Array{typeof(minT)}(Compat.undef, dim)
  total = dim[1]
  temp_range = collect(range(minT, stop = maxT, length = total))
  map(1:total) do seq
    M[seq, :] .= temp_range[seq]
  end
  ContinuousHab(M, size, HabitatUpdate{typeof(dimension(minT))}(TempChange, rate))
end

function raingrad(minR::Unitful.Length{Float64}, maxR::Unitful.Length{Float64}, size::Unitful.Length{Float64}, dim::Tuple{Int64, Int64}, rate::Quantity{Float64, 𝐋*𝐓^-1})
  dim[1] > 1 ||
  error("First dimension should be greater than 1 for temperature gradient")
  M = Array{typeof(minR)}(Compat.undef, dim)
  total = dim[1]
  rain_range = collect(range(minR, stop = maxR, length = total))
  map(1:total) do seq
    M[seq, :] .= rain_range[seq]
  end
  ContinuousHab(M, size, HabitatUpdate{typeof(dimension(minR))}(RainfallChange, rate))
end
