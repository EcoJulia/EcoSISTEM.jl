using JLSO
using SparseArrays
using HCubature

"""
    AbstractLookup

Abstract type for lookups.
"""
abstract type AbstractLookup end

"""
    Lookup

Lookup houses information on `x`, `y` grid locations and the probability of
occurrence at the location for the species in question `p`. `pnew` and `moves`
are initially empty storage and written over by the movement step in update!().
`pnew` is the recalculated probability based on which directions are available
and `moves` is the number of moves to that grid location in that step.
"""
mutable struct Lookup
  x::Vector{Int64}
  y::Vector{Int64}
  p::Vector{Float64}
  pnew::Vector{Float64}
  moves::Vector{Int64}
end
Lookup(df::DataFrame) = Lookup(df[!, :X], df[!, :Y], df[!, :Prob],
zeros(Float64, nrow(df)),zeros(Int64, nrow(df)))

"""
    SpeciesLookup  <: AbstractLookup

SpeciesLookup holds `Lookup` information for each species in a vector, `species`.
"""
mutable struct SpeciesLookup  <: AbstractLookup
    species::Vector{Lookup}
end

function _getlookup(lookup::SpeciesLookup, sp::Int64)
    return lookup.species[sp]
end

function _symmetric_grid(grid::DataFrame)
   for x in 1:nrow(grid)
     if grid[x, 1] != grid[x, 2]
       push!(grid, hcat(grid[x, 2], grid[x, 1] , grid[x, 3]))
     end
   end
   for x in 1:nrow(grid)
     if (grid[x, 1] > 0)
       push!(grid, hcat(-grid[x, 1], grid[x, 2] , grid[x, 3]))
     end
     if (grid[x, 2] > 0)
       push!(grid, hcat(grid[x, 1], -grid[x, 2] , grid[x, 3]))
     end
     if (grid[x, 1] > 0 && grid[x, 2] > 0)
       push!(grid, hcat(-grid[x, 1], -grid[x, 2] , grid[x, 3]))
     end
   end
   grid
 end

 # Define gaussian kernel function
function _gaussian_disperse(r)
  exp(-((r[3]-r[1])^2+(r[4]-r[2])^2)) / π
end

function _2Dt_disperse(r, b)
    return((b - 1)/(π)) * (1 + ((r[3]-r[1])^2+(r[4]-r[2])^2))^-b
end

"""
    genlookups(hab::AbstractHabitat, mov::GaussianMovement)

Function to generate lookup tables, which hold information on the probability
of moving to neighbouring squares.
"""
function genlookups(hab::AbstractHabitat, mov::GaussianKernel)
  sd = (2 * mov.dist) / sqrt(pi)
  relsize =  _getgridsize(hab) ./ sd
  m = maximum(_getdimension(hab))
  p = mov.thresh
  return Lookup(_lookup(relsize, m, p, _gaussian_disperse))
end
function genlookups(hab::AbstractHabitat, mov::LongTailKernel)
    sd = (2 * mov.dist) / sqrt(pi)
    relsize =  _getgridsize(hab) ./ sd
    m = maximum(_getdimension(hab))
    p = mov.thresh
    b = mov.shape
    return EcoSISTEM.Lookup(EcoSISTEM._lookup(relsize, m, p, b, EcoSISTEM._2Dt_disperse))
end

function _lookup(relSquareSize::Float64, maxGridSize::Int64,
                pThresh::Float64, dispersalfn::F) where {F<:Function}
  # Create empty array
  lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])

  # Loop through directions until probability is below threshold
  k = 0
  m = 0
  count = 0
  while (k <= maxGridSize && m <= maxGridSize)
    count = count + 1
    calc_prob = hcubature(r -> dispersalfn(r),
      [0, 0, k*relSquareSize, m*relSquareSize],
      [relSquareSize, relSquareSize, (k+1)*relSquareSize, (m+1)*relSquareSize],
      maxevals= 10000)[1] / relSquareSize^2
    if m == 0 && calc_prob < pThresh
      break
    end
    if count == 1
      push!(lookup_tab, [k m calc_prob])
      k = k + 1
    elseif (calc_prob > pThresh && m <= k)
      push!(lookup_tab, [k m calc_prob])
      m = m + 1
    else
      m = 0
      k = k + 1
    end
  end
  # If no probabilities can be calculated, threshold is too high
  nrow(lookup_tab) != 0 || error("probability threshold too high")
  # Find all other directions
  lookup_tab = _symmetric_grid(lookup_tab)
  #info(sum(lookup_tab[:, 3]))
  # Normalise
  lookup_tab[!, :Prob] = lookup_tab[!, :Prob]/sum(lookup_tab[!, :Prob])
  lookup_tab
end


function _lookup(relSquareSize::Float64, maxGridSize::Int64,
                pThresh::Float64, b::Float64, dispersalfn::F
                ) where {F<:Function}
  # Create empty array
  lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])

  # Loop through directions until probability is below threshold
  k = 0
  m = 0
  count = 0
  while (k <= maxGridSize && m <= maxGridSize)
    count = count + 1
    calc_prob = hcubature(r -> dispersalfn(r, b),
      [0, 0, k*relSquareSize, m*relSquareSize],
      [relSquareSize, relSquareSize, (k+1)*relSquareSize, (m+1)*relSquareSize],
      maxevals=10000)[1] / relSquareSize^2
    if m == 0 && calc_prob < pThresh
      break
    end
    if count == 1
      push!(lookup_tab, [k m calc_prob])
       k = k + 1
    elseif (calc_prob > pThresh && m <= k)
      push!(lookup_tab, [k m calc_prob])
      m = m + 1
    else
      m = 0
      k = k + 1
    end
  end
  # If no probabilities can be calculated, threshold is too high
  nrow(lookup_tab) != 0 || error("probability threshold too high")
  # Find all other directions
  lookup_tab = _symmetric_grid(lookup_tab)
  #info(sum(lookup_tab[:, 3]))
  # Normalise
  lookup_tab[!, :Prob] = lookup_tab[!, :Prob]/sum(lookup_tab[!, :Prob])
  lookup_tab
end

@enum MovementType localMovement regionMovement

"""
    EpiLookup  <: AbstractLookup

EpiLookup holds sparse matrices of local (`locallookup`) and
commuting (`regionlookup`) moves for epi simulations.
"""
struct EpiLookup <: AbstractLookup
  locallookup::SparseMatrixCSC{Float64, Int32}
  regionlookup::SparseMatrixCSC{Float64, Int32}
  function EpiLookup(locallookup::SparseMatrixCSC{Float64, Int32}, regionlookup::SparseMatrixCSC{Float64, Int32})
      all(0 .<= locallookup.nzval .<= 1) || error("Home lookup values must be between 0 and 1")
      all(0 .<= regionlookup.nzval .<= 1) || error("Work lookup values must be between 0 and 1")
      return new(locallookup, regionlookup)
  end
end

function _getlookup(lookup::EpiLookup, id::Int64)
    return lookup.locallookup[id, :], lookup.regionlookup[id, :]
end

function genlookups(epienv::AbstractEpiEnv, mov::LongDistance, pop_size)
  total_size = (size(epienv.active, 1) * size(epienv.active, 2))
  # Column access so Js should be source grid cells
  Js = Int32.(mov.move_record[!, :from])
  # Is should be destination grid cells
  Is = Int32.(mov.move_record[!, :to])
  Vs = mov.move_record[!, :count]
  regionmoves = sparse(Is, Js, Vs, total_size, total_size)
  # Divide through by total population size
  regionmoves.nzval ./= pop_size[Is]
  regionmoves.nzval[isnan.(regionmoves.nzval)] .= 0
  # Make sure each row adds to one (probability of movement)
  summed = map(j -> sum(regionmoves[:, j]), unique(Js))
  summed[summed .== 0] .= 1.0
  regionmoves.nzval ./= summed
  return regionmoves
end

function genlookups(epienv::GridEpiEnv, mov::AlwaysMovement)
    total_size = (size(epienv.active, 1) * size(epienv.active, 2))
    # Generate grid ids and x,y coords for active cells only
    grid_locs = 1:total_size
    activity = epienv.active[1:end]
    grid_locs = grid_locs[activity]
    xys = convert_coords.(grid_locs, size(epienv.active, 1))

    # Collate all movement related parameters
    grid_size = _getgridsize(epienv.habitat)
    sd = [(2 .* k.dist) ./ sqrt(pi) for k in mov.kernels][activity]
    relsize =  grid_size ./ sd
    thresh = [k.thresh for k in mov.kernels][activity]
    grid_size /= unit(grid_size)

    # Calculate lookup probabilities for each grid location
    res = map((i, r, t) -> EcoSISTEM.genlookups(i, grid_locs, xys, grid_size, r, t, epienv), grid_locs, relsize, thresh)

    # Column vectors are source grid cells (repeated for each destination calculated)
    Js = vcat([fill(grid_locs[r], length(res[r][1])) for r in eachindex(res)]...)
    # Row vectors are destination grid cells
    Is = vcat([r[1] for r in res]...)
    Vs = vcat([r[2] for r in res]...)
    return sparse(Int32.(Is), Int32.(Js), Vs, total_size, total_size)
end

function genlookups(from::Int64, to::Vector{Int64}, xys::Array{Tuple{Int64,Int64},1}, grid_size::Float64, relsize::Float64, thresh::Float64, epienv::GridEpiEnv)
    x, y = xys[to .== from][1]
    maxX = ceil(Int64, x + 1/relsize); minX = ceil(Int64, x - 1/relsize)
    maxY = floor(Int64, y + 1/relsize); minY = floor(Int64, y - 1/relsize)
    keep = [(i[1] <= maxX) & (i[2] <= maxY) & (i[1] >= minX) & (i[2] >= minY) for i in xys]
    to = to[keep]
    probs = [_lookup((x = x, y = y), (x = i[1], y = i[2]), relsize, _gaussian_disperse) for i in xys[keep]]
    keep = probs .> thresh
    probs = probs[keep]
    probs ./= sum(probs)
    return to[keep], probs
end

function _lookup(from::NamedTuple, to::NamedTuple, relSquareSize::Float64, dispersalfn::F) where {F<:Function}
    return calc_prob = hcubature(r -> dispersalfn(r),
      [from.y *relSquareSize - relSquareSize, from.x * relSquareSize - relSquareSize, to.y * relSquareSize - relSquareSize, to.x * relSquareSize - relSquareSize],
      [from.y * relSquareSize, from.x * relSquareSize, to.y * relSquareSize, to.x * relSquareSize],
      maxevals= 100, rtol = 0.01)[1] / relSquareSize^2
end
