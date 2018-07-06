using Diversity
using Cubature
using DataFrames
using Unitful
using myunitful

import Diversity: _calcabundance
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
"""
    Cache

Cache houses an integer array of moves made by all species in a timestep for the
update! function, `netmigration`.
"""
mutable struct Cache
  netmigration::Array{Int64, 2}
end

Lookup(df::DataFrame) = Lookup(df[:X], df[:Y], df[:Prob],
zeros(Float64, nrow(df)),zeros(Int64, nrow(df)))

function _mcmatch(m::AbstractMatrix, sim::SpeciesList, part::AbstractAbiotic)
    realm = _calcabundance(sim, m)
    return typematch(realm, sim, part) &&
    counttypes(sim) == size(realm, 1) &&
    countsubcommunities(part) == size(realm, 2)
end


"""
    tematch(sppl::SpeciesList, abenv::AbstractAbiotic)

Function to check that the types of a trait list and habitat list are
the same for a species list (`sppl`) and abiotic environment (`abenv`).
"""
function tematch(sppl::SpeciesList, abenv::AbstractAbiotic)
    (eltype(sppl.traits) == eltype(abenv.habitat)) &&
    (iscontinuous(sppl.traits) == iscontinuous(abenv.habitat))
end
"""
    trmatch(sppl::SpeciesList, traitrel::AbstractTraitRelationship)

Function to check that the types of a trait list and trait relationship list are
the same for a species list (`sppl`) and trait relationship (`traitrel`).
"""
function trmatch(sppl::SpeciesList, traitrel::AbstractTraitRelationship)
    eltype(sppl.traits) == eltype(traitrel) &&
    (iscontinuous(sppl.traits) == iscontinuous(traitrel))
end

importall Diversity.API
"""
    Ecosystem{Part <: AbstractAbiotic} <:
       AbstractMetacommunity{Float64, Matrix{Float64}, SpeciesList, Part}

Ecosystem houses information on species and their interaction with their
environment. For species, it holds abundances and locations, `abundances`,
as well as properties such as trait information, `spplist`, and movement types,
`lookup`. For environments, it provides information on environmental conditions
and available resources,`abenv`. Finally, there is a slot for the relationship
between the environment and the characteristics of the species, `relationship`.
"""
mutable struct Ecosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
    TR <: AbstractTraitRelationship} <:
   AbstractMetacommunity{Float64, Matrix{Int64}, Matrix{Float64}, SL, Part}
  abundances::GridLandscape
  spplist::SL
  abenv::Part
  ordinariness::Nullable{Matrix{Float64}}
  relationship::TR
  lookup::Vector{Lookup}
  cache::Cache

  function Ecosystem{Part, SL, TR}(abundances::GridLandscape,
    spplist::SL, abenv::Part, ordinariness::Nullable{Matrix{Float64}},
    relationship::TR, lookup::Vector{Lookup}, cache::Cache) where {Part <:
     AbstractAbiotic,
    SL <: SpeciesList, TR <: AbstractTraitRelationship}
    tematch(spplist, abenv) || error("Traits do not match habitats")
    trmatch(spplist, relationship) || error("Traits do not match trait functions")
    #_mcmatch(abundances.matrix, spplist, abenv) ||
    #  error("Dimension mismatch")
    new{Part, SL, TR}(abundances, spplist, abenv, ordinariness, relationship, lookup, cache)
  end
end
"""
    Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
        rel::AbstractTraitRelationship)

Function to create an `Ecosystem` given a species list, an abiotic environment
and trait relationship.
"""
function Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship)

  # Check there is enough energy to support number of individuals at set up
  all(getenergyusage(spplist) .<= getavailableenergy(abenv)) ||
    error("Environment does not have enough energy to support species")
  # Create matrix landscape of zero abundances
  ml = emptygridlandscape(abenv, spplist)
  # Populate this matrix with species abundances
  populate!(ml, spplist, abenv)
  # Create lookup table of all moves and their probabilities
  lookup_tab = genlookups(abenv.habitat, getkernel(spplist.movement))
  nm = zeros(Int64, size(ml.matrix))
  Ecosystem{typeof(abenv), typeof(spplist), typeof(rel)}(ml, spplist, abenv,
  Nullable{Matrix{Float64}}(), rel, lookup_tab, Cache(nm))
end

function _getabundance(eco::Ecosystem, input::Bool)
  relab = eco.abundances.matrix / sum(eco.abundances.matrix)
    return input ? relab : _calcabundance(eco.spplist, relab)
end

function _getmetaabundance(eco::Ecosystem)
  return sumoversubcommunities(eco, _getabundance(eco))
end
function _getpartition(eco::Ecosystem)
  return eco.abenv
end
function _gettypes(eco::Ecosystem)
    return eco.spplist
end

function _getordinariness!(eco::Ecosystem)
    if isnull(eco.ordinariness)
        eco.ordinariness = calcordinariness(eco.spplist, getabundance(eco, true))
    end
    get(eco.ordinariness)
end

"""
    gettraitrel(eco::Ecosystem)

Function to extract trait relationships.
"""
function gettraitrel(eco::Ecosystem)
  return eco.relationship
end

"""
    gethabitat(eco::Ecosystem)

Function to extract habitat from Ecosystem object.
"""
function gethabitat(eco::Ecosystem)
  return eco.abenv.habitat
end
"""
    getbudget(eco::Ecosystem)

Function to extract budget from Ecosystem object.
"""
function getbudget(eco::Ecosystem)
    return _getbudget(eco.abenv.budget)
end

"""
    getsize(eco::Ecosystem)

Function to extract size of habitat from Ecosystem object.
"""
function getsize(eco::Ecosystem)
  return _getsize(eco.abenv.habitat)
end

"""
    getgridsize(eco::Ecosystem)

Function to extract grid cell size of habitat from Ecosystem object.
"""
function getgridsize(eco::Ecosystem)
  return _getgridsize(eco.abenv.habitat)
end

"""
    getdimension(eco::Ecosystem)

Function to extract dimension of habitat from Ecosystem object.
"""
function getdimension(eco::Ecosystem)
    return _getdimension(eco.abenv.habitat)
end

"""
    getdispersaldist(eco::Ecosystem)

Function to extract average dispersal distance of species from Ecosystem object.
Returns a vector of distances, unless a specific species is provided as a String
or Integer.
"""
function getdispersaldist(eco::Ecosystem)
  dists = eco.spplist.movement.kernel.dist
  return dists
end
function getdispersaldist(eco::Ecosystem, spp::Int64)
  dists = eco.spplist.movement.kernel.dist
  return dists[spp]
end
function getdispersaldist(eco::Ecosystem, spp::String)
  num = find(eco.spplist.names.==spp)[1]
  getdispersaldist(eco, num)
end

"""
    getdispersalvar(eco::Ecosystem)

Function to extract dispersal varaince of species from Ecosystem object.
Returns a vector of distances, unless a specific species is provided as a String
or Integer.
"""
function getdispersalvar(eco::Ecosystem)
    vars = (eco.spplist.movement.kernel.dist).^2 .* pi ./ 4
    return vars
end
function getdispersalvar(eco::Ecosystem, spp::Int64)
    vars = (eco.spplist.movement.kernel.dist).^2 .* pi ./ 4
    return vars[spp]
end
function getdispersalvar(eco::Ecosystem, spp::String)
    num = find(eco.spplist.names.==spp)[1]
    getdispersalvar(eco, num)
end
"""
    resetrate!(eco::Ecosystem, rate::Quantity{Float64, typeof(𝐓^-1)})

Function to reset the rate of habitat change for a species.
"""
function resetrate!(eco::Ecosystem, rate::Quantity{Float64, typeof(𝐓^-1)})
    eco.abenv.habitat.change = HabitatUpdate{Unitful.Dimension{()}}(
    eco.abenv.habitat.change.changefun, rate)
end
function resetrate!(eco::Ecosystem, rate::Quantity{Float64, typeof(𝚯*𝐓^-1)})
    eco.abenv.habitat.change = HabitatUpdate{Unitful.Dimension{:Temperature}}(
    eco.abenv.habitat.change.changefun, rate)
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
 function _gaussian_disperse(r, Θ)
   t = r[1]
   exp(-(t/(1-t^2))^2) / (Θ^2 * π)
 end

 # Define gaussian kernel function
function _gaussian_disperse(r)
  exp(-((r[3]-r[1])^2+(r[4]-r[2])^2)) / π
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
  return map(r -> Lookup(_lookup(r, m, p, _gaussian_disperse)), relsize)
end

function _lookup(relSquareSize::Float64, maxGridSize::Int64,
                pThresh::Float64, dispersalfn::Function)
  # Create empty array
  lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])

   # Loop through directions until probability is below threshold
    k = 0
    m = 0
    count = 0
    while (k <= maxGridSize && m <= maxGridSize)
      count = count + 1
      calc_prob = pcubature(r -> dispersalfn(r),
                            [0, 0, k*relSquareSize, m*relSquareSize],
      [relSquareSize, relSquareSize, (k+1)*relSquareSize, (m+1)*relSquareSize],
       maxevals= 1000000)[1] / relSquareSize^2
      if m == 0 && calc_prob < pThresh
        break
      end
        if count == 1
          push!(lookup_tab, [k m calc_prob])
           k = k + 1
        else
          if (calc_prob > pThresh && m <= k)
            push!(lookup_tab, [k m calc_prob])
            m = m + 1
          else m = 0; k = k + 1
          end
        end
    end
    # If no probabilities can be calculated, threshold is too high
    nrow(lookup_tab) != 0 || error("probability threshold too high")
    # Find all other directions
    lookup_tab = _symmetric_grid(lookup_tab)
    #info(sum(lookup_tab[:, 3]))
    # Normalise
    lookup_tab[:Prob] = lookup_tab[:Prob]/sum(lookup_tab[:Prob])
    lookup_tab
end
