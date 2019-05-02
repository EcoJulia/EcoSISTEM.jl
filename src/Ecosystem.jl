using Diversity
if VERSION > v"1.0.0"
    using HCubature
else
    using Cubature
end
using DataFrames
using Unitful
using MyUnitful
using Missings
using Compat
using RecipesBase

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
  totalE::Matrix{Float64}
  valid::Bool
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

"""
    AbstractEcosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
        TR <: AbstractTraitRelationship} <: AbstractMetacommunity{Float64,
            Matrix{Int64}, Matrix{Float64}, SL, Part}

Abstract supertype for all ecosystem types and a subtype of AbstractMetacommunity.
"""
abstract type
    AbstractEcosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
        TR <: AbstractTraitRelationship} <: AbstractMetacommunity{Float64, Matrix{Int64},
                                        Matrix{Float64}, SL, Part}
end
"""
    Ecosystem{Part <: AbstractAbiotic} <:
       AbstractEcosystem{Part, SL, TR}

Ecosystem houses information on species and their interaction with their
environment. For species, it holds abundances and locations, `abundances`,
as well as properties such as trait information, `spplist`, and movement types,
`lookup`. For environments, it provides information on environmental conditions
and available resources,`abenv`. Finally, there is a slot for the relationship
between the environment and the characteristics of the species, `relationship`.
"""
mutable struct Ecosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
    TR <: AbstractTraitRelationship} <: AbstractEcosystem{Part, SL, TR}
  abundances::GridLandscape
  spplist::SL
  abenv::Part
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::TR
  lookup::Vector{Lookup}
  cache::Cache

  function Ecosystem{Part, SL, TR}(abundances::GridLandscape,
    spplist::SL, abenv::Part, ordinariness::Union{Matrix{Float64}, Missing},
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

@recipe function f(::AbstractMovement, eco::Ecosystem, spp::Int64)
    l = eco.lookup[spp]
    maxX = maximum(l.x)
    maxY = maximum(l.y)
    x, y = round(Int64, maxX/2), round(Int64, maxY/2)
    # Can't go over maximum dimension
    valid = findall((l.x .> -x) .& (l.y .> -y) .&
     (l.x .<= (maxX - x)) .& (l.y .<= (maxY - y)))
    probs = l.p[valid]
    probs ./= sum(probs)
    xs = (l.x[valid] .+ x)
    ys = (l.y[valid] .+ y)
    A = zeros(maxX, maxY)
    for i in eachindex(xs)
      A[xs[i], ys[i]] .= probs[i]
    end
    seriestype  :=  :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> "Movement kernel (km)"
    xrange(gethabitat(eco)), yrange(gethabitat(eco)), A
end
"""
    Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
        rel::AbstractTraitRelationship)

Function to create an `Ecosystem` given a species list, an abiotic environment
and trait relationship.
"""
function Ecosystem(popfun::Function, spplist::SpeciesList{T, Req}, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship) where {T, Req}

  # Check there is enough energy to support number of individuals at set up
  #all(getenergyusage(spplist) .<= getavailableenergy(abenv)) ||
    #error("Environment does not have enough energy to support species")
  # Create matrix landscape of zero abundances
  ml = emptygridlandscape(abenv, spplist)
  # Populate this matrix with species abundances
  popfun(ml, spplist, abenv)
  # Create lookup table of all moves and their probabilities
  lookup_tab = genlookups(abenv.habitat, getkernel(spplist.movement))
  nm = zeros(Int64, size(ml.matrix))
  totalE = zeros(Float64, (size(ml.matrix, 2), numrequirements(Req)))
  Ecosystem{typeof(abenv), typeof(spplist), typeof(rel)}(ml, spplist, abenv,
  missing, rel, lookup_tab, Cache(nm, totalE, false))
end

function Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship)
   return Ecosystem(populate!, spplist, abenv, rel)
end
GLOBAL_typedict["Ecosystem"] = Ecosystem

function addspecies!(eco::Ecosystem, abun::Int64)
    eco.abundances.matrix = vcat(eco.abundances.matrix, zeros(size(eco.abundances.matrix, 2)))
    eco.abundances.grid = reshape(eco.abundances.matrix, (counttypes(eco.spplist, true)+1, _getdimension(eco.abenv.habitat)...))
    repopulate!(eco, abun)
    push!(eco.spplist.names, string.(counttypes(eco.spplist, true)+1))
    append!(eco.spplist.abun, abun)
    append!(eco.spplist.native, true)
    addtraits!(eco.spplist.traits)
    addmovement!(eco.spplist.movement)
    addparams!(eco.spplist.params)
    addrequirement!(eco.spplist.requirement)
    addtypes!(eco.spplist.types)
end
function addtraits!(tr::GaussTrait)
    append!(tr.mean, tr.mean[end])
    append!(tr.var, tr.var[end])
end

addmovement!(mv::AbstractMovement) = append!(mv.kernel.dist, mv.kernel.dist[end])

function addparams!(pr::AbstractParams)
    append!(pr.birth, pr.birth[end])
    append!(pr.death, pr.death[end])
end

addrequirement!(rq::AbstractRequirement) = append!(rq.energy, rq.energy[end])

function addtypes!(ut::UniqueTypes)
    ut = UniqueTypes(ut.num+1)
end

"""
    CachedEcosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
        TR <: AbstractTraitRelationship} <: AbstractEcosystem{Part, SL, TR}

CachedEcosystem houses the same information as Ecosystem (see ?Ecosystem), but
holds the time period abundances as a CachedGridLandscape, so that they may
be present or missing.
"""
mutable struct CachedEcosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
    TR <: AbstractTraitRelationship} <: AbstractEcosystem{Part, SL, TR}
  abundances::CachedGridLandscape
  spplist::SL
  abenv::Part
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::TR
  lookup::Vector{Lookup}
  cache::Cache
end

"""
    CachedEcosystem(eco::Ecosystem, outputfile::String, rng::StepRangeLen)

Function to create a CachedEcosystem given an existing ecosystem, `eco`,
output folder to which the simulations are saved, `outputfile`, and a range of
times over which to simulate, `rng`.
"""
function CachedEcosystem(eco::Ecosystem, outputfile::String, rng::StepRangeLen)
    size(eco.abenv.habitat, 3) == length(rng) || error("Time range does not match habitat")
    abundances = CachedGridLandscape(outputfile, rng)
    abundances.matrix[1] = eco.abundances
  CachedEcosystem{typeof(eco.abenv), typeof(eco.spplist), typeof(eco.relationship)}(abundances,
  eco.spplist, eco.abenv, eco.ordinariness, eco.relationship, eco.lookup, eco.cache)
end

GLOBAL_typedict["CachedEcosystem"] = CachedEcosystem

import Diversity.API: _getabundance
function _getabundance(eco::Ecosystem, input::Bool)
    if input
        return eco.abundances.matrix
    else
        return _calcabundance(_gettypes(eco), eco.abundances.matrix / sum(eco.abundances.matrix))[1]
    end
end


function _getabundance(cache::CachedEcosystem, input::Bool)
    if all(ismissing.(cache.abundances.matrix))
        error("Abundances are missing")
    else
        id = Compat.findall(.!ismissing.(cache.abundances.matrix))[end]
        abun = cache.abundances.matrix[id]
    end
    if input
        return abun.matrix
    else
        return abun.matrix / sum(abun.matrix)
    end
end
import Diversity.API: _getmetaabundance
function _getmetaabundance(eco::Ecosystem)
  return sumoversubcommunities(eco, _getabundance(eco))
end

function _getmetaabundance(eco::CachedEcosystem)
  return sumoversubcommunities(eco, _getabundance(eco))
end

import Diversity.API: _getpartition
function _getpartition(eco::Union{Ecosystem, CachedEcosystem})
  return eco.abenv
end
import Diversity.API: _gettypes
function _gettypes(eco::Union{Ecosystem, CachedEcosystem})
    return eco.spplist
end
import Diversity.API: _getordinariness!
function _getordinariness!(eco::Union{Ecosystem, CachedEcosystem})
    if ismissing(eco.ordinariness)
        relab = getabundance(eco, false)
        eco.ordinariness = _calcordinariness(eco.spplist, relab)
    end
    return eco.ordinariness
end

import Diversity.API._getscale
function _getscale(eco::Ecosystem)
    return _calcabundance(_gettypes(eco), getabundance(eco, false))[2]
end

function invalidatecaches!(eco::Union{Ecosystem, CachedEcosystem})
    eco.ordinariness = missing
    eco.cache.netmigration .= 0
    eco.cache.valid = false
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
  num = Compat.findall(eco.spplist.names.==spp)[1]
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
    num = Compat.findall(eco.spplist.names.==spp)[1]
    getdispersalvar(eco, num)
end
"""
    resetrate!(eco::Ecosystem, rate::Quantity{Float64, typeof(𝐓^-1)})

Function to reset the rate of habitat change for a species.
"""
function resetrate!(eco::Ecosystem, rate::Quantity{Float64, typeof(𝐓^-1)})
    eco.abenv.habitat.change = HabitatUpdate{Unitful.Dimensions{()}}(
    eco.abenv.habitat.change.changefun, rate)
end
function resetrate!(eco::Ecosystem, rate::Quantity{Float64, typeof(𝚯*𝐓^-1)})
    eco.abenv.habitat.change = HabitatUpdate{typeof(dimension(1K))}(
    eco.abenv.habitat.change.changefun, rate)
end
function resetrate!(eco::Ecosystem, rate::Quantity{Float64, 𝐓^-1})
    eco.abenv.habitat.change = HabitatUpdate{Unitful.Dimensions{()}}(
    eco.abenv.habitat.change.changefun, rate)
end
function resetrate!(eco::Ecosystem, rate::Quantity{Float64, 𝚯*𝐓^-1})
    eco.abenv.habitat.change = HabitatUpdate{typeof(dimension(1K))}(
    eco.abenv.habitat.change.changefun, rate)
end

function resettime!(eco::Ecosystem)
    _resettime!(eco.abenv.habitat)
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

function _2Dt_disperse(r, a, b)
    newfun = function(r) return((b - 1)/(π * a^2)) * (1 + ((r[3]-r[1])^2+(r[4]-r[2])^2)/a^2)^-b end
    return newfun
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
function genlookups(hab::AbstractHabitat, mov::LongTailKernel)
  relsize = _getgridsize(hab) ./ km
  m = maximum(_getdimension(hab))
  p = mov.thresh
  a = mov.dist
  b = mov.shape
  return map(r -> Lookup(_lookup(r, m, p, _2Dt_disperse)), relsize)
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
