using Diversity
using Cubature
using DataFrames
using Unitful

mutable struct Lookup
  x::Vector{Int64}
  y::Vector{Int64}
  p::Vector{Float64}
  pnew::Vector{Float64}
  moves::Vector{Int64}
end

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

function tematch(sppl::SpeciesList, abenv::AbstractAbiotic)
    (eltype(sppl.traits) == eltype(abenv.habitat)) &&
    (iscontinuous(sppl.traits) == iscontinuous(abenv.habitat))
end
function trmatch(sppl::SpeciesList, traitrel::AbstractTraitRelationship)
    eltype(sppl.traits) == eltype(traitrel)
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
   AbstractMetacommunity{Float64, Matrix{Float64}, SL, Part}
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
    Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv, traits)

Function to create an ecosystem given a species list and an abiotic environment.
When the landscape is populated, there is the option for the individuals to be
assigned locations based on their trait preferences.
"""
function Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship)

  # Check there is enough energy to support number of individuals at set up
  #sum(spplist.abun .* spplist.requirement.energy) <= sum(abenv.budget.matrix) ||
  #error("Environment does not have enough energy to support species")
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
    return input ? relab : calcabundance(eco.spplist, relab)
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

function gettraitrel(eco::Ecosystem)
  return eco.relationship
end
function gethabitat(eco::Ecosystem)
  return eco.abenv.habitat
end

function getenvtype(eco::Ecosystem)
  return typeof(eco.abenv.habitat)
end

function getsize(eco::Ecosystem)
  return _getsize(eco.abenv.habitat)
end

function getgridsize(eco::Ecosystem)
  return _getgridsize(eco.abenv.habitat)
end
function getdimension(eco::Ecosystem)
    return _getdimension(eco.abenv.habitat)
end

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
