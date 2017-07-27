using Diversity
using Cubature
using DataFrames
using Unitful

struct Lookup
  x::Vector{Int64}
  y::Vector{Int64}
  p::Vector{Float64}
end

Lookup(df::DataFrame) = Lookup(df[:X], df[:Y], df[:Prob])

function _mcmatch(m::AbstractMatrix, sim::SpeciesList, part::AbstractAbiotic)
    realm = _calcabundance(sim, m)
    return typematch(realm, sim, part) &&
    counttypes(sim) == size(realm, 1) &&
    countsubcommunities(part) == size(realm, 2)
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
mutable struct Ecosystem{Part <: AbstractAbiotic, SL <: SpeciesList} <:
   AbstractMetacommunity{Float64, Matrix{Float64}, SL, Part}
  abundances::GridLandscape
  spplist::SL
  abenv::Part
  ordinariness::Nullable{Matrix{Float64}}
  relationship::TraitRelationship
  lookup::Vector{Lookup}

  function Ecosystem{Part, SL}(abundances::GridLandscape,
    spplist::SL, abenv::Part, ordinariness::Nullable{Matrix{Float64}},
    relationship::TraitRelationship, lookup::Vector{Lookup}) where {Part <:
     AbstractAbiotic,
    SL <: SpeciesList}
    eltype(abenv.budget) == eltype(spplist.requirement) ||
      error("Environment and species energy not of the same type")
    #_mcmatch(abundances.matrix, spplist, abenv) ||
    #  error("Dimension mismatch")
    new{Part, SL}(abundances, spplist, abenv, ordinariness, relationship, lookup)
  end
end
"""
    Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv, traits)

Function to create an ecosystem given a species list and an abiotic environment.
When the landscape is populated, there is the option for the individuals to be
assigned locations based on their trait preferences.
"""
function Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
   rel::TraitRelationship)

  # Check there is enough energy to support number of individuals at set up
  sum(spplist.abun .* spplist.requirement.energy) <= sum(abenv.budget.matrix) ||
  error("Environment does not have enough energy to support species")
  # Create matrix landscape of zero abundances
  ml = emptygridlandscape(abenv, spplist)
  # Populate this matrix with species abundances
  populate!(ml, spplist, abenv)
  # Create lookup table of all moves and their probabilities
  lookup_tab = genlookups(abenv.habitat, getkernel(spplist.movement))

  Ecosystem{typeof(abenv), typeof(spplist)}(ml, spplist, abenv,
  Nullable{Matrix{Float64}}(), rel, lookup_tab)
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

function gethabitat(eco::Ecosystem, pos::Int64)
  x, y = convert_coords(pos, size(gethabitat(eco).matrix, 1))
  return gethabitat(eco).matrix[x, y]
end

function getenvtype(eco::Ecosystem)
  return typeof(eco.abenv.habitat)
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

"""
    genlookups(hab::AbstractHabitat, mov::GaussianMovement)

Function to generate lookup tables, which hold information on the probability
of moving to neighbouring squares.
"""
function genlookups(hab::AbstractHabitat, mov::GaussianKernel)
  relsize =  hab.size ./ mov.var
  m = maximum(size(hab.matrix))
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
