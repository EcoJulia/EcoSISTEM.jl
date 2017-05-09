using Diversity
using Cubature
using DataFrames
using Unitful

importall Diversity.API


type Ecosystem{FP, A, Sim <: SpeciesList, Part <: AbstractAbiotic,
   Rel <: TraitRelationship, Look <:AbstractArray, G <: GridLandscape} <:
   AbstractMetacommunity{AbstractFloat, A, Sim, Part}
  abundances::G
  spplist::Sim
  abenv::Part
  ordinariness::Nullable{A}
  relationship::Rel
  lookup::Look
end


function Ecosystem(spplist::SpeciesList, abenv::MatrixAbioticEnv, traits::Bool)

  # Check there is enough energy to support number of individuals at set up
  sum(spplist.abun.*spplist.energy.energy)<= sum(abenv.budget.matrix) ||
  error("Environment does not have enough energy to support species")
  # Create matrix landscape of zero abundances
  ml = emptygridlandscape(abenv, spplist)
  # Populate this matrix with species abundances
  populate!(ml, spplist, abenv, traits)
  # For now create an identity matrix for the species relationships
  rel = eye(length(spplist.traits.traits))
  # Create lookup table of all moves and their probabilities
  lookup_tab = map(x -> lookup(abenv.habitat.size, maximum(size(abenv.habitat.matrix)),
   x, spplist.movement.thresh), spplist.movement.var)
   A = typeof(ml.matrix)
  Ecosystem{AbstractFloat, A, SpeciesList, AbstractAbiotic,
  TraitRelationship, typeof(lookup_tab), GridLandscape}(ml, spplist, abenv, Nullable{A}(), TraitRelationship(rel),
   lookup_tab)
end

function _getabundance(eco::Ecosystem)
  return eco.abundances.matrix
end
function _getmetaabundance(eco::Ecosystem)
  return eco.abundances.matrix
end
function _getpartition(eco::Ecosystem)
  return eco.abenv
end
function _gettypes(eco::Ecosystem)
    return eco.spplist
end

function _getordinariness!(eco::Ecosystem)
    if isnull(eco.ordinariness)
        eco.ordinariness = _calcordinariness(_gettypes(eco), _getabundance(eco))
    end
    get(eco.ordinariness)
end

# Function to copy an Ecosystem type for modification - similar to array copy
function copy_eco(eco::Ecosystem)
  # Collects species list and abiotic environment
  spplist = eco.spplist
  abenv = eco.abenv
  # Create new ML with abundances from ecosystem
  ml= MatrixLandscape(abenv, spplist)
  ml.abundances = copy(eco.partition.abundances)
  # Create relationships same as before
  A = typeof(ml.abundances)
  rel = eye(length(spplist.traits.traits), size(abenv.habitat.matrix,3))
  lookup_tab = eco.lookup
  # Create a new ecosystem with the same components
  Ecosystem(ml, Nullable{A}(), spplist, abenv, TraitRelationship(rel), lookup)
end

function symmetric_grid(grid::DataFrame)
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


function lookup(squareSize::Float64, maxGridSize::Int64, dispersalSD::Float64,
                pThresh::Float64)
  # Create empty array
  lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])
  # Define gaussian kernel function
   function disperse(r::AbstractArray)
     1/(π * dispersalSD^2)*exp(-((r[3]-r[1])^2+(r[4]-r[2])^2)/(dispersalSD^2))
   end
   # Loop through directions until probability is below threshold
    k=0; m=0; count = 0
    while (k<= maxGridSize && m <=maxGridSize)
      count= count + 1
      calc_prob = pcubature(disperse, [0, 0, k*squareSize, m*squareSize],
      [squareSize, squareSize, (k+1)*squareSize, (m+1)*squareSize],
       maxevals= 1000000)[1] / squareSize^2
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
    lookup_tab = symmetric_grid(lookup_tab)
    #info(sum(lookup_tab[:, 3]))
    # Normalise
    lookup_tab[:Prob] = lookup_tab[:Prob]/sum(lookup_tab[:Prob])
    lookup_tab
end
