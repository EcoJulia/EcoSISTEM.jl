using Diversity
using Cubature
using DataFrames
using Unitful

importall Diversity.API



# Species trait types
abstract AbstractTraits

type StringTraits <: AbstractTraits
  traits::Vector{String}
end
type RealTraits <: AbstractTraits
  traits::Vector{Real}
end

# Species energy types
abstract AbstractEnergy

type RealEnergy <: AbstractEnergy
  energy::Vector{Real}
end


# Species list type - all info on species
type SpeciesList{FP <: AbstractFloat, M <: AbstractMatrix, T <: AbstractTraits,
                 E<: AbstractEnergy, TR <: AbstractTree, MO <: AbstractMovement} <:
                             AbstractTypes
  similarity::M
  names::Vector{String}
  traits::T
  abun::Vector{Int64}
  energy::E
  phylo::TR
  movement::MO
end

function SpeciesList{FP <: AbstractFloat, T <: AbstractTraits,
    E<: AbstractEnergy, TR <: AbstractTree, MO <: AbstractMovement}(similarity::AbstractMatrix{FP},
    names:: Vector{String}, traits::T, abun::Vector{Int64}, energy::E, phylo::TR,
    movement::MO)
    # Assign names
    names = map(x -> "$x", 1:size(similarity, 1))

    # Check similarity is square matrix
    size(similarity, 1) == size(similarity, 2) ||
    throw(DimensionMismatch("Similarity matrix is not square"))

    # Check dimensions of abundance and similarity match
    length(abun) == size(similarity, 1) ||
    throw(DimensionMismatch("Similarity matrix does not match abundances"))

    # Check similarity is bounded between 0 and 1
    minimum(similarity) ≥ 0 || throw(DomainError())
    maximum(similarity) ≤ 1 || warn("Similarity matrix has values above 1")
    SpeciesList{FP, typeof(similarity), T, E, TR, MO}(similarity, names, traits, abun, energy, phylo, movement)
end

function SpeciesList{E<: AbstractEnergy,
    MO <: AbstractMovement}(NumberSpecies::Int64,
    NumberTraits::Int64, abun_dist::Distribution, energy::E,
    movement::MO)

    names = map(x -> "$x", 1:NumberSpecies)
    # Create tree
    tree = jcoal(names, 100)
    # Create traits and assign to tips
    trts = map(string, 1:NumberTraits)
    assign_traits!(tree, 0.5, trts)
    # Get traits from tree
    sp_trt = StringTraits(get_traits(tree, names, true)[:,1])
    # Create similarity matrix (for now identity)
    similarity = eye(NumberSpecies)
    # Draw random set of abundances from distribution
    abun = rand(abun_dist)
    # error out when abun dist and NumberSpecies are not the same (same for energy dist)
    length(abun)==NumberSpecies || throw(DimensionMismatch("Abundance vector
                                          doesn't match number species"))
    length(energy.energy)==NumberSpecies || throw(DimensionMismatch("Energy vector
                                          doesn't match number species"))
    size(similarity)==(NumberSpecies,NumberSpecies) || throw(DimensionMismatch("
                                Similarity matrix doesn't match number species"))

  SpeciesList(similarity, names, sp_trt, abun, energy, tree, movement)
end

function _counttypes(sl::SpeciesList)
    return size(sl.similarity, 1)
end

function _calcsimilarity(sl::SpeciesList)
    return sl.similarity
end

function _floattypes{FP, M}(::SpeciesList{FP, M})
    return Set([FP])
end

function _getnames(sl::SpeciesList)
    return sl.names
end

function _calcordinariness(sl::SpeciesList, a::AbstractArray)
    _calcsimilarity(sl) * a
end


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
