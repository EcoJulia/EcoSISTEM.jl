using Diversity
using Diversity.AbstractPartition
using Diversity.AbstractMetacommunity
using Diversity.psmatch
using Diversity.AbstractSimilarity

## Habitat types
abstract AbstractHabitat


# Habitats : matrix of float values
type Habitats <: AbstractHabitat
  matrix::Matrix{Float64}
  size::Real
end

# Niches : matrix of string values
type Niches <: AbstractHabitat
  matrix::Matrix{String}
  size::Real
end

# Env budget types
abstract AbstractBudget

type Budget <: AbstractBudget
  matrix::Matrix{Float64}
end

# Species trait types
abstract AbstractTraits

type StringTraits <: AbstractTraits
  traits::Vector{String}
end
type RealTraits <: AbstractTraits
  traits::Vector{Real}
end

type TraitRelationship
  matrix::Matrix{Real}
end

# Species energy types
abstract AbstractEnergy

type RealEnergy <: AbstractEnergy
  energy::Vector{Real}
end

abstract AbstractMovement

type GaussianMovement <: AbstractMovement
  var::Vector{Real}
  thresh::Float64
end

function GaussianMovement(move_var, numSpecies, pThresh)
  GaussianMovement(repmat([move_var], numSpecies), pThresh)
end

# Species list type - all info on species
type SpeciesList{FP, M <: AbstractMatrix, T <: AbstractTraits, A <: AbstractVector,
                 E<: AbstractEnergy, TR <: Tree, MO <: AbstractMovement} <:
                             AbstractSimilarity{FP, M}
  similarity::M
  traits::T
  abun::A
  energy::E
  phylo::TR
  movement::MO
end

function SpeciesList{FP <: AbstractFloat}(M::AbstractMatrix{FP}, T::AbstractTraits,
                      A::AbstractVector, E::AbstractEnergy, TR::Tree, MO::AbstractMovement)
  SpeciesList{FP, typeof(M), typeof(T), typeof(A), typeof(E), typeof(TR), typeof(MO)}(M, T, A, E, TR, MO)
end

function SpeciesList(NumberSpecies::Int64, NumberTraits::Int64,
                      abun_dist::Distribution, energy::AbstractVector,
                      movement::GaussianMovement)
  # Create tree
  tree = jcoal(NumberSpecies, 100)
  # Create traits and assign to tips
  trts = map(string, 1:NumberTraits)
  assign_traits!(tree, 0.5, trts)
  # Get traits from tree
  sp_trt = get_traits(tree, true)
  # Create similarity matrix (for now identity)
  similarity = eye(NumberSpecies)
  # Draw random set of abundances from distribution
  abun = rand(abun_dist)
  # error out when abun dist and NumberSpecies are not the same (same for energy dist)
  length(abun)==NumberSpecies || throw(DimensionMismatch("Abundance vector
                                        doesn't match number species"))
  length(energy)==NumberSpecies || throw(DimensionMismatch("Energy vector
                                        doesn't match number species"))
  size(similarity)==(NumberSpecies,NumberSpecies) || throw(DimensionMismatch("
                              Similarity matrix doesn't match number species"))

  SpeciesList(similarity, StringTraits(sp_trt), abun, RealEnergy(energy), tree, movement)
end


# Abiotic environment types- all info about habitat and relationship to species
# traits
abstract AbstractAbiotic{H<: AbstractHabitat, B<:AbstractBudget}

type MatrixAbioticEnv{H, B} <: AbstractAbiotic{H, B}
  habitat::H
  budget::B
end

function MatrixAbioticEnv(NumberNiches::Int64, dimension::Tuple, maxBud::Real)
  # Create niches
  niches = map(string, 1:NumberNiches)
  # Create niche-like environment
  hab = random_habitat(dimension, niches, 0.5, repmat([0.5], NumberNiches))
  # Create empty budget and for now fill with one value
  bud = zeros(dimension)
  fill!(bud, maxBud)

  MatrixAbioticEnv(Niches(hab), Budget(bud))
end

# Matrix Landscape types - houses abundances (initially empty)
abstract AbstractStructuredPartition{A} <: AbstractPartition{Float64, A}

type MatrixLandscape{A} <: AbstractStructuredPartition{A}
  abundances::A
end

function MatrixLandscape(abenv::AbstractAbiotic, spplist::SpeciesList)
  # Create an array of zero abundances according to the size of the habitat
  # and the number of species
  abundances=zeros(length(spplist.abun),size(abenv.habitat.matrix,1),
                   size(abenv.habitat.matrix,2))
  MatrixLandscape(abundances)
end

# Ecosystem type - holds all information and populates ML
abstract AbstractEcosystem{A, Part <: AbstractStructuredPartition,
          S <: SpeciesList, AB <: AbstractAbiotic, R<: TraitRelationship
          } <: AbstractMetacommunity{Float64, A, Part}

type Ecosystem{A, Part, S, AB, R} <: AbstractEcosystem{A, Part, S, AB, R}
  partition::Part
  ordinariness::Nullable{A}
  spplist::S
  abenv::AB
  relationship::R
end

function Ecosystem(spplist::SpeciesList, abenv::AbstractAbiotic, traits::Bool)
  sum(spplist.abun.*spplist.energy.energy)<= sum(abenv.budget.matrix) ||
  error("Environment does not have enough energy to support species")
  # Create matrix landscape of zero abundances
  ml = MatrixLandscape(abenv, spplist)
  # Populate this matrix with species abundances
  species = length(spplist.abun)
  populate!(ml, spplist, abenv, traits)
  # For now create an identity matrix for the species relationships
  A = typeof(ml.abundances)
  rel = eye(length(spplist.traits.traits), size(abenv.habitat.matrix,3))

  Ecosystem(ml, Nullable{A}(), spplist, abenv, TraitRelationship(rel))
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
  # Create a new ecosystem with the same components
  Ecosystem(ml, Nullable{A}(), spplist, abenv, TraitRelationship(rel))
end
