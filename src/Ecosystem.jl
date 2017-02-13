using Diversity
using Diversity.AbstractPartition
using Diversity.AbstractMetacommunity
using Diversity.psmatch
using Diversity.AbstractSimilarity

## Habitat types
abstract AbstractHabitat

type Habitats <: AbstractHabitat
  matrix::Matrix{Float64}
end
type Niches <: AbstractHabitat
  matrix::Matrix{String}
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

# Species list type - all info on species

type SpeciesList{FP, M <: AbstractMatrix, A <: AbstractVector, N<: Any, B <: Any,
                             T <: AbstractTraits, E<: AbstractEnergy} <:
                             AbstractSimilarity{FP, M}
  similarity::M
  traits::T
  abun::A
  energy::E
  phylo::Tree{N,B}
end

function SpeciesList{FP <: AbstractFloat, A <: AbstractVector, N <: Any, B <: Any,
          T <: AbstractTraits, E<: AbstractEnergy}(similarity::AbstractMatrix{FP},
          traits::T, abun::A, energy::E, phylo::Tree{N,B})
  SpeciesList{FP, typeof(similarity), A, N, B, T, E}(similarity, traits, abun,
                                                     energy, tree)
end
function SpeciesList(NumberSpecies::Int64, NumberTraits::Int64,
                      abun_dist::Distribution, energy::AbstractVector)
  # error out when abun dist and NumberSpecies are not the same (same for energy dist)
  tree = jcoal(NumberSpecies, 100)
  trts = map(string, 1:NumberTraits)
  assign_traits!(tree, 0.2, trts)
  sp_trt = get_traits(tree, true)
  similarity = eye(NumberSpecies)
  abun = rand(abun_dist)
  length(abun)==NumberSpecies || throw(DimensionMismatch("Abundance vector
                                        doesn't match number species"))
  length(energy)==NumberSpecies || throw(DimensionMismatch("Energy vector
                                        doesn't match number species"))
  size(similarity)==(NumberSpecies,NumberSpecies) || throw(DimensionMismatch("
                              Similarity matrix doesn't match number species"))
  SpeciesList(similarity, StringTraits(sp_trt), abun,
                            RealEnergy(energy), tree)
end
# Class of Tree traits eventually, but for now just pass in a fully annotated tree.

# Abiotic environment types- all info about habitat and relationship to species
# traits
abstract AbstractAbiotic{H<: AbstractHabitat, B<:AbstractBudget}

type MatrixAbioticEnv{H, B} <: AbstractAbiotic{H, B}
  habitat::H
  budget::B
end

function MatrixAbioticEnv(NumberNiches::Int64, dimension::Tuple, maxBud::Real)
  niches = map(string, 1:NumberNiches)
  hab = random_habitat(dimension, niches, 0.5, repmat([0.5], NumberNiches))
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

function Ecosystem(spplist::SpeciesList, abenv::AbstractAbiotic)
  ml = MatrixLandscape(abenv, spplist)
  species = length(spplist.abun)
  populate!(ml, spplist, abenv)
  A = typeof(ml.abundances)
  rel = eye(length(spplist.traits.traits), size(abenv.habitat.matrix,3))
  Ecosystem(ml, Nullable{A}(), spplist, abenv, TraitRelationship(rel))
end

function copy_eco(eco::Ecosystem)
  spplist = eco.spplist
  abenv = eco.abenv
  ml= MatrixLandscape(abenv, spplist)
  ml.abundances = copy(eco.partition.abundances)
  A = typeof(ml.abundances)
  rel = eye(length(spplist.traits.traits), size(abenv.habitat.matrix,3))
  Ecosystem(ml, Nullable{A}(), spplist, abenv, TraitRelationship(rel))
end
# Abstract species list subclass of abstract similarity

tree = jcoal(2, 100)
trts = map(string, 1:2)
assign_traits!(tree, 0.2, trts)
