using Diversity
using Diversity.AbstractPartition
using Diversity.AbstractMetacommunity
using Diversity.psmatch


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
abstract AbstractSpeciesList{A, N, B, T <: AbstractTraits, Sim, E<: AbstractEnergy}

type SpeciesList{A, N, B, T, Sim, E} <: AbstractSpeciesList{A, N, B, T, Sim, E}
  traits::T
  similarity::Sim
  abun::Vector{A}
  energy::E
  phylo::Tree{N,B}
end

function SpeciesList{T, Sim, E}(traits::T, similarity::Sim, abun, energy::E, phylo)
  SpeciesList{A, N, B, T, Sim, E}(traits, similarity, abun, energy, phylo)
end

function SpeciesList(NumberSpecies::Int64, NumberTraits::Int64,
                      abun_dist::Distribution, energy_dist::Distribution)
  tree = jcoal(NumberSpecies, 100)
  trts = map(string, 1:NumberTraits)
  assign_traits!(tree, 0.2, trts)
  sp_trt = get_traits(tree, true)
  similarity = eye(NumberSpecies)
  abun = rand(abun_dist)
  energy = rand(energy_dist)
  SpeciesList(StringTraits(sp_trt), similarity, abun, RealEnergy(energy), tree)
end


# Abiotic environment types- all info about habitat and relationship to species
# traits
abstract AbstractAbiotic{H<: AbstractHabitat, R<: TraitRelationship, B<:AbstractBudget}

type AbioticEnv{H, R, B} <: AbstractAbiotic{H, R, B}
  habitat::H
  relationship::R
  budget::B
end
function AbioticEnv(NumberNiches::Int64, dimension::Tuple,
                    spplist::AbstractSpeciesList)
  niches = map(string, 1:NumberNiches)
  hab = random_habitat(dimension, niches, 0.5, [0.5,0.5])
  rel = eye(length(spplist.traits.traits), NumberNiches)
  bud = zeros(dimension)
  fill!(bud, 100)
  AbioticEnv(Niches(hab), TraitRelationship(rel), Budget(bud))
end

# Matrix Landscape types - houses abundances (initially empty)
abstract AbstractStructuredPartition{A} <: AbstractPartition{Float64, A}

type MatrixLandscape{A} <: AbstractStructuredPartition{A}
  abundances::A
end


function MatrixLandscape(abenv::AbstractAbiotic, spplist::AbstractSpeciesList)
  abundances=zeros(size(abenv.habitat.matrix,1),size(abenv.habitat.matrix,2),
             length(spplist.abun))
  MatrixLandscape(abundances)
end

# Ecosystem type - holds all information and populates ML
abstract AbstractEcosystem{A, Part <: AbstractStructuredPartition,
          S <: AbstractSpeciesList, AB <: AbstractAbiotic} <:
                AbstractMetacommunity{Float64, A, Part}

type Ecosystem{A, Part, S, AB} <: AbstractEcosystem{A, Part, S, AB}
  partition::Part
  ordinariness::Nullable{A}
  ssplist::S
  abenv::AB
end

function Ecosystem(spplist::AbstractSpeciesList, abenv::AbstractAbiotic)
  ml = MatrixLandscape(abenv, spplist)
  species = length(spplist.abun)
  populate!(ml, spplist, abenv)
  A = typeof(ml.abundances)
  Ecosystem(ml, Nullable{A}(), spplist, abenv)
end
