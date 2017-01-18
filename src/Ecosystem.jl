using Diversity
using Diversity.AbstractPartition
using Diversity.AbstractMetacommunity
using Diversity.psmatch

abstract AbstractHabitat

type Habitats <: AbstractHabitat
  matrix::Matrix{Float64}
end
type Niches <: AbstractHabitat
  matrix::Matrix{String}
end

abstract AbstractBudget

type Budget <: AbstractBudget
  matrix::Matrix{Float64}
end

#abstract AbstractStructuredPartition{A, H <: AbstractHabitat, B<: AbstractBudget} <:
# AbstractPartition{Float64, A}

#type MatrixLandscape{A, H, B} <: AbstractStructuredPartition{A, H, B}
#  abundances::A
#  habitat::H
#  budget::B
#end

#function MatrixLandscape{A, H, B}(abundances::A, habitats::H, budget::B)
#  dima = size(abundances)
#  dimh = size(habitats.matrix)
#  dimb = size(budget.matrix)
#  length(dima) == 3 || error("Abundances must be three dimensional!")
#  ((dima[2] == dimh[1]) && (dima[3] == dimh[2])) ||
#   error("Dimension mismatch between abundances and habitats")
#  ((dimh[1] == dimb[1]) && (dimh[2] == dimb[2])) ||
#   error("Dimension mismatch between budgets and habitats")
#  relative = abundances #/ sum(abundances)
#  MatrixLandscape{typeof(relative), H, B}(relative, habitats, budget)
#end

abstract AbstractTraits

type StringTraits <: AbstractTraits
  traits::Vector{String}
end
type RealTraits <: AbstractTraits
  traits::Vector{Real}
end

abstract AbstractEnergy

type RealEnergy <: AbstractEnergy
  energy::Vector{Real}
end


#abstract AbstractEcosystem{A, Part <: AbstractStructuredPartition, Sim, T <: AbstractTraits, E <: AbstractEnergy} <:
# AbstractMetacommunity{Float64, A, Part, Sim}

#type Ecosystem{A, Part, Sim, T, E} <: AbstractEcosystem{A, Part, Sim, T, E}
#  partition::Part
#  similarity::Sim
#  ordinariness::Nullable{A}
#  traits::T
#  energy::E
#end

#function Ecosystem{Part, Sim, T, E}(partition::Part, similarity::Sim, traits::T, energy::E)
#  psmatch(partition, similarity) || error("Type mismatch between partition and similarity")
#  size(partition.abundances, 1) == length(traits.traits) ||
#   error("Type mismatch between partition and traits")
#  size(partition.abundances, 1) == length(energy.energy) ||
#   error("Type mismatch between partition and energy")
#  A = typeof(partition.abundances)
#  Ecosystem{A, Part, Sim, T, E}(partition, similarity, Nullable{A}(), traits, energy)
#end


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

type TraitRelationship
  matrix::Matrix{Real}
end

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

