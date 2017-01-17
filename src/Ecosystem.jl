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

abstract AbstractStructuredPartition{A, H <: AbstractHabitat, B<: AbstractBudget} <:
 AbstractPartition{Float64, A}

type MatrixLandscape{A, H, B} <: AbstractStructuredPartition{A, H, B}
  abundances::A
  habitat::H
  budget::B
end

function MatrixLandscape{A, H, B}(abundances::A, habitats::H, budget::B)
  dima = size(abundances)
  dimh = size(habitats.matrix)
  dimb = size(budget.matrix)
  length(dima) == 3 || error("Abundances must be three dimensional!")
  ((dima[2] == dimh[1]) && (dima[3] == dimh[2])) ||
   error("Dimension mismatch between abundances and habitats")
  ((dimh[1] == dimb[1]) && (dimh[2] == dimb[2])) ||
   error("Dimension mismatch between budgets and habitats")
  relative = abundances #/ sum(abundances)
  MatrixLandscape{typeof(relative), H, B}(relative, habitats, budget)
end

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


abstract AbstractEcosystem{A, Part <: AbstractStructuredPartition, Sim, T <: AbstractTraits, E <: AbstractEnergy} <:
 AbstractMetacommunity{Float64, A, Part, Sim}

type Ecosystem{A, Part, Sim, T, E} <: AbstractEcosystem{A, Part, Sim, T, E}
  partition::Part
  similarity::Sim
  ordinariness::Nullable{A}
  traits::T
  energy::E
end

function Ecosystem{Part, Sim, T, E}(partition::Part, similarity::Sim, traits::T, energy::E)
  psmatch(partition, similarity) || error("Type mismatch between partition and similarity")
  size(partition.abundances, 1) == length(traits.traits) ||
   error("Type mismatch between partition and traits")
  size(partition.abundances, 1) == length(energy.energy) ||
   error("Type mismatch between partition and energy")
  A = typeof(partition.abundances)
  Ecosystem{A, Part, Sim, T, E}(partition, similarity, Nullable{A}(), traits, energy)
end


abstract AbstractSpeciesList{A, N, B, T <: AbstractTraits, Sim, E<: AbstractEnergy}

type SpeciesList{A, N, B, T, Sim, E} <: AbstractSpeciesList{A, N, B, T, Sim, E}
  traits::T
  similarity::Sim
  abun::Vector{A}
  energy::E
  phylo::Tree{N,B}
end
