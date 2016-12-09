using Diversity
using Diversity.AbstractPartition
using Diversity.AbstractMetacommunity
using Diversity.psmatch

abstract AbstractHabitat

type Habitats <: AbstractHabitat
  matrix::Matrix{Float64}
end

abstract AbstractStructuredPartition{A, H <: AbstractHabitat} <:
 AbstractPartition{Float64, A}

type MatrixLandscape{A, H} <: AbstractStructuredPartition{A, H}
  abundances::A
  habitat::H
end

function MatrixLandscape{A, H}(abundances::A, habitats::H)
  dima = size(abundances)
  dimh = size(habitats.matrix)
  length(dima) == 3 || error("Abundances must be three dimensional!")
  ((dima[2] == dimh[1]) && (dima[3] == dimh[2])) ||
   error("Dimension mismatch between abundances and habitats")
  relative = abundances / sum(abundances)
  MatrixLandscape{typeof(relative), H}(relative, habitats)
end

abstract AbstractTraits

type StringTraits <: AbstractTraits
  traits::Vector{String}
end


abstract AbstractEcosystem{A, Part <: AbstractStructuredPartition, Sim, T <: AbstractTraits} <:
 AbstractMetacommunity{Float64, A, Part, Sim}
type Ecosystem{A, Part, Sim, T} <: AbstractEcosystem{A, Part, Sim, T}
  partition::Part
  similarity::Sim
  ordinariness::Nullable{A}
  traits::T
end

function Ecosystem{Part, Sim, T}(partition::Part, similarity::Sim, traits::T)
  psmatch(partition, similarity) || error("Type mismatch between partition and similarity")
  size(partition.abundances, 1) == length(traits.traits) ||
   error("Type mismatch between partition and traits")
  A = typeof(partition.abundances)
  Ecosystem{A, Part, Sim, T}(partition, similarity, Nullable{A}(),traits)
end
