"""
    TraitRelationship

The relationship between a trait and its environment, represented as a Matrix
of Floats.

"""
mutable struct TraitRelationship
  matrix::Matrix{Float64}
end
