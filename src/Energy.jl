"""
    AbstractRequirement{Energy}

Abstract supertype for all species energy requirement types, parameterised by
the type(s) of energy required `Energy`.
"""
abstract AbstractRequirement{Energy}
"""
    SimpleRequirement <: AbstractRequirement{Float64}

A simple energy requirement is a single float for each species.
"""
type SimpleRequirement <: AbstractRequirement{Float64}
  energy::Vector{Float64}
end

"""
    AbstractBudget

Abstract supertype for all budget types
"""
abstract AbstractBudget{Requirement}

function countsubcommunities(ab::AbstractBudget)
  return _countsubcommunities(ab)
end

"""
    SimpleBudget <: AbstractBudget{Float64}

This budget type has a matrix of floats, representing the energy budget of each
subcommunity in the abiotic environment.
"""
type SimpleBudget <: AbstractBudget{Float64}
  matrix::Matrix{Float64}
end

function _countsubcommunities(bud::SimpleBudget)
  return length(bud.matrix)
end
