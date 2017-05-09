using Diversity
importall Diversity.API


"""
    AbstractBudget

Abstract supertype for all budget types
"""
abstract AbstractBudget

function countsubcommunities(ab::AbstractBudget)
  return _countsubcommunities(ab)
end

"""
    Budget <: AbstractBudget

This budget type has a matrix of floats, representing the energy budget of each
subcommunity in the abiotic environment.
"""
type Budget <: AbstractBudget
  matrix::Matrix{Float64}
end

function _countsubcommunities(bud::Budget)
  return length(bud.matrix)
end

"""
    AbstractAbiotic{H <: AbstractHabitat, B <:AbstractBudget} <: AbstractPartition

Abstract supertype for all abiotic environment types and a subtype of
AbstractPartition
"""
abstract AbstractAbiotic{H <: AbstractHabitat, B <:AbstractBudget} <: AbstractPartition
"""
    GridAbioticEnv{H, B} <: AbstractAbiotic{H, B}

This abiotic environment type holds a habitat and budget, as well as a string
of subcommunity names.
"""
type GridAbioticEnv{H, B} <: AbstractAbiotic{H, B}
  habitat::H
  budget::B
  names::Vector{String}
  function (::Type{GridAbioticEnv{H, B}}){H, B}(habitat::H, budget::B,
                                names::Vector{String} =
                                map(x -> "$x", 1:countsubcommunities(habitat)))

    countsubcommunities(habitat) == countsubcommunities(budget) ||
      error("Habitat and budget must have same dimensions")
    countsubcommunities(habitat) == length(names) ||
      error("Number of subcommunities must match subcommunity names")
    return new(habitat, budget, names)
  end
end
"""
    simplenicheAE(numniches::Int64, dimension::Tuple,
                        maxBud::Float64, gridsquaresize::Float64)

Function to create a simple `Niches`, `Budget` type abiotic environment. Given a
number of niche types `numniches`, it creates a `Niches` environment with
dimensions `dimension` and a grid squaresize `gridsquaresize`. It also creates a
`Budget` type filled with the maximum budget value `maxbud`.
"""
function simplenicheAE(numniches::Int64, dimension::Tuple,
                        maxbud::Float64, gridsquaresize::Float64)
  # Create niches
  niches = map(string, 1:numniches)
  # Create niche-like environment
  hab = randomniches(dimension, niches, 0.5, repmat([0.5], numniches),
  gridsquaresize)
  # Create empty budget and for now fill with one value
  bud = zeros(dimension)
  fill!(bud, maxbud)

  return GridAbioticEnv(hab, Budget(bud))
end

function _countsubcommunities(gae::GridAbioticEnv)
  return countsubcommunities(gae.habitat)
end

function _getnames(gae::GridAbioticEnv)
    return gae.names
end
