# Species list type - all info on species
type SpeciesList{T <: AbstractTraits,
                 E <: AbstractRequirement, TR <: AbstractTree,
                 MO <: AbstractMovement} <: AbstractTypes
  similarity::Matrix{Float64}
  names::Vector{String}
  traits::T
  abun::Vector{Int64}
  energy::E
  phylo::TR
  movement::MO

  function (::Type{SpeciesList{T, E, TR, MO}}){T <: AbstractTraits,
                   E <: AbstractRequirement, TR <: AbstractTree,
                   MO <: AbstractMovement}(similarity::Matrix{Float64},
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
      new{T, E, TR, MO}(similarity, names, traits, abun, energy, phylo, movement)
  end
end

function SpeciesList{E <: AbstractRequirement,
    MO <: AbstractMovement}(numspecies::Int64,
    numtraits::Int64, abun_dist::Distribution, energy::E,
    movement::MO)

    names = map(x -> "$x", 1:numspecies)
    # Create tree
    tree = jcoal(names, 100)
    # Create traits and assign to tips
    trts = map(string, 1:numtraits)
    assign_traits!(tree, 0.5, trts)
    # Get traits from tree
    sp_trt = BasicTrait(get_traits(tree, names, true)[:,1])
    # Create similarity matrix (for now identity)
    similarity = eye(numspecies)
    # Draw random set of abundances from distribution
    abun = rand(abun_dist)
    # error out when abun dist and NumberSpecies are not the same (same for energy dist)
    length(abun)==numspecies || throw(DimensionMismatch("Abundance vector
                                          doesn't match number species"))
    length(energy.energy)==numspecies || throw(DimensionMismatch("Requirement vector
                                          doesn't match number species"))
    size(similarity)==(numspecies, numspecies) || throw(DimensionMismatch("
                                Similarity matrix doesn't match number species"))

  SpeciesList{typeof(sp_trt), typeof(energy),
              typeof(tree), typeof(movement)}(similarity, names,
                                              sp_trt, abun, energy,
                                              tree, movement)
end

function _counttypes(sl::SpeciesList)
    return size(sl.similarity, 1)
end

function _calcsimilarity(sl::SpeciesList)
    return sl.similarity
end

function _floattypes(::SpeciesList)
    return Set([Float64])
end

function _getnames(sl::SpeciesList)
    return sl.names
end

function _calcordinariness(sl::SpeciesList, a::AbstractArray)
    _calcsimilarity(sl) * a
end
