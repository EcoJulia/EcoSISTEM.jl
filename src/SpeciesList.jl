using Diversity.Phylogenetics
import Diversity.getnames
"""
    simmatch(sim::AbstractTypes)

Checks for size and value incompatibility in similarity objects
"""
function simmatch end

function _simmatch(sim::Phylogeny)
  # Check similarity is square matrix
  size(sim.Zmatrix, 1) == size(sim.Zmatrix, 2) ||
    throw(DimensionMismatch("Similarity matrix is not square"))

  # Check similarity is bounded between 0 and 1
  minimum(sim.Zmatrix) ≥ 0 || throw(DomainError())
  maximum(sim.Zmatrix) ≤ 1 || warn("Similarity matrix has values above 1")
end

function _simmatch(sim::UniqueTypes)
  length(getnames(sim)) == counttypes(sim) ||
    throw(DimensionMismatch("Names do not match number of types"))
end
"""
    SpeciesList{T <: AbstractTraits, R <: AbstractRequirement,
                 MO <: AbstractMovement} <: AbstractTypes
Species list houses all species-specific information including trait information,
phylogenetic relationships, requirement for energy and movement types.
"""
type SpeciesList{T <: AbstractTraits,
                 R <: AbstractRequirement,
                 MO <: AbstractMovement,
                 P <: AbstractTypes} <: AbstractTypes
  names::Vector{String}
  traits::T
  abun::Vector{Int64}
  requirement::R
  phylo::P
  movement::MO

  function (::Type{SpeciesList{T, R, MO, P}}){T <: AbstractTraits,
                   R <: AbstractRequirement,
                   MO <: AbstractMovement,
                   P <: AbstractTypes}(names:: Vector{String},
      traits::T, abun::Vector{Int64}, req::R,
      phylo::P, movement::MO)
      # Check dimensions
      _simmatch(phylo)
      new{T, R, MO, P}(names, traits, abun, req, phylo, movement)
  end
  function (::Type{SpeciesList{T, R, MO, P}}){T <: AbstractTraits,
                   R <: AbstractRequirement,
                   MO <: AbstractMovement,
                   P <: AbstractTypes}(
      traits::T, abun::Vector{Int64}, req::R,
      phylo::P, movement::MO)
      # Check dimensions
      _simmatch(phylo)
      # Assign names
      names = map(x -> "$x", 1:length(abun))
      new{T, R, MO, P}(names, traits, abun, req, phylo, movement)
  end
end
"""
    SpeciesList{R <: AbstractRequirement,
      MO <: AbstractMovement}(numspecies::Int64,
      numtraits::Int64, abun_dist::Distribution, req::R,
      movement::MO)
Function to create a SpeciesList given a number of species, the number of traits
they possess, their abundances, requirement from the environment and their
movement kernel.
"""
function SpeciesList{R <: AbstractRequirement,
    MO <: AbstractMovement}(numspecies::Int64,
    numtraits::Int64, abun_dist::Distribution, req::R,
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
    phy = Phylogeny(tree)
    # Draw random set of abundances from distribution
    abun = rand(abun_dist)
    # error out when abun dist and NumberSpecies are not the same (same for energy dist)
    length(abun)==numspecies || throw(DimensionMismatch("Abundance vector
                                          doesn't match number species"))
    length(req.energy)==numspecies || throw(DimensionMismatch("Requirement vector
                                          doesn't match number species"))
  SpeciesList{typeof(sp_trt), typeof(req),
              typeof(movement), typeof(phy)}(names, sp_trt, abun, req,
                                              phy, movement)
end
"""
    SpeciesList{R <: AbstractRequirement,
      MO <: AbstractMovement, P <: AbstractTypes}(numspecies::Int64,
      numtraits::Int64, abun_dist::Distribution, req::R,
      movement::MO, phy::P)
Function to create a SpeciesList given a number of species, the number of traits
they possess, their abundances, requirement from the environment and their
movement kernel and any type of AbstractTypes.
"""
function SpeciesList{R <: AbstractRequirement,
    MO <: AbstractMovement, P <: AbstractTypes}(numspecies::Int64,
    numtraits::Int64, abun_dist::Distribution, req::R,
    movement::MO, phy::P)

    names = map(x -> "$x", 1:numspecies)
    # Create tree
    tree = jcoal(names, 100)
    # Create traits and assign to tips
    trts = map(string, 1:numtraits)
    assign_traits!(tree, 0.5, trts)
    # Get traits from tree
    sp_trt = BasicTrait(get_traits(tree, names, true)[:,1])
    # Draw random set of abundances from distribution
    abun = rand(abun_dist)
    # error out when abun dist and NumberSpecies are not the same (same for energy dist)
    length(abun)==numspecies || throw(DimensionMismatch("Abundance vector
                                          doesn't match number species"))
    length(req.energy)==numspecies || throw(DimensionMismatch("Requirement vector
                                          doesn't match number species"))
  SpeciesList{typeof(sp_trt), typeof(req),
              typeof(movement), typeof(phy)}(names, sp_trt, abun, req,
                                              phy, movement)
end

function _simmatch(sim::SpeciesList)
  _simmatch(sim.phylo)
end

function _calcsimilarity(ut::UniqueTypes)
  return eye(ut.num)
end

function _calcsimilarity(ph::Phylogeny)
return ph.Zmatrix
end

function _getnames(sl::SpeciesList)
    return getnames(sl.phylo)
end

function _counttypes(sl::SpeciesList)
    return length(getnames(sl))
end

function _calcsimilarity(sl::SpeciesList)
    return _calcsimilarity(sl.phylo)
end

function _floattypes(::SpeciesList)
    return Set([Float64])
end

function _calcordinariness(sl::SpeciesList, a::AbstractArray)
    _calcsimilarity(sl) * a
end
