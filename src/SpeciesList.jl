using Diversity.Phylogenetics
importall Diversity.API

"""
    simmatch(sim::AbstractTypes)

Checks for size and value incompatibility in similarity objects
"""
function simmatch end

function _simmatch(sim::PhyloTypes)
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
    SpeciesList{TR <: AbstractTraits, R <: AbstractRequirement,
                MO <: AbstractMovement, T <: AbstractTypes,
                P <: AbstractParams} <: AbstractTypes
Species list houses all species-specific information including trait information,
phylogenetic relationships, requirement for energy and movement types.
"""
mutable struct SpeciesList{TR <: AbstractTraits,
                 R <: AbstractRequirement,
                 MO <: AbstractMovement,
                 T <: AbstractTypes,
                 P <: AbstractParams} <: AbstractTypes
  names::Vector{String}
  traits::TR
  abun::Vector{Int64}
  requirement::R
  types::T
  movement::MO
  params::P

  function SpeciesList{TR, R, MO, T, P}(names:: Vector{String},
      traits::TR, abun::Vector{Int64}, req::R,
      types::T, movement::MO, params::P) where {TR <: AbstractTraits,
                       R <: AbstractRequirement,
                       MO <: AbstractMovement,
                       T <: AbstractTypes,
                       P <: AbstractParams}
      # Check dimensions
      _simmatch(types)
      new{TR, R, MO, T, P}(names, traits, abun, req, types, movement, params)
  end
  function SpeciesList{TR, R, MO, T, P}(
      traits::TR, abun::Vector{Int64}, req::R,
      types::T, movement::MO, params::P) where {TR <: AbstractTraits,
                       R <: AbstractRequirement,
                       MO <: AbstractMovement,
                       T <: AbstractTypes,
                       P <: AbstractParams}
      # Check dimensions
      _simmatch(types)
      # Assign names
      names = map(x -> "$x", 1:length(abun))
      new{TR, R, MO, T, P}(names, traits, abun, req, types, movement, params)
  end
end
"""
    SpeciesList{R <: AbstractRequirement,
      MO <: AbstractMovement, P <: AbstractParams}(numspecies::Int64,
      numtraits::Int64, abun_dist::Distribution, req::R,
      movement::MO, params::P)
Function to create a SpeciesList given a number of species, the number of traits
they possess, their abundances, requirement from the environment and their
movement kernel.
"""
function SpeciesList{R <: AbstractRequirement,
    MO <: AbstractMovement, P <: AbstractParams}(numspecies::Int64,
    numtraits::Int64, abun_dist::Distribution, req::R,
    movement::MO, params::P)

    names = map(x -> "$x", 1:numspecies)
    # Create tree
    tree = jcoal(names, 100.0)
    # Create traits and assign to tips
    trts = map(string, 1:numtraits)
    assign_traits!(tree, 0.5, trts)
    # Get traits from tree
    sp_trt = BasicTrait(get_traits(tree, names, true)[:,1])
    # Create similarity matrix (for now identity)
    phy = PhyloTypes(tree)
    # Draw random set of abundances from distribution
    abun = rand(abun_dist)
    # error out when abun dist and NumberSpecies are not the same (same for energy dist)
    length(abun)==numspecies || throw(DimensionMismatch("Abundance vector
                                          doesn't match number species"))
    length(req.energy)==numspecies || throw(DimensionMismatch("Requirement vector
                                          doesn't match number species"))
  SpeciesList{typeof(sp_trt), typeof(req),
              typeof(movement), typeof(phy),typeof(params)}(names, sp_trt, abun,
              req, phy, movement, params)
end

function SpeciesList{R <: AbstractRequirement,
    MO <: AbstractMovement}(numspecies::Int64,
    numtraits::Int64, abun_dist::Distribution, req::R,
    movement::MO, params::EqualPop)

    equal_params = equalpop(params, numspecies)
    return SpeciesList(numspecies, numtraits, abun_dist, req, movement, equal_params)
end
"""
    SpeciesList{R <: AbstractRequirement, MO <: AbstractMovement,
      T <: AbstractTypes, P <: AbstractParams}(numspecies::Int64,
      numtraits::Int64, abun_dist::Distribution, req::R,
      movement::MO, phy::T, params::P)
Function to create a SpeciesList given a number of species, the number of traits
they possess, their abundances, requirement from the environment and their
movement kernel and any type of AbstractTypes.
"""
function SpeciesList{R <: AbstractRequirement, MO <: AbstractMovement,
    T <: AbstractTypes, P <: AbstractParams}(numspecies::Int64,
    numtraits::Int64, abun_dist::Distribution, req::R,
    movement::MO, phy::T, params::P)

    names = map(x -> "$x", 1:numspecies)
    # Create tree
    tree = jcoal(names, 100.0)
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
              typeof(movement), typeof(phy), typeof(params)}(names, sp_trt, abun,
               req, phy, movement, params)
end

function SpeciesList{R <: AbstractRequirement, MO <: AbstractMovement,
    T <: AbstractTypes}(numspecies::Int64,
    numtraits::Int64, abun_dist::Distribution, req::R,
    movement::MO, phy::T, params::EqualPop)

    equal_params = equalpop(params, numspecies)
    return SpeciesList(numspecies, numtraits, abun_dist, req, movement, phy, equal_params)
end

function _simmatch(sim::SpeciesList)
  _simmatch(sim.types)
end

#function _calcsimilarity(ut::UniqueTypes)
#  return eye(ut.num)
#end

#function _calcsimilarity(ph::PhyloTypes)
# return ph.Zmatrix
#end

function _getnames(sl::SpeciesList, input::Bool)
    return _getnames(sl.types, input)
end

function _counttypes(sl::SpeciesList, input::Bool)
    return _counttypes(sl.types, input)
end

function _calcsimilarity(sl::SpeciesList, a::AbstractArray)
    return _calcsimilarity(sl.types, a)
end

function _floattypes(::SpeciesList)
    return Set([Float64])
end

function _calcordinariness(sl::SpeciesList, a::AbstractArray)
    _calcordinariness(sl.types, a)
end

function _calcabundance(sl::SpeciesList, a::AbstractArray)
  return _calcabundance(sl.types, a)
end
