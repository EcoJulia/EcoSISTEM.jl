using Diversity
using Phylo
using Compat


abstract type AbstractSpeciesList <: AbstractTypes end

abstract type AbstractSpeciesTypes <: AbstractTypes end

abstract type AbstractPathogenTypes <: AbstractTypes end

mutable struct SpeciesList{S <: AbstractSpeciesTypes, P <: AbstractPathogenTypes, PA <: AbstractParams} <: AbstractSpeciesList
    species::S
    pathogens::P
    params::PA
end

"""
    SpeciesTypes{TR <: AbstractTraits, R <: AbstractRequirement,
                MO <: AbstractMovement, T <: AbstractTypes,
                P <: AbstractParams} <: AbstractTypes
Species list houses all species-specific information including trait information,
phylogenetic relationships, requirement for energy and movement types.
"""
mutable struct SpeciesTypes{TR <: AbstractTraits,
                 R <: AbstractRequirement,
                 MO <: AbstractMovement,
                 T <: AbstractTypes} <: AbstractSpeciesTypes
  names::Vector{String}
  traits::TR
  abun::Vector{Int64}
  requirement::R
  types::T
  movement::MO
  native::Vector{Bool}
  susceptible::Vector{Union{Missing, Float64}}

  function SpeciesTypes{TR, R, MO, T}(names:: Vector{String},
      traits::TR, abun::Vector{Int64}, req::R,
      types::T, movement::MO, native::Vector{Bool}) where {
                       TR <: AbstractTraits,
                       R <: AbstractRequirement,
                       MO <: AbstractMovement,
                       T <: AbstractTypes}
      # Check dimensions
      sus = Vector{Union{Missing, Float64}}(Compat.undef, length(names))
      new{TR, R, MO, T}(names, traits, abun, req, types,
       movement, native, sus)
  end
  function SpeciesTypes{TR, R, MO, T}(
      traits::TR, abun::Vector{Int64}, req::R,
      types::T, movement::MO, native::Vector{Bool}) where {
                       TR <: AbstractTraits,
                       R <: AbstractRequirement,
                       MO <: AbstractMovement,
                       T <: AbstractTypes}
      # Assign names
      names = map(x -> "$x", 1:length(abun))
      sus = Vector{Union{Missing, Float64}}(Compat.undef, length(names))
      new{TR, R, MO, T}(names, traits, abun, req, types,
       movement, native, sus)
  end
end

mutable struct NoPathogen <: AbstractPathogenTypes end
function _counttypes(virus::NoPathogen, input::Bool)
    return 0
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
function SpeciesList(numspecies::Int64,
    numtraits::Int64, abun::Vector{Int64}, req::R,
    movement::MO, params::PA, native::Vector{Bool}, switch::Vector{Float64},
    pathogens::P = NoPathogen()) where {R <: AbstractRequirement,
        MO <: AbstractMovement, PA <: AbstractParams, P <: AbstractPathogenTypes}

    names = map(x -> "$x", 1:numspecies)
    # Create tree
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(names))
    # Create traits and assign to tips
    trts = DataFrame(trait1 = collect(1:numtraits))
    assign_traits!(tree, switch, trts)
    # Get traits from tree
    sp_trt = DiscreteTrait(Array(get_traits(tree, true)[:, 1]))
    # Create similarity matrix (for now identity)
    phy = PhyloBranches(tree)
    # Draw random set of abundances from distribution
    if length(abun) < numspecies
        abun = vcat(abun, repmat([0], numspecies - length(abun)))
    end
    # error out when abun dist and NumberSpecies are not the same (same for energy dist)
    length(abun)==numspecies || throw(DimensionMismatch("Abundance vector
                                          doesn't match number species"))
    length(req)==numspecies || throw(DimensionMismatch("Requirement vector
                                          doesn't match number species"))

    params = equalpop(params, length(names))
    species = SpeciesTypes{typeof(sp_trt), typeof(req),
              typeof(movement), typeof(phy)}(names, sp_trt, abun,
              req, phy, movement, native)
    return SpeciesList{typeof(species), typeof(pathogens), typeof(params)}(species, pathogens, params)
end
function SpeciesList(numspecies::Int64,
    numtraits::Int64, abun::Vector{Int64}, req::R,
    movement::MO, params::P, native::Vector{Bool}) where {R <: AbstractRequirement,
        MO <: AbstractMovement, P <: AbstractParams}
        return SpeciesList(numspecies, numtraits, abun, req, movement,
                params, native, [0.5])
end

function SpeciesList(numspecies::Int64,
    numtraits::Int64, pop_mass::Float64, mean::Float64, var::Float64,
    area::Unitful.Area, movement::MO, params::PA, native::Vector{Bool},
    switch::Vector{Float64},
    pathogens::P = NoPathogen()) where {MO <: AbstractMovement, PA <: AbstractParams,
    P <: AbstractPathogenTypes}

    names = map(x -> "$x", 1:numspecies)
    # Create tree
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(names))
    # Create traits and assign to tips
    trts = DataFrame(trait1 = collect(1:numtraits))
    assign_traits!(tree, switch, trts)
    # Get traits from tree
    sp_trt = DiscreteTrait(Array(get_traits(tree, true)[:, 1]))
    # Evolve size as a trait along the tree
    EcoSISTEM.resettraits!(tree)
    energy = abs.(ContinuousEvolve(mean, var, tree).mean)
    req = SizeRequirement(energy, pop_mass, area)
    # Calculate density from size and relationship
    density = exp.(log.(energy) * pop_mass)./ km^2
    # Multiply density by area to get final population sizes
    abun = round.(Int64, density * area)
    # Create similarity matrix (for now identity)
    phy = PhyloBranches(tree)
    # Draw random set of abundances from distribution
    # error out when abunance and NumberSpecies are not the same (same for energy dist)
    length(abun)==numspecies || throw(DimensionMismatch("Abundance vector
                                          doesn't match number species"))
    length(req)==numspecies || throw(DimensionMismatch("Requirement vector
                                          doesn't match number species"))
    params = equalpop(params, length(names))
    species = SpeciesTypes{typeof(sp_trt), typeof(req),
              typeof(movement), typeof(phy)}(names, sp_trt, abun,
              req, phy, movement, native)
    return SpeciesList{typeof(species), typeof(pathogens), typeof(params)}(species, pathogens, params)
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
function SpeciesList(numspecies::Int64,
    numtraits::Int64, abun::Vector{Int64}, req::R,
    movement::MO, phy::T, params::PA, native::Vector{Bool},
    pathogens::P = NoPathogen()) where
    {R <: AbstractRequirement, MO <: AbstractMovement,
        T <: AbstractTypes, PA <: AbstractParams, P <: AbstractPathogenTypes}

    names = map(x -> "$x", 1:numspecies)
    # Create tree
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(names))
    # Create traits and assign to tips
    trts = DataFrame(trait1 = collect(1:numtraits))
    assign_traits!(tree, 0.5, trts)
    # Get traits from tree
    sp_trt = DiscreteTrait(Array(get_traits(tree, true)[:, 1]))
    # Draw random set of abundances from distribution
    if length(abun) < numspecies
        abun = vcat(abun, repmat([0], numspecies - length(abun)))
    end
    # error out when abun dist and NumberSpecies are not the same (same for energy dist)
    length(abun)==numspecies || throw(DimensionMismatch("Abundance vector
                                          doesn't match number species"))
    length(req)==numspecies || throw(DimensionMismatch("Requirement vector
                                          doesn't match number species"))
   params = equalpop(params, length(names))
   species = SpeciesTypes{typeof(sp_trt), typeof(req),
              typeof(movement), typeof(phy)}(names, sp_trt, abun,
               req, phy, movement, native)
    return SpeciesList{typeof(species), typeof(pathogens), typeof(params)}(species, pathogens, params)
end



function SpeciesList(numspecies::Int64,
    traits::TR, abun::Vector{Int64}, req::R,
    movement::MO, params::PA, native::Vector{Bool},
    pathogens::P = NoPathogen()) where
    {TR<: AbstractTraits, R <: AbstractRequirement,
        MO <: AbstractMovement, PA <: AbstractParams,
        P <: AbstractPathogenTypes}

    names = map(x -> "$x", 1:numspecies)
    # Create similarity matrix (for now identity)
    ty = UniqueTypes(numspecies)
    # Draw random set of abundances from distribution
    if length(abun) < numspecies
        abun = vcat(abun, fill(0, numspecies - length(abun)))
    end
    # error out when abun dist and NumberSpecies are not the same (same for energy dist)
    length(abun)==numspecies || throw(DimensionMismatch("Abundance vector
                                          doesn't match number species"))
    length(req)==numspecies || throw(DimensionMismatch("Requirement vector
                                          doesn't match number species"))
  params = equalpop(params, length(names))
  species = SpeciesTypes{typeof(traits), typeof(req),
              typeof(movement), typeof(ty)}(names, traits, abun,
              req, ty, movement, native)
  return SpeciesList{typeof(species), typeof(pathogens), typeof(params)}(species, pathogens, params)
end


function getenergyusage(sppl::SpeciesList)
    return _getenergyusage(sppl.species.abun, sppl.species.requirement)
end

function getnames(sppl::SpeciesList)
    return sppl.species.names
end

function _simmatch(sim::SpeciesList)
  _simmatch(sim.species.types)
end

#function _calcsimilarity(ut::UniqueTypes)
#  return eye(ut.num)
#end

#function _calcsimilarity(ph::PhyloTypes)
# return ph.Zmatrix
#end
import Diversity.API: _gettypenames
function _gettypenames(sl::SpeciesList, input::Bool)
    return _gettypenames(sl.species.types, input)
end
import Diversity.API: _counttypes
function _counttypes(sl::SpeciesList{A, B, C}, input::Bool) where {A <: SpeciesTypes, B <: AbstractPathogenTypes, C <: AbstractParams}
    return _counttypes(sl.species.types, input)
end

import Diversity.API: _calcsimilarity
function _calcsimilarity(sl::SpeciesList, a::AbstractArray)
    return _calcsimilarity(sl.species.types, a)
end
import Diversity.API: floattypes
function floattypes(::SpeciesList)
    return Set([Float64])
end
import Diversity.API: _calcordinariness
function _calcordinariness(sl::SpeciesList, a::AbstractArray)
    _calcordinariness(sl.species.types, a, one(eltype(a)))
end
import Diversity.API: _calcabundance
function _calcabundance(sl::SpeciesList, a::AbstractArray)
  return _calcabundance(sl.species.types, a)
end
import Diversity.API._getdiversityname
function _getdiversityname(sl::SpeciesList)
    return _getdiversityname(sl.species.types)
end
