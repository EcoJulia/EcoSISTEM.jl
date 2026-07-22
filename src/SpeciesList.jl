# SPDX-License-Identifier: LGPL-3.0-or-later

using Diversity
using Phylo

"""
    SpeciesList{TR <: AbstractTolerance, R <: AbstractDemand,
                MO <: AbstractMovement, T <: AbstractTypes,
                P <: AbstractParams} <: AbstractTypes

Species list housing all species-specific information. `names` holds species
names, `tolerance` encodes niche preferences, `abun` holds current abundances,
`demand` encodes resource needs, `types` holds the similarity structure,
`movement` describes dispersal, `params` holds demographic parameters, `native`
flags whether each species is native, and `susceptible` holds optional disease
susceptibility values.
"""
mutable struct SpeciesList{TR <: AbstractTolerance,
                           R <: AbstractDemand,
                           MO <: AbstractMovement,
                           T <: AbstractTypes,
                           P <: AbstractParams} <: AbstractTypes
    names::Vector{String}
    tolerance::TR
    abun::Vector{Int64}
    demand::R
    types::T
    movement::MO
    params::P
    native::Vector{Bool}
    susceptible::Vector{Union{Missing, Float64}}

    function SpeciesList{TR, R, MO, T, P}(names::Vector{String},
                                          tolerance::TR,
                                          abun::Vector{Int64},
                                          demand::R,
                                          types::T,
                                          movement::MO,
                                          params::P,
                                          native::Vector{Bool}) where {TR <:
                                                                       AbstractTolerance,
                                                                       R <:
                                                                       AbstractDemand,
                                                                       MO <:
                                                                       AbstractMovement,
                                                                       T <:
                                                                       AbstractTypes,
                                                                       P <:
                                                                       AbstractParams}
        # Check dimensions
        equal_param = equalpop(params, length(names))
        sus = Vector{Union{Missing, Float64}}(undef, length(names))
        return new{TR, R, MO, T, typeof(equal_param)}(names,
                                                      tolerance,
                                                      abun,
                                                      demand,
                                                      types,
                                                      movement,
                                                      equal_param,
                                                      native,
                                                      sus)
    end
    function SpeciesList{TR, R, MO, T, P}(tolerance::TR,
                                          abun::Vector{Int64},
                                          demand::R,
                                          types::T,
                                          movement::MO,
                                          params::P,
                                          native::Vector{Bool}) where {TR <:
                                                                       AbstractTolerance,
                                                                       R <:
                                                                       AbstractDemand,
                                                                       MO <:
                                                                       AbstractMovement,
                                                                       T <:
                                                                       AbstractTypes,
                                                                       P <:
                                                                       AbstractParams}
        # Assign names
        names = map(x -> "$x", eachindex(abun))
        equal_param = equalpop(params, length(names))
        sus = Vector{Union{Missing, Float64}}(undef, length(names))
        return new{TR, R, MO, T, typeof(equal_param)}(names,
                                                      tolerance,
                                                      abun,
                                                      demand,
                                                      types,
                                                      movement,
                                                      equal_param,
                                                      native,
                                                      sus)
    end
end

"""
    SpeciesList(numspecies::Int64, numtraits::Int64, abun::Vector{Int64},
      demand::R, movement::MO, params::P, native::Vector{Bool},
      switch::Vector{Float64})

Create a `SpeciesList` for `numspecies` species with `numtraits` discrete niche
tolerance evolved along a random ultrametric phylogeny. `switch` controls the rate
of trait change along branches. A `PhyloBranches` similarity structure is
computed from the tree. Abundances are provided via `abun` and resource
demands via `demand`.
"""
function SpeciesList(numspecies::Int64,
                     numtraits::Int64,
                     abun::Vector{Int64},
                     demand::R,
                     movement::MO,
                     params::P,
                     native::Vector{Bool},
                     switch::Vector{Float64}) where {R <: AbstractDemand,
                                                     MO <: AbstractMovement,
                                                     P <: AbstractParams}
    names = map(x -> "$x", 1:numspecies)
    # Create tree
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(names))
    # Create tolerance and assign to tips
    traits = DataFrame(trait1 = collect(1:numtraits))
    assign_traits!(tree, switch, traits)
    # Get tolerance from tree
    sp_trt = DiscreteTolerance(Array(get_traits(tree, true)[:, 1]))
    # Create similarity matrix (for now identity)
    phy = PhyloBranches(tree)
    # Draw random set of abundances from distribution
    if length(abun) < numspecies
        abun = vcat(abun, repmat([0], numspecies - length(abun)))
    end
    # error out when abun dist and NumberSpecies are not the same (same for resource dist)
    length(abun) == numspecies || throw(DimensionMismatch("Abundance vector
                                            doesn't match number species"))
    length(demand) == numspecies || throw(DimensionMismatch("Demand vector
                                            doesn't match number species"))
    return SpeciesList{typeof(sp_trt),
                       typeof(demand),
                       typeof(movement),
                       typeof(phy),
                       typeof(params)}(names,
                                       sp_trt,
                                       abun,
                                       demand,
                                       phy,
                                       movement,
                                       params,
                                       native)
end
function SpeciesList(numspecies::Int64,
                     numtraits::Int64,
                     abun::Vector{Int64},
                     demand::R,
                     movement::MO,
                     params::P,
                     native::Vector{Bool}) where {R <: AbstractDemand,
                                                  MO <: AbstractMovement,
                                                  P <: AbstractParams}
    return SpeciesList(numspecies, numtraits, abun, demand, movement, params,
                       native, [0.5])
end
@doc (@doc SpeciesList) SpeciesList(::Int64,
                                    ::Int64,
                                    ::Vector{Int64},
                                    ::R,
                                    ::MO,
                                    ::P,
                                    ::Vector{Bool}) where {R <:
                                                           AbstractDemand,
                                                           MO <:
                                                           AbstractMovement,
                                                           P <: AbstractParams}

"""
    SpeciesList(numspecies::Int64, numtraits::Int64, pop_mass::Float64,
      mean::Float64, var::Float64, area::Unitful.Area, movement::MO,
      params::P, native::Vector{Bool}, switch::Vector{Float64})

Create a `SpeciesList` where body size is evolved as a continuous trait along
the phylogeny via Brownian motion with mean `mean` and variance `var`.
Abundances and resource demands are derived from body size and population
mass `pop_mass` scaled to `area`. Trait switching rates along the tree are
controlled by `switch`.
"""
function SpeciesList(numspecies::Int64,
                     numtraits::Int64,
                     pop_mass::Float64,
                     mean::Float64,
                     var::Float64,
                     area::Unitful.Area,
                     movement::MO,
                     params::P,
                     native::Vector{Bool},
                     switch::Vector{Float64}) where {MO <: AbstractMovement,
                                                     P <: AbstractParams}
    names = map(x -> "$x", 1:numspecies)
    # Create tree
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(names))
    # Create tolerance and assign to tips
    traits = DataFrame(trait1 = collect(1:numtraits))
    assign_traits!(tree, switch, traits)
    # Get tolerance from tree
    sp_trt = DiscreteTolerance(Array(get_traits(tree, true)[:, 1]))
    # Evolve size as a trait along the tree
    EcoSISTEM.resettraits!(tree)
    resource = abs.(_nichemeans(ContinuousEvolve(mean, var, tree)))
    demand = SizeDemand(resource, pop_mass, area)
    # Calculate density from size and nichefit
    density = exp.(log.(resource) * pop_mass) ./ km^2
    # Multiply density by area to get final population sizes
    abun = round.(Int64, density * area)
    # Create similarity matrix (for now identity)
    phy = PhyloBranches(tree)
    # Draw random set of abundances from distribution
    # error out when abunance and NumberSpecies are not the same (same for resource dist)
    length(abun) == numspecies || throw(DimensionMismatch("Abundance vector
                                            doesn't match number species"))
    length(demand) == numspecies || throw(DimensionMismatch("Demand vector
                                            doesn't match number species"))
    return SpeciesList{typeof(sp_trt),
                       typeof(demand),
                       typeof(movement),
                       typeof(phy),
                       typeof(params)}(names,
                                       sp_trt,
                                       abun,
                                       demand,
                                       phy,
                                       movement,
                                       params,
                                       native)
end

"""
    SpeciesList(numspecies::Int64, numtraits::Int64, abun::Vector{Int64},
      demand::R, movement::MO, phy::T, params::P, native::Vector{Bool})

Create a `SpeciesList` with an explicitly supplied similarity structure `phy` of
type `AbstractTypes`, rather than computing a `PhyloBranches` similarity
internally. Discrete tolerance are still evolved along a random ultrametric tree.
"""
function SpeciesList(numspecies::Int64,
                     numtraits::Int64,
                     abun::Vector{Int64},
                     demand::R,
                     movement::MO,
                     phy::T,
                     params::P,
                     native::Vector{Bool}) where {R <: AbstractDemand,
                                                  MO <: AbstractMovement,
                                                  T <: AbstractTypes,
                                                  P <: AbstractParams}
    names = map(x -> "$x", 1:numspecies)
    # Create tree
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(names))
    # Create tolerance and assign to tips
    traits = DataFrame(trait1 = collect(1:numtraits))
    assign_traits!(tree, 0.5, traits)
    # Get tolerance from tree
    sp_trt = DiscreteTolerance(Array(get_traits(tree, true)[:, 1]))
    # Draw random set of abundances from distribution
    if length(abun) < numspecies
        abun = vcat(abun, repmat([0], numspecies - length(abun)))
    end
    # error out when abun dist and NumberSpecies are not the same (same for resource dist)
    length(abun) == numspecies || throw(DimensionMismatch("Abundance vector
                                            doesn't match number species"))
    length(demand) == numspecies || throw(DimensionMismatch("Demand vector
                                            doesn't match number species"))
    return SpeciesList{typeof(sp_trt),
                       typeof(demand),
                       typeof(movement),
                       typeof(phy),
                       typeof(params)}(names,
                                       sp_trt,
                                       abun,
                                       demand,
                                       phy,
                                       movement,
                                       params,
                                       native)
end

"""
    SpeciesList(numspecies::Int64, tolerance::TR, abun::Vector{Int64}, demand::R,
      movement::MO, params::P, native::Vector{Bool})

Create a `SpeciesList` from an explicitly supplied trait object `tolerance` of type
[`AbstractTolerance`](@ref). Uses `UniqueTypes` as the similarity structure,
treating all species as maximally distinct.
"""
function SpeciesList(numspecies::Int64,
                     tolerance::TR,
                     abun::Vector{Int64},
                     demand::R,
                     movement::MO,
                     params::P,
                     native::Vector{Bool}) where {TR <: AbstractTolerance,
                                                  R <: AbstractDemand,
                                                  MO <: AbstractMovement,
                                                  P <: AbstractParams}
    names = map(x -> "$x", 1:numspecies)
    # Create similarity matrix (for now identity)
    ty = UniqueTypes(numspecies)
    # Draw random set of abundances from distribution
    if length(abun) < numspecies
        abun = vcat(abun, fill(0, numspecies - length(abun)))
    end
    # error out when abun dist and NumberSpecies are not the same (same for resource dist)
    length(abun) == numspecies || throw(DimensionMismatch("Abundance vector
                                            doesn't match number species"))
    length(demand) == numspecies || throw(DimensionMismatch("Demand vector
                                            doesn't match number species"))
    return SpeciesList{typeof(tolerance),
                       typeof(demand),
                       typeof(movement),
                       typeof(ty),
                       typeof(params)}(names,
                                       tolerance,
                                       abun,
                                       demand,
                                       ty,
                                       movement,
                                       params,
                                       native)
end

"""
    getdemand(sppl::SpeciesList)

Return the total resource usage for the species list `sppl`, combining abundances
with per-species resource demands.
"""
function getdemand(sppl::SpeciesList)
    return _getdemand(sppl.abun, sppl.demand)
end

function _simmatch(sim::SpeciesList)
    return _simmatch(sim.types)
end

#function _calcsimilarity(ut::UniqueTypes)
#  return eye(ut.num)
#end

#function _calcsimilarity(ph::PhyloBranches)
# return ph.Zmatrix
#end

import Diversity.API: _gettypenames
function _gettypenames(sl::SpeciesList, input::Bool)
    return _gettypenames(sl.types, input)
end

import Diversity.API: _counttypes
function _counttypes(sl::SpeciesList, input::Bool)
    return _counttypes(sl.types, input)
end

import Diversity.API: _calcsimilarity
function _calcsimilarity(sl::SpeciesList, a::AbstractArray)
    return _calcsimilarity(sl.types, a)
end

import Diversity.API: floattypes
function floattypes(::SpeciesList)
    return Set([Float64])
end

import Diversity.API: _calcordinariness
function _calcordinariness(sl::SpeciesList, a::AbstractArray)
    return _calcordinariness(sl.types, a, one(eltype(a)))
end

import Diversity.API: _calcabundance
function _calcabundance(sl::SpeciesList, a::AbstractArray)
    return _calcabundance(sl.types, a)
end

import Diversity.API._getdiversityname
function _getdiversityname(sl::SpeciesList)
    return _getdiversityname(sl.types)
end
