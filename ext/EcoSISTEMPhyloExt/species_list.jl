# SPDX-License-Identifier: LGPL-3.0-or-later

# ---------------------------------------------------------------------------
# PARKED - the size-based `SpeciesList` constructor
# ---------------------------------------------------------------------------
# **Here rather than in `src/Species.jl` because uncommenting it would need `Phylo`**: it
# evolves body mass as a continuous trait along a phylogeny, as the deprecated tree-drawing
# constructors in `deprecations.jl` evolve niche classes. It is *parked*, not dead - see `SizeDemand`
# in `src/Demand.jl`. Keeping the two halves of that decision together is the point.

# The size-based `SpeciesList` constructor - COMMENTED OUT with `SizeDemand`
# ---------------------------------------------------------------------------
# **Parked, not dead.** This was `SizeDemand`'s only caller and cannot work without it, so the
# two are commented out together. It evolves body size along the phylogeny and derives both the
# initial abundances and the demand from it.
#
# **The part worth not losing** is the allometry on the `density` line below: the mass-abundance
# (self-thinning) relation `size^pop_mass / km^2`. It is computed here, used once for abundance, and
# then thrown away - which is exactly the material a metabolic (or footprint) reading of
# `SizeDemand` would need. See the note above `SizeDemand` in `src/Demand.jl`.
#
# `continuous_evolve`/`assigntraits!`/`resettraits!` are NOT affected - all are public and
# used elsewhere; checked before commenting this out.
#
# """
#     SpeciesList(numspecies::Int64, numtraits::Int64, pop_mass::Float64,
#       mean::Float64, var::Float64, area::Unitful.Area, movement::MO,
#       params::P, native::Vector{Bool}, switch::Vector{Float64})
#
# Create a `SpeciesList` where body size is evolved as a continuous trait along
# the phylogeny via Brownian motion with mean `mean` and variance `var`.
# Abundances and resource demands are derived from body size and population
# mass `pop_mass` scaled to `area`. Trait switching rates along the tree are
# controlled by `switch`.
# """
# function SpeciesList(numspecies::Int64,
#                      numtraits::Int64,
#                      pop_mass::Float64,
#                      mean::Float64,
#                      var::Float64,
#                      area::Unitful.Area,
#                      movement::MO,
#                      params::P,
#                      native::Vector{Bool},
#                      switch::Vector{Float64}) where {MO <: AbstractVector{<:SpeciesDispersal},
#                                                      P <: AbstractParams}
#     names = map(x -> "$x", 1:numspecies)
#     # Create tree
#     tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(names))
#     # Create tolerance and assign to tips
#     traits = DataFrame(trait1 = collect(1:numtraits))
#     assigntraits!(tree, switch, traits)
#     # Get tolerance from tree
#     sp_trt = SimpleCategoricalTolerance(Array(gettraits(tree, true)[:, 1]),
#                                         penalty = 0.5)
#     # Evolve size as a trait along the tree
#     EcoSISTEM.resettraits!(tree)
#     resource = abs.(_nichemeans(continuous_evolve(mean, var, tree)))
#     demand = SizeDemand(resource, pop_mass, area)
#     # Calculate density from size and nichefit
#     density = exp.(log.(resource) * pop_mass) ./ km^2
#     # Multiply density by area to get final population sizes
#     abun = round.(Int64, density * area)
#     # Create similarity matrix (for now identity)
#     phy = PhyloBranches(tree)
#     # Draw random set of abundances from distribution
#     # error out when abunance and NumberSpecies are not the same (same for resource dist)
#     length(abun) == numspecies || throw(DimensionMismatch("Abundance vector
#                                             doesn't match number species"))
#     counttypes(demand) == numspecies || throw(DimensionMismatch("Demand vector
#                                             doesn't match number species"))
#     return SpeciesList{typeof(sp_trt),
#                        typeof(demand),
#                        typeof(movement),
#                        typeof(phy),
#                        typeof(params)}(names,
#                                        sp_trt,
#                                        abun,
#                                        demand,
#                                        phy,
#                                        movement,
#                                        params,
#                                        native)
# end
