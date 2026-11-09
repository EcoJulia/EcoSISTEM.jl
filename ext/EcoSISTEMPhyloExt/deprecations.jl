# SPDX-License-Identifier: LGPL-3.0-or-later

# --- Deprecations whose bodies need Phylo ---------------------------------------------------------
#
# A deprecated method that draws a tree cannot be written without `Phylo` loaded, so these are here
# rather than in `src/deprecations.jl`, labelled and ordered as that file is. Their docstring is
# there, on the form taking no `switch`, since no docstring can live in an extension.

# ---------------------------------------------------------------------------
# A species list with a phylogeny is built from a pool given the tree (`add_phylogeny!`) and the
# traits evolved along it (`add_trait!`). These draw the tree, and the niche classes evolved along it,
# as they did, from the same stream, so a seed gives the same species; the first's types follow the
# species' order, as a pool's do, so each species' abundances meet its own tip. The second keeps the
# tolerances' tree apart from the similarity it is given, as it always did.
#
# Deprecated in v0.11.0.
# ---------------------------------------------------------------------------
function EcoSISTEM.SpeciesList(numspecies::Int64,
                               numtraits::Int64,
                               demand::DM,
                               movement::MO,
                               params::P,
                               switch::Vector{Float64};
                               seed = nothing) where {DM <: AbstractDemand,
                                                      MO <:
                                                      AbstractVector{<:SpeciesDispersal},
                                                      P <: AbstractParams}
    Base.depwarn(_treelistwarning(), :SpeciesList)
    names = map(x -> "$x", 1:numspecies)
    movement = _dispersalvector(movement)
    rng = _phylogenyrng(seed)
    # Create tree
    tree = rand(rng,
                Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(names))
    # Create tolerance and assign to tips
    traits = DataFrame(trait1 = collect(1:numtraits))
    assigntraits!(rng, tree, switch, traits)
    # Get tolerance from tree
    # `penalty = 0.5` is **soft** exclusion, and must be stated rather than left to the default:
    # this evolves a *niche* label, and a species outside its niche should do worse there rather than
    # be unable to live there at all. The tolerance defaults to `0.0` (hard), which is right for a
    # habitat filter and wrong for this.
    sp_trt = SimpleCategoricalTolerance(Array(gettraits(tree, true)[:, 1]),
                                        axis = EcoSISTEM.TypologyAxis,
                                        penalty = 0.5)
    # The types in the species' order, so that each species' abundances meet its own tip.
    phy = PhyloBranches(tree, names)
    counttypes(demand) == numspecies ||
        throw(DimensionMismatch("Demand vector
 doesn't match number species"))
    return SpeciesList{typeof(sp_trt),
                       typeof(demand),
                       typeof(movement),
                       typeof(phy),
                       typeof(params)}(names,
                                       sp_trt,
                                       demand,
                                       phy,
                                       movement,
                                       params)
end

function EcoSISTEM.SpeciesList(numspecies::Int64,
                               numtraits::Int64,
                               demand::DM,
                               movement::MO,
                               phy::T,
                               params::P;
                               seed = nothing) where {DM <: AbstractDemand,
                                                      MO <:
                                                      AbstractVector{<:SpeciesDispersal},
                                                      T <: AbstractTypes,
                                                      P <: AbstractParams}
    Base.depwarn(_treelistwarning(), :SpeciesList)
    names = map(x -> "$x", 1:numspecies)
    movement = _dispersalvector(movement)
    rng = _phylogenyrng(seed)
    # Create tree
    tree = rand(rng,
                Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(names))
    # Create tolerance and assign to tips
    traits = DataFrame(trait1 = collect(1:numtraits))
    assigntraits!(rng, tree, 0.5, traits)
    # Get tolerance from tree
    # `penalty = 0.5` is **soft** exclusion, and must be stated rather than left to the default:
    # this evolves a *niche* label, and a species outside its niche should do worse there rather than
    # be unable to live there at all. The tolerance defaults to `0.0` (hard), which is right for a
    # habitat filter and wrong for this.
    sp_trt = SimpleCategoricalTolerance(Array(gettraits(tree, true)[:, 1]),
                                        axis = EcoSISTEM.TypologyAxis,
                                        penalty = 0.5)
    counttypes(demand) == numspecies ||
        throw(DimensionMismatch("Demand vector
 doesn't match number species"))
    return SpeciesList{typeof(sp_trt),
                       typeof(demand),
                       typeof(movement),
                       typeof(phy),
                       typeof(params)}(names,
                                       sp_trt,
                                       demand,
                                       phy,
                                       movement,
                                       params)
end

# The warning both forms give, with the route that replaces them.
function _treelistwarning()
    return "`SpeciesList(numspecies, numtraits, ...)`, drawing a tree and evolving niche classes " *
           "along it, is deprecated: draw the tree with Phylo (`rand(Ultrametric(names))`), give " *
           "it to a pool with `add_phylogeny!`, the classes evolved along it with `add_trait!`, " *
           "and build the species from the pool with `build_species`."
end
