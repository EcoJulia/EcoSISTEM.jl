# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Conformance to Diversity.jl, not our own API — every method this package defines on a `Diversity`,
# `Diversity.API` or `Diversity.Ecology` generic, ordered **by type**, so that "what does Diversity
# see when it is handed one of these?" can be answered by reading down. Its siblings are
# `BaseInterface.jl` and `EcoBaseInterface.jl`.
#
# Three rules hold throughout, and each has a failure mode that is silent:
#
#   * **Every definition is qualified.** Written unqualified in a file with no matching `import`,
#     `function _calcsimilarity(...)` defines `EcoSISTEM._calcsimilarity` instead of extending
#     Diversity's. The package loads, every gate passes, and the interface is broken.
#   * **`SpeciesList` is an `AbstractTypes` that HOLDS the real types**, in `sl.types`, so every
#     `AbstractTypes` hook has to delegate there. Diversity's defaults are not neutral —
#     `_hassimilarity` defaults to `true` — so a missing delegation answers plausibly and wrongly. A
#     test sweeps `Diversity.API` rather than listing the hooks, so a new one cannot be forgotten.
#   * **A call graph cannot see any of this.** A hook's only other mention in the repo is the
#     `import` line, so every one of them looks dead. Check `docs/overloads.md` before believing a
#     dead-code report.
#
# `EcoSISTEMMPIExt` and `EcoSISTEMPhyloExt` carry six more, which cannot live here: an extension is
# not loaded until its trigger package is.

using Diversity

# ══ Demand and SpeciesRequirementCollection ════════════════════════════════════════════════════════
# `length` on a collection means its **member** count (`src/collections.jl`), so the species count has
# to be asked for by name, and the name used is Diversity's rather than one of ours. A `Demand` is a
# per-species vector, so "how many types does this cover" is exactly `counttypes`.
#
# Diversity's own dispatch cannot reach a `Demand`, which is not an `EcoBase.AbstractThings`, so this
# is a method on their generic for our type — the same shape as `countsubcommunities` on a regime and
# a supply below, and `EcoBase.placenames` on a habitat.
#
# The `_counttypes` hook is deliberately not used: it takes `(t, raw::Bool)` everywhere in Diversity,
# and the raw-versus-processed distinction means nothing for a plain per-species vector, so nothing
# would ever call it.
Diversity.counttypes(demand::Demand) = length(demand.resource)
# Every member covers the same species, so the count of the whole is its first member's.
function Diversity.counttypes(demand::SpeciesRequirementCollection{Resource})
    return counttypes(first(values(demand)))
end

# ══ The materialised layers ════════════════════════════════════════════════════════════════════════
function Diversity.countsubcommunities(regime::ContinuousRegime)
    return length(regime.matrix)
end
function Diversity.countsubcommunities(regime::CategoricalRegime)
    return length(regime.matrix)
end
function Diversity.countsubcommunities(regime::LayerCollection{Condition})
    return countsubcommunities(first(values(regime)))
end
# On the public generic rather than the `_` hook: a supply is not an `AbstractPlaces`, so the hook
# would never be reached. The note in `Layer.jl` says the same of a regime.
function Diversity.countsubcommunities(supply::ContinuousLayer{Resource})
    return length(supply.matrix)
end
# Every sub-supply of a collection is on the same grid, so the count of the whole is the first's.
function Diversity.countsubcommunities(supply::LayerCollection{Resource})
    return countsubcommunities(first(values(supply)))
end

# ══ GridHabitat — the subcommunities ═══════════════════════════════════════════════════════════════
function Diversity.API._countsubcommunities(habitat::GridHabitat)
    return countsubcommunities(habitat.regime)
end
# A cell is named by **where it is** — its own half-open extent — so a row of a diversity table says
# which ground it describes rather than only which position it held.
#
# Computed on demand rather than stored. On a 1.2 million-cell grid the eager form costs about 33 MB
# against 8 bytes for this, and Diversity only ever `reshape`s or broadcasts over the names, which on
# a lazy vector gives a view rather than a copy. Nothing in the simulation itself reads a cell name.
function Diversity.API._getsubcommunitynames(habitat::GridHabitat)
    return CellNames(getcoords(habitat))
end

# ══ SpeciesList — the AbstractTypes hooks, every one delegating to `sl.types` ══════════════════════
function Diversity.API._gettypenames(sl::SpeciesList, input::Bool)
    return _gettypenames(sl.types, input)
end
function Diversity.API._counttypes(sl::SpeciesList, input::Bool)
    return _counttypes(sl.types, input)
end
# The second argument is a **scale**, not an array: the contract is set by `_calcordinariness`, which
# does `abundance, scale = _calcabundance(t, a)` and then `_calcsimilarity(t, scale)`. Annotating it
# `AbstractArray` matches nothing and falls silently through to Diversity's "not implemented".
function Diversity.API._calcsimilarity(sl::SpeciesList, scale::Real)
    return _calcsimilarity(sl.types, scale)
end
# What `EcoBase.view(eco; species = ...)` needs, and with it SpatialEcology's `groupspecies`.
# Delegating keeps the *kind* of the wrapped types where Diversity can preserve it, so a
# `UniqueTypes` subsets to a `UniqueTypes` rather than collapsing to a `GeneralTypes`. A **phylogeny**
# cannot be preserved: Diversity's generic rebuilds it as a `GeneralTypes` from the similarity matrix,
# so subsetting species costs the tree. Subsetting **sites** does not, because `view` leaves the types
# untouched when species are unrestricted, which is why `groupsites` keeps the whole `SpeciesList`.
function Diversity.API._subsettypes(sl::SpeciesList, idx, scale::Real)
    return _subsettypes(sl.types, idx, scale)
end
# A `SpeciesList` imposes no float constraint of its own: the `Float64` in
# `AbstractMetacommunity{Float64, ...}` is the **ecosystem's** and belongs there, while the wrapped
# types may well allow more.
Diversity.API.floattypes(sl::SpeciesList) = floattypes(sl.types)
Diversity.API._hassimilarity(sl::SpeciesList) = _hassimilarity(sl.types)
# Delegated rather than implemented, and that is what makes a phylogeny work. Diversity's generic
# `_calcordinariness(t::AbstractThings, a, ::Real)` recomputes `_calcabundance(t, a)`, but
# `Metacommunity` passes abundances that have already been through it, so re-applying them is a shape
# error for `PhyloBranches`. `DiversityPhyloExt` supplies the method that handles that, and it can
# only be reached by delegating.
function Diversity.API._calcordinariness(sl::SpeciesList, a::AbstractArray,
                                         scale::Real)
    return _calcordinariness(sl.types, a, scale)
end
# The extra output columns a wrapped phylogeny contributes. Reached as
# `getaddedoutput(m) = _getaddedoutput(_gettypes(m))`, so without these the `AbstractThings` defaults
# — `nothing`, and an empty `Dict` — silently drop them.
Diversity.API._getaddedoutput(sl::SpeciesList) = _getaddedoutput(sl.types)
Diversity.API._addedoutputcols(sl::SpeciesList) = _addedoutputcols(sl.types)
function Diversity.API._calcabundance(sl::SpeciesList, a::AbstractArray)
    return _calcabundance(sl.types, a)
end
function Diversity.API._getdiversityname(sl::SpeciesList)
    return _getdiversityname(sl.types)
end

# ══ AbstractEcosystem — the AbstractMetacommunity hooks ════════════════════════════════════════════
# **Diversity must be handed plain numeric data, never a labelled array.** Its internals are not
# robust to an arbitrary `DimArray`'s dims — reducing over `Dim{:location}` can fail inside its own
# reduction machinery. `GridLandscape.matrix` is a plain `Matrix{Int64}` precisely so that it can be
# passed straight through here; the labelled views are `dimmatrix`/`dimgrid` and must not be used at
# this boundary.
function Diversity.API._getabundance(eco::AbstractEcosystem, raw::Bool)
    abun = eco.abundances.matrix
    if raw
        return abun
    else
        return _calcabundance(_gettypes(eco), abun / sum(abun))[1]
    end
end
# The cached ecosystem answers from its most recent computed slice, and otherwise exactly as the
# method above: the same plain matrix, the same `_calcabundance`. It must forward through
# `_gettypes` for the same reason every other hook does — `SpeciesList` holds its `AbstractTypes` in
# `sl.types`, so a phylogeny's branch abundances only exist once `_calcabundance` has mapped the
# species onto them. Passing the species abundances straight through is invisible for a
# `UniqueTypes` list, which `_calcabundance` returns unchanged, and wrong for every phylogeny.
function Diversity.API._getabundance(cache::CachedEcosystem, raw::Bool)
    if all(ismissing.(cache.abundances.matrix))
        error("Abundances are missing")
    else
        id = findall(.!ismissing.(cache.abundances.matrix))[end]
        abun = cache.abundances.matrix[id]
    end

    abun = abun.matrix
    if raw
        return abun
    else
        return _calcabundance(_gettypes(cache), abun / sum(abun))[1]
    end
end
function Diversity.API._getpartition(eco::AbstractEcosystem)
    return eco.habitat
end
function Diversity.API._gettypes(eco::AbstractEcosystem)
    return eco.spplist
end
# Recomputed on every call, into preallocated scratch, rather than cached. A cached ordinariness has
# to be invalidated by everything that mutates the ecosystem, and a matrix whose species count is
# silently one short is not an error anything notices.
#
# Filling `cache.relativeabun` in place keeps almost all of the saving anyway: the chain costs 32
# bytes rather than gigabytes at UK scale, because `_calcabundance` and `_calcordinariness` pass a
# `UniqueTypes` array straight through instead of copying it.
#
# The shape check is what makes adding a species safe, and it is deliberately a check rather than a
# hook into `_addspecies!` — nothing then has to remember to keep the two in step.
function Diversity.API._getordinariness!(eco::AbstractEcosystem)
    abun = eco.abundances.matrix
    size(eco.cache.relativeabun) == size(abun) ||
        (eco.cache.relativeabun = similar(abun, Float64))
    eco.cache.relativeabun .= abun ./ sum(abun)
    processed, scale = _calcabundance(_gettypes(eco), eco.cache.relativeabun)
    return _calcordinariness(eco.spplist, processed, scale)
end
# The same call `_getabundance(eco, false)` makes, taking the scale instead of the abundances, which
# is how Diversity's own `Metacommunity` works: one `_calcabundance` gives both. It must be handed the
# **normalised** matrix rather than the raw counts, or the scale does not match the abundances it is
# paired with.
function Diversity.API._getscale(eco::AbstractEcosystem)
    abun = eco.abundances.matrix
    return _calcabundance(_gettypes(eco), abun / sum(abun))[2]
end
