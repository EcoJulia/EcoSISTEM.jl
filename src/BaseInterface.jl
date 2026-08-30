# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Conformance to `Base`, not our own API — every method this package defines on a `Base` or `Core`
# generic, gathered here rather than beside the type each one is about, and ordered **by type** so
# that "what does `Base` see when it is handed one of these?" can be answered by reading down. Its
# siblings are `DiversityInterface.jl` and `EcoBaseInterface.jl`.
#
# Every definition is **qualified**, `Base.eltype(...)`, so this file needs no `import` lines.
# `docs/overloads.md` is the cross-reference, listing every foreign generic this package extends and
# where.
#
# Five groups of `Base` methods are deliberately **not** here, each because moving it would split a
# conformance that belongs together:
#
#   - **`show` lives beside the type it prints**, in the file that declares it. How a type prints is
#     part of that type rather than a separate conformance: each method reads the struct's own fields
#     and its prose explains what those fields mean, so it is unreadable away from them.
#   - **`Dist.jl`** keeps `rand`, `minimum` and `maximum` on `Trapezoid` beside its `pdf` and
#     `params`. Those five *are* the `Distributions` interface, and which module owns each generic is
#     not the useful grouping.
#   - **`deprecations.jl`** keeps `eltype` for the shims it defines.
#   - **`ClimatePref`** keeps its own — `size`, `length`, `eltype`, `show`, `read`.
#   - **`ext/`** keeps `Base.read`'s dataset methods; an extension cannot be loaded early.
#
# Included last, which is safe rather than merely tidy: method definitions are order-independent, and
# nothing in `src/` calls one of these at load time.

# ═══════════════════════════════════════════════════════════════════════════════════════════════
# The member-container interface — `AbstractLayer`, `AbstractSpeciesRequirement`, `AbstractNicheFit`
# ═══════════════════════════════════════════════════════════════════════════════════════════════
#
# Every one of these families is "one member, or several members named by their niche axis" — which
# is exactly what a `NamedTuple` is. So rather than a bespoke accessor trio per family, the
# collections **forward the standard container interface to their backing** and the leaves implement
# the same interface as one-member containers. A caller who knows `NamedTuple` already knows this one,
# and a further family costs nothing.
#
# Two members of the interface are deliberately **absent**, and neither is an oversight:
#
#   - **`collect`** — not defined here, and not thereby prevented: `iterate` gives Base's generic
#     `collect` for free, so `collect(lc)` works and allocates into an abstractly-typed `Vector`.
#     There is no way to have iteration without it, so the position is *not encouraged, not defined,
#     still reachable*, and a caller who wants one is better served by `collect(values(x))`, which
#     says what it is doing.
#   - **`eltype`** — a **leaf** property, not a container one. A leaf's `eltype` is the unit frame its
#     data is in, and it is what supplies a nichefit's frame parameter; a collection has no single
#     one, so asking is an error rather than a guess. See the note in `SpeciesRequirement.jl`.
#
# The two sides of the interface. `_SingletonsAndCollections` is **not** just the singletons — it is
# the three abstract supertypes, so a `LayerCollection` matches it too. `_JustCollections` is a
# **subtype** of it, so a collection's methods are the more specific ones and win, which is what lets
# the methods over the wider union read as the plain one-member case.
const _JustCollections = Union{LayerCollection, SpeciesRequirementCollection,
                               CombiningFit}

const _SingletonsAndCollections = Union{AbstractLayer,
                                        AbstractSpeciesRequirement,
                                        AbstractNicheFit}

# --- The two base cases, and the only methods that need a leaf/collection split ------------------
# A leaf is named by its own axis, exactly as it is inside a collection, so a leaf and a collection of
# one answer identically and nothing downstream has to know which it holds.
Base.keys(x::_SingletonsAndCollections) = (nameof(axisof(x)),)

Base.keys(c::_JustCollections) = keys(_backing(c))

Base.values(x::_SingletonsAndCollections) = (x,)

Base.values(c::_JustCollections) = values(_backing(c))

# `length` keeps a split only to avoid building a `NamedTuple` in order to count one member.
Base.length(::_SingletonsAndCollections) = 1

Base.length(c::_JustCollections) = length(_backing(c))

Base.length(raster::ClimateRaster) = length(raster.array)

Base.length(c::LayerCache) = length(c.reads)

# --- Everything else: ONE generic method, over the backing ----------------------------------------
Base.iterate(x::_SingletonsAndCollections) = iterate(_backing(x))

Base.iterate(x::_SingletonsAndCollections, i::Int) = iterate(_backing(x), i)

function Base.getindex(x::_SingletonsAndCollections, i::Integer)
    return getindex(_backing(x), i)
end

function Base.getindex(x::_SingletonsAndCollections, name::Symbol)
    return getindex(_backing(x), name)
end

Base.getindex(names::CellNames, i::Int) = _cellname(names.grid, i)

function Base.haskey(x::_SingletonsAndCollections, name::Symbol)
    return haskey(_backing(x), name)
end

Base.haskey(c::LayerCache, spec::SourceSpec) = haskey(c.reads, ReadKey(spec))

function Base.get(x::_SingletonsAndCollections, name::Symbol, default)
    return get(_backing(x), name, default)
end

Base.pairs(x::_SingletonsAndCollections) = pairs(_backing(x))

Base.firstindex(::_SingletonsAndCollections) = 1

Base.lastindex(x::_SingletonsAndCollections) = length(x)

# The public spelling: `NamedTuple(lc)` reads as the conversion it is, and a leaf and a collection of
# one give the same answer.
Base.NamedTuple(x::_SingletonsAndCollections) = _backing(x)

# Named member access, routed through `getindex` so there is one implementation rather than two. That
# makes it depend on `_backing` using `getfield` — see the warning there, since with `getproperty`
# this line would recurse — and it keeps a member named `:nt` reachable. `propertynames` agrees with
# `keys`, because for a collection the properties **are** the members.
#
# Deliberately collection-only: for a leaf the properties are its own fields, `matrix`, `size`,
# `change`, and claiming otherwise would make `propertynames` lie and break tab-completion. `keys`
# already answers the container question for both.
Base.getproperty(c::_JustCollections, name::Symbol) = c[name]

Base.propertynames(c::_JustCollections) = keys(c)

Base.hasproperty(c::_JustCollections, name::Symbol) = haskey(c, name)

# `merge` rebuilds through each family's own sole constructor, so the result is checked exactly as a
# freshly written collection is and the role is re-read off the members.
#
# A repeated axis **replaces** rather than erroring, which is `NamedTuple`'s own `merge` semantics and
# the right ones here: merging is how one layer is swapped for another. Both arguments may be a leaf,
# so `merge(temperature, rainfall)` builds the two-layer collection while `merge(collection, extra)`
# extends one; the family comes from the **first** argument, and a second from another family is
# refused by the constructor rather than here.
function Base.merge(a::AbstractLayer, b::_SingletonsAndCollections)
    return _mergemembers(LayerCollection, a, b)
end

function Base.merge(a::AbstractSpeciesRequirement, b::_SingletonsAndCollections)
    return _mergemembers(SpeciesRequirementCollection, a, b)
end

function Base.merge(a::AbstractNicheFit, b::_SingletonsAndCollections)
    return CombiningFit(nichefitcombine(a), merge(_backing(a), _backing(b)))
end

# ══ ClimateRaster — the array and broadcast interface ══════════════════════════════════════════════
# The private helpers these call — `_unwrapped`, `_rastersof`, `_rewrap` — stay beside `ClimateRaster`
# itself: they are about rebuilding a raster, not about conforming to `Base`.
Base.size(raster::ClimateRaster) = size(raster.array)

function Base.size(regime::Union{CategoricalRegime, ContinuousRegime}, d)
    return size(regime.matrix, d)
end

function Base.size(regime::LayerCollection{Condition}, d)
    return size(first(values(regime)), d)
end

function Base.size(names::CellNames)
    return (length(names.grid.y) * length(names.grid.x),)
end

Base.eltype(::ClimateRaster{S, C, A}) where {S, C, A} = eltype(A)

# ══ The materialised layers ════════════════════════════════════════════════════════════════════════
function Base.eltype(regime::ContinuousRegime{C}) where {C}
    return C
end

function Base.eltype(regime::CategoricalRegime{D}) where {D}
    return D
end

Base.eltype(::ContinuousLayer{Resource, A, V}) where {A, V} = V

# ══ The species requirements ═══════════════════════════════════════════════════════════════════════
# The species-response type of a leaf. Declared on the two leaf-side abstracts rather than on
# `AbstractSpeciesRequirement`, which a collection would also match — see the note there.
Base.eltype(::AbstractCategoricalTolerance{A, V}) where {A, V} = V

Base.eltype(::ContinuousTolerance{A, V}) where {A, V} = V

Base.eltype(::Demand{A, V}) where {A, V} = V

# ══ The nichefits ══════════════════════════════════════════════════════════════════════════════════
function Base.eltype(nichefit::NicheSuitability{A, V}) where {A, V}
    return V
end

function Base.eltype(nichefit::CategoricalSuitability{A, V}) where {A, V}
    return V
end

function Base.eltype(nichefit::NoFitContinuous{A, V}) where {A, V}
    return V
end

function Base.eltype(nichefit::NoFitCategorical{A, V}) where {A, V}
    return V
end

# `Base.`-qualified in the bodies too: several of this module's dependencies export `axes`, which
# leaves the bare name ambiguous and therefore undefined here.
Base.axes(raster::ClimateRaster) = Base.axes(raster.array)

Base.ndims(::Type{<:ClimateRaster{S, C, A}}) where {S, C, A} = Base.ndims(A)

function Base.copy(bc::Broadcast.Broadcasted{ClimateRasterStyle})
    return _rewrap(copy(Broadcast.instantiate(_unwrapped(bc))), _rastersof(bc))
end

Base.hash(k::ReadKey, h::UInt) = hash(k.readkw,
                                      hash(k.code, hash(k.source, h)))

function Base.hash(fate::AbstractLayerFate, h::UInt)
    return foldl((acc, f) -> hash(getfield(fate, f), acc),
                 fieldnames(typeof(fate)), init = hash(typeof(fate), h))
end

Base.IndexStyle(::Type{<:CellNames}) = IndexLinear()

Base.BroadcastStyle(::Type{<:ClimateRaster}) = ClimateRasterStyle()

# A raster combined with anything — a scalar, a plain array, another raster — stays a raster. Only
# this ordering is needed: Julia tries the arguments the other way round before giving up.
function Base.BroadcastStyle(::ClimateRasterStyle, ::Broadcast.BroadcastStyle)
    return ClimateRasterStyle()
end

# Required, and easy to miss: a type that declares its own style must also say that it enters a
# broadcast **as itself**. Without this the default treats a raster as a scalar to be wrapped in a
# `Ref`, and the style is never consulted at all.
Base.broadcastable(raster::ClimateRaster) = raster

# `sum(bands)` over a tuple of rasters is the shape a multi-band combine wants, and it reduces with
# `+`, so defining `+` and its `-` partner here is what lets a spec write `sum(bands)` rather than
# reaching for `.array`. Both are broadcasts, so they inherit the style rule above rather than
# restating it.
Base.:+(a::ClimateRaster, b::ClimateRaster) = a .+ b
Base.:-(a::ClimateRaster, b::ClimateRaster) = a .- b

# ══ AbstractChangeSpec ═════════════════════════════════════════════════════════════════════════════
Base.:+(a::AbstractChangeSpec, b::AbstractChangeSpec) = CombinedChange(a, b)

# ══ The study-area report and its cache ════════════════════════════════════════════════════════════
# Value equality and hashing are defined explicitly rather than relying on the `===` fallback: the
# `readkw` `NamedTuple` may hold heap values (a month range, say), for which the default `objectid`
# hash is identity-based and would miss every cache hit.
function Base.:(==)(a::ReadKey, b::ReadKey)
    return a.source == b.source && a.code == b.code && a.readkw == b.readkw
end

# Value equality, because a fate is a value: two are the same if they are the same kind and say the
# same thing. Not the default — Julia falls back to `===` for a struct, which on `LayerResampled`
# compares its `reason` by **reference**, so two separately built reports would disagree about layers
# they in fact treat identically. `hash` follows, to keep the pair consistent.
function Base.:(==)(a::AbstractLayerFate, b::AbstractLayerFate)
    typeof(a) === typeof(b) || return false
    return all(f -> getfield(a, f) == getfield(b, f), fieldnames(typeof(a)))
end

# ══ DiversitySet ═══════════════════════════════════════════════════════════════════════════════════
"""
    append!(diversityset::DiversitySet, dat::DataFrame)

Append a `DataFrame` of diversity results to the data a [`DiversitySet`](@ref) already holds.

# Arguments

  - `diversityset`: the set to append to.
  - `dat`: the results, one row per subcommunity per timepoint.
"""
function Base.append!(diversityset::DiversitySet, dat::DataFrame)
    return append!(diversityset.data, dat)
end

# ---------------------------------------------------------------------------
# The backing
# ---------------------------------------------------------------------------
# The members as a `NamedTuple` — the **one** place the two sides of the interface differ. A
# collection *stores* one, reached with `getfield` so that no member name is reserved; a leaf *is* one
# member, so its backing is built from `keys` and `values`. Everything else above is then a single
# generic method over `_SingletonsAndCollections`, forwarding to whichever of these returns.
#
# Not circular: the leaf's `keys` and `values` are defined directly, and the collection's more
# specific method never reaches this one.
_backing(x) = NamedTuple{keys(x)}(values(x))

# **`getfield`, not `getproperty`**, and that is load-bearing for more than naming.
# `Base.getproperty(c::_JustCollections, ...)` above is defined as `c[name]`, which routes through
# `getindex` to **this** method and then to the backing `NamedTuple`. Were this to use `getproperty`,
# that path would call straight back into `getproperty` and recurse forever. The recursion is
# invisible from either site on its own — it exists only because the two are wired together — so
# change neither without re-reading the other.
_backing(c::_JustCollections) = getfield(c, :nt)

# Shared body of the two role-carrying `merge`s — the third needs its combining function too, so it
# cannot use this.
_mergemembers(::Type{C}, a, b) where {C} = C(merge(_backing(a), _backing(b)))
