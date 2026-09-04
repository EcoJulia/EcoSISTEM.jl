# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The materialised, hot-loop grid layers: what a cell *is like* (a `Condition` regime) or what it
# *provides* (a `Resource` supply), on one niche axis each, and several of them over one grid.

using Unitful

using DimensionalData

using RecipesBase

using Unitful.DefaultSymbols

using DimensionalData.Lookups: NoLookup

using StatsBase

using Diversity

using Diversity.API

using EcoSISTEM.Units

using EcoBase

using Measures: AbsoluteLength

using Dates: Dates

"""
    ContinuousLayer{R <: Role, A <: NicheAxis, V <: Number, Arr <: DimensionalData.AbstractDimArray{V, 2}, S <: Unitful.Quantity}

A continuous (numeric) grid layer of role `R` on niche axis `A`, holding a value `matrix`
(eltype `V`) as a `DimArray` over `(Y, X)` - real `Projected`/`Sampled` lookups for a
data-driven source, `NoLookup` for a synthetic one (see `data/architecture.md`).

A layer holds one grid of values, never a stack of them: a time-varying layer holds the values
current *now* and carries a [`SeriesLayerChange`](@ref) as its `change`, which decides from elapsed time
which stored slice that is. The array bound says `2` to keep it that way - the stack and the
cursor that walked it were what made a layer's own state ambiguous.
"""
mutable struct ContinuousLayer{R <: Role, A <: NicheAxis, V <: Number,
                               Arr <:
                               DimensionalData.AbstractDimArray{V, 2},
                               S <: Unitful.Quantity} <:
               AbstractLayer{R, A}
    matrix::Arr
    size::S
    change::AbstractLayerChange
end

"""    ContinuousRegime{V} - a static continuous regime (a `Condition`-role [`ContinuousLayer`](@ref) over a `(Y, X)` `DimArray{V}`). """
const ContinuousRegime{V} = ContinuousLayer{Condition, A, V, Arr,
                                            S} where {A, S,
                                                      Arr <:
                                                      DimArray{V, 2,
                                                               <:Tuple{<:Y,
                                                                       <:X}}}

# ---------------------------------------------------------------------------
# Back-compat aliases + constructors (regime role)
# ---------------------------------------------------------------------------

# Old positional constructors -> new layer types, defaulting to the `NicheAxis` axis (the
# axis-aware `materialise` path constructs `ContinuousLayer{Condition, A}` with the real axis
# directly).
function ContinuousRegime(matrix::DimArray{C, 2, <:Tuple{<:Y, <:X}},
                          change::AbstractLayerChange) where {C}
    _checkhascoords(matrix)
    size = _derivecellsize(matrix)
    return ContinuousLayer{Condition, NicheAxis, C, typeof(matrix),
                           typeof(size)}(matrix, size, change)
end

# A plain `Matrix` (no real CRS to attach - the synthetic case) is wrapped onto a fresh
# `NoLookup` `(Y, X)` grid, mirroring the `Supply{A}` constructors below.
# **`Unitful.Quantity`, not `Unitful.Length`, and only because a GEOGRAPHIC grid's cell size is
# genuinely an angle** (`0.1666...°`). Nothing is being converted: a degree grid really has degree
# cells, and forcing a metric equivalent here would fabricate a number the grid does not have - the
# area-preserving size exists to let the *simulator* pretend cells are uniform, which is exactly the
# pretence `build_ecosystem` refuses to make. The layer's `S` parameter keeps the field concrete
# and dispatchable, so anything that needs a metric grid can require `S <: Unitful.Length`.
function ContinuousRegime(matrix::Matrix{C}, size::Unitful.Quantity,
                          change::AbstractLayerChange) where {C}
    return ContinuousRegime(_sizedyx(matrix, size), change)
end

"""
    Supply{A}

A supply of the resource measured on niche axis `A` - a `Resource`-role
[`ContinuousLayer`](@ref) over a `(Y, X)` grid: `Supply{SolarRadiation}` (`kJ/day` per cell),
`Supply{Precipitation}` (`L/day`), `Supply{CarbonFlux}` (`g/day`).

The axis is the *only* declaration of what a supply measures: its unit comes from
`canonicalunit(Resource, A)`, and an axis that declares no resource cannot be built as a supply at
all. Construct one from a matrix of per-cell rates, `Supply{SolarRadiation}(fill(200.0kJ/day, 10,
10))`, or let [`GridHabitat`](@ref) do it from a spec that names the axis.

A supply that varies in time is the same type, carrying a [`SeriesLayerChange`](@ref).
"""
const Supply{A} = ContinuousLayer{Resource, A, V, Arr,
                                  S} where {V, S,
                                            Arr <:
                                            DimArray{V, 2, <:Tuple{<:Y, <:X}}}

# There is deliberately no `_getavailablesupply` here. Totalling a supply needs the habitat's
# `active` mask - an inactive cell's resource is not available to anything - so it belongs with the
# habitat, and `totalsupply` (`GridHabitat.jl`) is the whole of it. A supply-only version
# existed, was called by nothing, and quietly became `sum(supply.matrix)` tested against
# `sum(supply.matrix)`.

# --- Constructors ---------------------------------------------------------------------
# **One pair of methods for every resource axis**, because the axis says everything that differs:
# `_canonicalresource` refuses an axis that is not a resource and states the values in that axis's
# canonical unit, so `MJ/day` is accepted and stored as `kJ/day`, and a wrong dimension is refused
# **here**, at the call that made the mistake, rather than as a `DimensionError` in the hot loop. A
# new resource axis needs no constructor of its own.
#
# There is no free or dimensionless supply, because such a thing could not name an axis, and an axis
# is where a layer's meaning comes from.
#
# The `<:Tuple{<:Y, <:X}` in the signature is load-bearing, not decoration: it is what rejects a
# wrongly-ordered `(X, Y)` array outright at construction, and so what stops the historical
# dimension-order bug recurring.
#
# **A bare `Matrix` is not accepted here.** It carries no dims, so there is nothing for a supply's
# cell size to be derived *from*, and the alternative is a placeholder every supply reports whatever
# its grid. A `Matrix` method exists for the released `SolarBudget(matrix)` path and lives in
# `src/deprecations.jl`; new code passes a `DimArray` that states its own grid.
function Supply{A}(mat::DimArray{V, 2,
                                 <:Tuple{<:Y, <:X}}) where {A <:
                                                            NicheAxis,
                                                            V}
    _checkhascoords(mat)
    canon = _canonicalresource(mat, A)
    # Derived from the array's own coordinates, never stamped with a placeholder: a supply that
    # reports a cell size it did not measure is worse than one that cannot answer.
    size = _derivecellsize(canon)
    return ContinuousLayer{Resource, A, eltype(canon), typeof(canon),
                           typeof(size)}(canon, size, NoLayerChange())
end

"""
    CategoricalLayer{A <: NicheAxis, V, Arr <: DimensionalData.AbstractDimArray{V, 2}, S <: Unitful.Quantity}

A categorical (class-code) grid layer on niche axis `A` - always a `Condition` (there is no
categorical supply). `matrix` is a `DimArray` over `(Y, X)`.
"""
mutable struct CategoricalLayer{A <: NicheAxis, V,
                                Arr <: DimensionalData.AbstractDimArray{V, 2},
                                S <: Unitful.Quantity} <:
               AbstractLayer{Condition, A}
    matrix::Arr
    size::S
    change::AbstractLayerChange
end

"""    CategoricalRegime{V} - a categorical regime (a [`CategoricalLayer`](@ref), e.g. land cover), a `(Y, X)` `DimArray{V}`. """
const CategoricalRegime{V} = CategoricalLayer{A, V, Arr,
                                              S} where {A, S,
                                                        Arr <:
                                                        DimArray{V, 2,
                                                                 <:Tuple{<:Y,
                                                                         <:X}}}

function CategoricalRegime(matrix::DimArray{D, 2, <:Tuple{<:Y, <:X}},
                           change::AbstractLayerChange) where {D}
    _checkhascoords(matrix)
    size = _derivecellsize(matrix)
    return CategoricalLayer{NicheAxis, D, typeof(matrix),
                            typeof(size)}(matrix, size, change)
end

function CategoricalRegime(matrix::Matrix{D}, size::Unitful.Quantity,
                           change::AbstractLayerChange) where {D}
    return CategoricalRegime(_sizedyx(matrix, size), change)
end

"""
    LayerCollection{R <: Role, A, C <: NamedTuple}

Several layers of the same role `R` over one grid (e.g. temperature + rainfall) and over the
**axis structure** `A`, held by name in a `NamedTuple`, `(Temperature = ..., Precipitation = ...)`.
Sub-layer types stay concrete in the backing's own type, so the hot loop is exactly as type-stable
as the fixed-arity named fields it replaces, and the arity is limited only by what the caller
writes.

**Layers are named by their axis**, so a `Tuple` is accepted and named for you; two layers on
the **same** axis cannot be told apart that way and are refused, asking for explicit names.

Layers are reached through the standard container interface - `lc.Precipitation`, `lc[1]`,
`lc[:Precipitation]`, `keys`, `values`, `pairs`, `iterate`, `length`, `merge` and `NamedTuple(lc)`,
all forwarded to the backing (`src/collections.jl`). A **single layer answers identically**, as a
one-member container.
"""
struct LayerCollection{R <: Role, A, C <: NamedTuple} <:
       AbstractLayer{R, A}
    nt::C

    # The sole constructor, so there is exactly one spelling and `R` cannot be got wrong: the role
    # is *read off* the layers rather than chosen, so leaving the generated positional constructors
    # alive would let a caller assert a role its own layers contradict.
    # Takes a `Tuple` or a `NamedTuple` as it always did - a tuple's names are **derived** here
    # (`_asnamedtuple`, `collections.jl`) rather than the two backings being kept apart downstream.
    # `A` is a third parameter rather than being computed in the supertype expression, because
    # Julia refuses that: `<: AbstractLayer{R, first(fieldtypes(C))}` is a `MethodError` on a
    # `TypeVar`. So the constructor derives it, which keeps the "read off the members" guarantee.
    function LayerCollection(layers::Union{Tuple, NamedTuple})
        nt = _asnamedtuple(layers)
        return new{_sharedrole(AbstractLayer, values(nt), "layer"),
                   _axisstructure(values(nt)), typeof(nt)}(nt)
    end
end

# == Functions ==================================================================================

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# A layer is routinely nested - inside a habitat, inside a collection, inside a vector - so it needs
# the compact one-liner as much as the multi-line form. Without one the default prints the whole
# matrix: measured at 50 016 characters for a 60 x 100 regime and 122 007 for a two-member
# collection, on a single line.
#
# Every method here is constant in the grid size, and none reduces over the data. A `show` runs
# every time a value is printed, so a `sum` or an `extrema` over the matrix would make displaying a
# layer cost as much as using it.
#
# The role decides the spelling, because it decides what a *user* writes. A supply is named
# `Supply{SolarRadiation}`, so its axis belongs in the braces; a regime is reached through the
# `ContinuousRegime` alias, so its axis is information rather than part of the name.
function Base.show(io::IO, l::ContinuousLayer{Condition})
    ny, nx = size(l.matrix)
    return print(io,
                 "ContinuousRegime($(nameof(axisof(l))), $(ny) × $(nx), $(unit(eltype(l))))")
end

function Base.show(io::IO, l::ContinuousLayer{Resource})
    ny, nx = size(l.matrix)
    return print(io,
                 "Supply{$(nameof(axisof(l)))}($(ny) × $(nx), $(unit(eltype(l))))")
end

# No unit: a categorical layer holds class codes, and the axis is the whole of what they mean.
function Base.show(io::IO, l::CategoricalLayer)
    ny, nx = size(l.matrix)
    return print(io, "CategoricalRegime($(nameof(axisof(l))), $(ny) × $(nx))")
end

# One member, as it appears inside a collection's summary: its axis, and its unit where it has one.
# The stored name is shown only when it differs from the axis, since a derived name repeats it.
function _membersummary(name::Symbol, l::AbstractLayer)
    axis = nameof(axisof(l))
    u = l isa CategoricalLayer ? "" : " $(unit(eltype(l)))"
    return name === axis ? "$(axis)$(u)" : "$(name): $(axis)$(u)"
end

function Base.show(io::IO, c::LayerCollection{R}) where {R}
    ny, nx = size(first(values(c)).matrix)
    members = join(map(_membersummary, keys(c), values(c)), " + ")
    return print(io, "LayerCollection{$(nameof(R))}($(members), $(ny) × $(nx))")
end

function Base.show(io::IO, ::MIME"text/plain", l::AbstractLayer)
    ny, nx = size(l.matrix)
    println(io, sprint(show, l))
    println(io, "  axis      ", axisof(l))
    l isa CategoricalLayer || println(io, "  unit      ", unit(eltype(l)))
    println(io, "  grid      $(ny) × $(nx) cells")
    print(io, "  change    ", nameof(typeof(l.change)))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", c::LayerCollection)
    println(io, sprint(show, c))
    for (name, l) in pairs(c)
        println(io, "  ", rpad(name, 14), sprint(show, l))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Supplies - `Resource`-role layers (folded onto `ContinuousLayer{Resource}`)
# ---------------------------------------------------------------------------
# `Supply{A}` (defined in Layer.jl) is an alias over `ContinuousLayer{Resource, A, V, Arr}`, always
# over `(Y, X)` and with the value type left free - the axis says what the unit is. The
# constructors below fill the (unused) `size` and the per-timestep `change` rule and zero
# NaNs, reproducing the old supply structs. A supply's `size` is never read
# (geometry/dispersal use the regime), so a placeholder is stored. A supply built here never
# changes; one read from a monthly stack carries a `SeriesLayerChange` installed by `_setseries!`.
const _SUPPLY_SIZE = 1.0m

const px = AbsoluteLength(0.254)

"""
    setchange!(layer::AbstractLayer, spec)

Install `spec` as `layer`'s per-timestep change, converting and validating its values against the
layer exactly once. `spec` is a change recipe (e.g. `IncrementBy(0.02K/yr)`) or an already
materialised [`AbstractLayerChange`](@ref).
"""
function setchange!(layer::AbstractLayer, spec)
    return layer.change = _attachchange(spec, layer)
end

# ---------------------------------------------------------------------------
# Coverage: does the run fit the series driving it?
# ---------------------------------------------------------------------------
# Reported before the first step rather than discovered during it. `ErrorAtEnd` would fail anyway,
# but at the step it happens - three years into a fifty-year run, after three years of compute - so
# saying so up front is the whole value. `HoldAtEnd` never fails, which is exactly why it earns a
# warning instead: a series pinned at its last slice for most of a run is rarely what was meant, and
# the *proportion* is the number that tells someone whether they meant it. `RepeatAtEnd` cycles by
# design and has nothing to report.

"""
    checkcoverage(eco::AbstractEcosystem, duration::Unitful.Time, timestep::Unitful.Time)

Check that every stored series in `eco`'s environment covers a run of `duration` in steps of
`timestep`, erroring or warning as each series' own [`AbstractSeriesEnd`](@ref) policy dictates.
Called by [`simulate!`](@ref) before the first timestep.
"""
function checkcoverage(eco::AbstractEcosystem, duration::Unitful.Time,
                       timestep::Unitful.Time)
    final = _finalelapsed(eco, duration, timestep)
    foreach(l -> _checkcoverage(l, final),
            vcat(_coveredlayers(eco.habitat.regime),
                 _coveredlayers(eco.habitat.supply)))
    return nothing
end

# ---------------------------------------------------------------------------
# Bounds, checked before the run rather than during it
# ---------------------------------------------------------------------------
# Most out-of-bounds conditions are **predictable**: a steady rate over a known duration ends at
# `initial + rate × duration`, and an absolute series' values are all in hand. So they can be reported
# before the first timestep, exactly as `checkcoverage` reports a series that will not cover the run
# - which is the difference between "shorten your run" and a simulation that dies three years in.
#
# Deliberately **partial, and silent where it cannot be sure.** A `PatternedLayerChange` takes an
# arbitrary user function of elapsed time, so its reach is not bounded by its amplitude in general; a
# rate-valued series integrates. Those are left to `_enforcebounds!` at the write site rather than
# guessed at, because a *wrong* pre-flight refusal is worse than a late true one.

"""
    check_bounds(eco::AbstractEcosystem, duration::Unitful.Time, timestep::Unitful.Time)

Check that no change will drive a condition layer outside its axis' physical range
([`EcoSISTEM.bounds`](@ref)) during a run of `duration` in steps of `timestep`, erroring before the
first timestep if one would. Called by [`simulate!`](@ref).

Only changes whose reach can be computed exactly are checked; anything else is caught at the moment
it happens instead.
"""
function check_bounds(eco::AbstractEcosystem, duration::Unitful.Time,
                      timestep::Unitful.Time)
    final = _finalelapsed(eco, duration, timestep)
    foreach(_coveredlayers(eco.habitat.regime)) do layer
        return _checkreach(layer, final)
    end
    return nothing
end

# The resource available in each cell. A supply holds one grid of values whether or not it varies
# in time - a time-varying one carries a `SeriesLayerChange` that writes the current slice into that
# same matrix each step - so there is nothing here to choose between.
function _getsupply(supply::ContinuousLayer{Resource})
    return supply.matrix
end

function _getsupply(supply::LayerCollection{Resource}, name::Symbol)
    return _getsupply(getproperty(supply, name))
end

# --- Recipes --------------------------------------------------------------------------
@recipe function f(B::ContinuousLayer{Resource, A, V}) where {A, V}
    b = ustrip.(B.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> _resourcetitle(A)
    clims --> (minimum(b) * 0.99, maximum(b) * 1.01)
    return b
end

@recipe function f(H::LayerCollection{Resource})
    layers = values(H)
    layout := length(layers)
    for (i, s) in enumerate(layers)
        @series begin
            subplot := i
            s
        end
    end
end

# One recipe for both kinds of regime, where there were three. A categorical layer's codes and a
# continuous layer's values are both `ustrip`ped onto the same heatmap; nothing here needed to
# distinguish them, and the split is what let the categorical case rot unnoticed.
# The axis comes from `_layeraxis(H)` rather than a signature parameter: the `*Regime` aliases pin
# the *value* type and leave the axis free, so it cannot be destructured here.
@recipe function f(H::Union{ContinuousRegime, CategoricalRegime})
    h = ustrip.(H.matrix)
    seriestype := :heatmap
    grid --> false
    right_margin --> 0.1px
    margin --> 10px
    aspect_ratio --> 1
    title --> _regimetitle(_layeraxis(H))
    clims --> (minimum(h) * 0.99, maximum(h) * 1.01)
    ax = _plotaxes(H)
    return ax.x, ax.y, h
end

@recipe function f(H::LayerCollection{Condition})
    layers = values(H)
    layout := length(layers)
    for (i, r) in enumerate(layers)
        @series begin
            subplot := i
            r
        end
    end
end

# The methods below dispatch on the `Condition`-role layer types - `AbstractLayer{Condition}`,
# `ContinuousLayer`, `CategoricalLayer` and `LayerCollection` - which the `*Regime` aliases name. The
# released `*Hab` and `HabitatCollection2` spellings are deprecated aliases over the same types.

# **One level, not two.** A layer is **not** an `EcoBase.AbstractPlaces`, so Diversity's own
# `countsubcommunities` could never dispatch here and its `_countsubcommunities` hook would never be
# called - the indirection would buy nothing. These are methods on Diversity's *public* generic for
# our type, which is the whole extension mechanism. Contrast `GridHabitat`, which genuinely **is**
# an `AbstractPartition`: there the `_` hook is right, and the public name comes free.
# **A layer is not a grid**, and there are deliberately no `EcoBase.AbstractGrid` methods on one. It
# holds values *on* a grid, and cannot answer where its cells are without inventing an origin or
# dividing its stored size by a unit it may not be in - a geographic layer asked for a metric cell
# size can only answer `1.0 ° km^-1`, which is neither a length nor an angle. [`StudyGrid`](@ref)
# answers those from the grid's own dimensions instead, and a habitat reaches it through its
# [`StudyArea`](@ref).

iscontinuous(::ContinuousRegime) = true

iscontinuous(regime::CategoricalRegime) = false

# A plot title for a regime, keyed on its **axis** - the mirror of `_resourcetitle` for the other
# role, and the same rule as everything else on this path. The fallback names any axis from its own
# declaration, so a layer on a newly declared axis plots sensibly with no entry added anywhere.
#
# **This replaces a `Dict(K => "Temperature (K)", mm => "Rainfall (mm)", ...)` keyed on the value's
# unit** - the exact inference this subproject exists to remove, surviving in the one place nobody
# looked. It was broken in three separate ways, all of which the axis form cannot reproduce: the
# categorical recipe read `unitdict` from a *sibling recipe's local scope*, so plotting any
# categorical layer threw `UndefVarError`; the continuous one threw `KeyError` for any unit outside
# those four, which is most of them; and the temperature one converted to °C and then titled the
# result "Temperature (K)" whenever the stored unit was K.
# No `Temperature` entry: the fallback already produces exactly "Temperature (K)" from the axis's
# own name and canonical unit, so one would be a second copy of a fact the declaration already holds.
_regimetitle(::Type{<:Precipitation}) = "Rainfall (mm/day)"

_regimetitle(::Type{<:SolarRadiation}) = "Solar Radiation (kJ/m^2/day)"

function _regimetitle(::Type{A}) where {A <: NicheAxis}
    u = canonicalunit(A)
    return (isnothing(u) || u == NoUnits) ? "$(nameof(A))" :
           "$(nameof(A)) ($u)"
end

# The `(x, y)` coordinate vectors a heatmap wants, read off the layer's own dimensions - `x` named
# first because that is the order `plot(x, y, z)` takes, while `z` stays `(Y, X)` as everywhere else.
#
# **This replaces `xrange(H), yrange(H)`**, EcoBase derivations that worked only while a layer
# claimed to be a grid. They were built from a fabricated origin of `0` and a cell size divided by
# `km`, so a geographic layer was plotted against axes that started nowhere in particular and were
# labelled in `° km^-1`. The dims have said where the cells are all along.
function _plotaxes(layer::AbstractLayer)
    yx = _yx(layer)
    return (x = parent(DimensionalData.lookup(yx[2])),
            y = parent(DimensionalData.lookup(yx[1])))
end

# The one cell side dispersal is expressed against. Reads the grid rather than any stored field, and
# in the grid's own length units, so a projected grid keeps `2500.0 m` rather than being converted to
# kilometres behind the caller's back.
#
# `s.y` alone is not a shortcut: `_checksimulatable` refuses a non-square grid before any ecosystem
# exists, so `s.y` and `s.x` are equal here by construction. Dispersal's one-uniform-cell assumption
# is enforced at the boundary rather than papered over with a geometric mean at the point of use,
# which `sqrt(dy*dx)` would do invisibly.
function _uniformcellside(regime::Union{CategoricalRegime, ContinuousRegime})
    s = getcellsizes(regime)
    isnothing(s) &&
        error("this regime has no cell size, so dispersal cannot be set up on it - " *
              "`_checksimulatable` should have refused it first.")
    side = first(s.y)
    side isa Unitful.Length ||
        error("this regime's cells are $(Unitful.unit(side)) across, not a length, so dispersal " *
              "cannot be set up on it - `_checksimulatable` should have refused it first.")
    return side
end

# Every sub-layer of a collection is on the same grid (`_yx` checks it at construction), so the
# geometry of the whole is the geometry of the first.
function _uniformcellside(regime::LayerCollection{Condition})
    return _uniformcellside(first(values(regime)))
end

# Function to create a regime from a categorical set of types according to the
# Saura-Martinez-Millan algorithm (2000)
function _percolate!(M::AbstractMatrix, clumpiness::Real)
    for i in eachindex(M)
        if rand(Uniform(0, 1)) < clumpiness
            M[i] = 1
        end
    end
end

# Function to create clusters from percolated grid
#
# **The loop runs `(y, x)` - row, then column - and `_getneighbours` is called in that order.**
# Handing it the column first clusters each cell with the neighbours of its transpose, which a square
# grid cannot distinguish and any other rejects with *"Coordinates outside grid"*. The inner
# `mapslices` closures take `n` rather than `x`, because an argument named `x` here would shadow the
# column index - the same confusion one level down.
function _identifyclusters!(M::AbstractMatrix)
    # Begin cluster count
    count = 1
    # Loop through each grid square in M
    for y in Base.axes(M, 1)
        for x in Base.axes(M, 2)

            # If square is marked as 1, then apply cluster finding algorithm
            if M[y, x] == 1.0
                # Find neighbours of M at this location
                neighbours = _getneighbours(M, y, x)
                # Find out if any of the neighbours also have a value of 1, thus, have
                # not been assigned a cluster yet
                cluster = vcat(mapslices(n -> M[n[1], n[2]] .== 1, neighbours,
                                         dims = 2)...)
                # Find out if any of the neighbours have a value > 1, thus, have already
                # been assigned a cluster
                already = vcat(mapslices(n -> M[n[1], n[2]] .> 1, neighbours,
                                         dims = 2)...)
                # If any already assigned neighbours, then assign the grid square to this
                # same type
                if any(already)
                    neighbours = neighbours[already, :]
                    M[y, x] = M[neighbours[1, 1], neighbours[1, 2]]
                    # If none are assigned yet, then create a new cluster
                else
                    count = count + 1
                    neighbours = neighbours[cluster, :]
                    M[y, x] = count
                    map(i -> M[neighbours[i, 1], neighbours[i, 2]] = count,
                        Base.axes(neighbours, 1))
                end
            end
        end
    end
end

# `(y, x)` - row then column - for the same reason as `_identifyclusters!` above, and it had the
# same transposed call.
# **`assigned` is an EXPLICIT mask, and it must be** - this used `isassigned(T, y, x)`, which
# for an `Array{Int64}` is **always `true`**. `isassigned` reports whether a *reference* slot has been
# filled, so it means something for an `Array{Any}` and nothing at all for a concrete bits eltype:
# every index of `Array{Int64}(undef, ...)` is "assigned" the moment the array exists.
# Used that way, **every neighbour counts, including cells never written**, and the tally is taken
# over **uninitialised memory**: the same seed, single-threaded, gives several different maps across
# a dozen runs while the RNG stream itself stays byte-identical, so the irreproducibility is not in
# the seeding.
#
# **A sweep over grid shapes cannot catch this**, because a corrupted result is still a member of the
# class list and so still passes any "every value is a valid class" assertion. The property that
# actually breaks is determinism from a seed, which is what the regression testset asserts.
# It *looks* like a load problem, since it grows commoner when other test sets run beside it
# - but a single controlled run reproduced it with no load on one thread. The correlation was real
# and the causation was not.
function _fillin!(T, M, types, wv, assigned::AbstractMatrix{Bool})
    # Loop through grid of clusters
    for y in Base.axes(M, 1)
        for x in Base.axes(M, 2)
            # If square is zero then it is yet to be assigned
            if M[y, x] == 0
                # Find neighbours of square on string grid
                neighbours = _getneighbours(T, y, x, 8)
                # Check if they have already been assigned
                already = vcat(mapslices(n -> assigned[n[1], n[2]],
                                         neighbours, dims = 2)...)
                # If any already assigned then sample from most frequent neighbour types
                if any(already)
                    neighbours = neighbours[already, :]
                    # Find all neighbour types
                    neighbour_traits = map(i -> T[neighbours[i, 1],
                                                  neighbours[i, 2]],
                                           Base.axes(neighbours, 1))
                    # Find which one is counted most often
                    ind = argmax(map(t -> sum(neighbour_traits .== t), types))
                    # Assign this type to the grid square in T
                    T[y, x] = types[ind]
                    # If none are assigned in entire grid already,
                    # sample randomly from types
                elseif all(M .<= 1)
                    T[y, x] = sample(types, wv)
                    # If some are assigned in grid, sample from these
                else
                    T[y, x] = sample(T[M .> 1])
                end
                # Whichever branch wrote it, the cell is now assigned - this is what the old
                # `isassigned` was trying and failing to say.
                assigned[y, x] = true
            end
        end
    end
end

# The class codes themselves: a `dimension` grid of integer niche `types` with relative `weights`
# and spatial clumpiness `clumpiness`. **Values only, and deliberately no cell size** - the
# percolation and clustering depend on the grid's *shape* alone, so a cell size cannot affect the
# pattern. It was previously threaded in only to build the `CategoricalRegime` below, which
# `_specfield` then discarded via `.matrix`; that vestigial argument was annotated
# `::Unitful.Length` and so refused a geographic grid's angular cell size once `[GEO-SIZE]` made it
# honest. Splitting the field out is what lets both materialise paths reach it.
function _nichefield(dimension::Tuple,
                     types::Vector{Int64},
                     clumpiness::Float64,
                     weights::Vector)
    # Check that the proportion of coverage for each type matches the number
    # of types and that they add up to 1
    length(weights) == length(types) ||
        error("There must be a weight for each type")
    sum(weights) == 1 || error("Weights of regimes must sum to 1")
    # Create weighting from proportion regimes
    wv = Weights(weights)

    # Create an empty grid of the right dimension
    M = zeros(dimension)

    # If the dimensions are too small for the algorithm, just use a weighted sample
    if dimension[1] <= 2 || dimension[2] <= 2
        return sample(types, Weights(weights), dimension)
    end
    # Percolation step
    _percolate!(M, clumpiness)
    # Select clusters and assign types
    _identifyclusters!(M)
    # Create a string grid of the same dimensions
    T = Array{Int64}(undef, dimension)
    # Fill in T with clusters already created
    map(x -> T[M .== x] .= sample(types, wv), 1:maximum(M))
    # The loop above wrote exactly the cells the clustering labelled, i.e. `M >= 1`; everything
    # else is still uninitialised. Saying so explicitly is what `isassigned` could not.
    assigned = M .>= 1
    # Fill in undefined squares with most frequent neighbour
    _fillin!(T, M, types, wv, assigned)
    return T
end

# A `CategoricalRegime` of dimension `dimension`, made up of integer niche `types` with relative
# `weights` and spatial clumpiness `clumpiness`, on cells of `gridsquaresize`. Its only callers
# now are the two deprecated builders in `deprecations.jl`, which want the wrapper and supply a
# genuine metric cell size; the live paths take `_nichefield` and place it on real dims themselves.
function _randomniches(dimension::Tuple,
                       types::Vector{Int64},
                       clumpiness::Float64,
                       weights::Vector,
                       gridsquaresize::Unitful.Length)
    return CategoricalRegime(_nichefield(dimension, types, clumpiness,
                                         weights),
                             gridsquaresize, NoLayerChange())
end

"""
    simpleregime(val::Unitful.Quantity, size::Unitful.Length,
    dim::Tuple{Int64, Int64}, axis::Type{<:NicheAxis})

Create a [`ContinuousRegime`](@ref) regime of dimension `dim`, with cell `size`
and filled value, `val`, on niche axis `axis`.
"""
function simpleregime(val::Unitful.Quantity, size::Unitful.Length,
                      dim::Tuple{Int64, Int64},
                      axis::Type{A}) where {A <: NicheAxis}
    M = _sizedyx(fill(val, dim), size)
    return ContinuousLayer{Condition, A, typeof(val), typeof(M),
                           typeof(size)}(M, size, NoLayerChange())
end

"""
    simpleregime(val::Float64, size::Unitful.Length, dim::Tuple{Int64, Int64},
    axis::Type{<:NicheAxis})

Create a dimensionless [`ContinuousRegime`](@ref) filled with `val`.
"""
function simpleregime(val::Float64, size::Unitful.Length,
                      dim::Tuple{Int64, Int64},
                      axis::Type{A}) where {A <: NicheAxis}
    M = _sizedyx(fill(val, dim), size)
    return ContinuousLayer{Condition, A, typeof(val), typeof(M),
                           typeof(size)}(M, size, NoLayerChange())
end

# == Functions ==================================================================================

# ---------------------------------------------------------------------------
# Driving the layers
# ---------------------------------------------------------------------------

# Apply one layer's own change, and recurse into a collection. Role-free - it dispatches on the
# layer's shape, not on whether it is a regime or a supply - so the same two methods drive both.
function _layerupdate!(layer::AbstractLayer, elapsed::Unitful.Time,
                       timestep::Unitful.Time)
    return _applychange!(layer.change, layer, elapsed, timestep)
end

function _layerupdate!(layer::LayerCollection, elapsed::Unitful.Time,
                       timestep::Unitful.Time)
    foreach(l -> _layerupdate!(l, elapsed, timestep), values(layer))
    return nothing
end

# The slice a series-carrying layer starts on: the first, which is what elapsed time zero selects.
_firstslice(stack) = stack[:, :, 1]

# Install a read stack as `layer`'s own time series. The layer holds one slice at a time - the first
# to begin with - and the series writes whichever slice is current from then on, so a stored series
# is a layer change like any other rather than a second, parallel mechanism.
#
# `RepeatAtEnd`, because a stored stack cycles: a twelve-month climatology is meant to repeat. The
# cycle is a true modulus, so the overshoot is kept rather than discarded.
function _setseries!(layer::AbstractLayer, stack)
    setchange!(layer,
               ReplaceWith(SeriesChange(stack, atend = RepeatAtEnd())))
    return layer
end

# The two ecosystem-level layer-update entry points. They are about layer change, and they live
# here with the rest of the layer machinery because `AbstractEcosystem` is declared in
# `Ecology.jl`, the package's first include, rather than eleven files after this one.
"""
    regimeupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)

Update the regime of an ecosystem for one timestep.
"""
function regimeupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return _layerupdate!(eco.habitat.regime, simulationtime(eco), timestep)
end

"""
    supplyupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)

Update the supply of an ecosystem for one timestep.
"""
function supplyupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return _layerupdate!(eco.habitat.supply, simulationtime(eco), timestep)
end

# The elapsed time the run actually finishes at - which is *not* `duration`, and the difference is a
# whole timestep. `simulate!` takes `length(0s:timestep:duration)` steps, a range that includes both
# ends, so a twelve-month run in one-month steps advances the clock thirteen times. It also starts
# from wherever the clock already is rather than from zero, since `simulate!` does not reset it.
# Checking against `duration` instead would pass runs that then failed mid-flight, which is the
# exact failure this check exists to pre-empt.
function _finalelapsed(eco::AbstractEcosystem, duration::Unitful.Time,
                       timestep::Unitful.Time)
    steps = length((zero(duration)):timestep:duration)
    return simulationtime(eco) + steps * uconvert(s, float(timestep))
end

# The layers of a habitat, flattened, so a collection is walked member by member.
_coveredlayers(layer::AbstractLayer) = [layer]

function _coveredlayers(layer::LayerCollection)
    return mapreduce(_coveredlayers, vcat, values(layer),
                     init = AbstractLayer[])
end

# A layer covers a run if its change does, so the question is forwarded rather than answered here.
# The layer method exists so a caller can ask the thing it holds, not the thing inside it.
function _checkcoverage(layer::AbstractLayer, duration)
    return _checkcoverage(layer.change,
                          duration)
end

_checkcoverage(::AbstractLayerChange, _) = nothing

function _checkcoverage(change::SumOfLayerChanges, final)
    foreach(p -> _checkcoverage(p, final), change.parts)
    return nothing
end

function _checkcoverage(change::SeriesLayerChange, final)
    return _reportcoverage(change.atend, change, final)
end

# A cycling series always covers the run, and one that reverts to its layer has been told explicitly
# what to do past the end. Neither has anything to report.
_reportcoverage(::RepeatAtEnd, _, _) = nothing

_reportcoverage(::RevertToLayer, _, _) = nothing

function _reportcoverage(::ErrorAtEnd, change::SeriesLayerChange, final)
    _seriesspare(change, final) >= zero(final) && return nothing
    return error("this series does not cover the run: it reaches " *
                 "$(uconvert(month_mean_duration, _seriesreach(change) - change.origin)) of " *
                 "elapsed time from where the run starts, but the run reaches " *
                 "$(uconvert(month_mean_duration, final)). Use `atend = HoldAtEnd()` to " *
                 "keep the last slice, `atend = RepeatAtEnd()` to cycle, a longer series, or a " *
                 "shorter run. (An epoch later than the series' own start shortens what is left of " *
                 "it - see `build_ecosystem`.)")
end

function _reportcoverage(::HoldAtEnd, change::SeriesLayerChange, final)
    spare = _seriesspare(change, final)
    spare >= zero(final) && return nothing
    duration = final
    covered = duration + spare
    @warn "This series ends $(round(typeof(1.0 * month_mean_duration), uconvert(month_mean_duration, covered), digits = 1)) " *
          "into a $(round(typeof(1.0 * month_mean_duration), uconvert(month_mean_duration, duration), digits = 1)) run, " *
          "so `atend = HoldAtEnd()` holds its last slice for " *
          "$(round(100 * ustrip(NoUnits, -spare / duration), digits = 1))% of the simulation. " *
          "That is rarely intended: use `atend = RepeatAtEnd()` to cycle, a longer series, or a " *
          "shorter run."
    return nothing
end

# How much series is left over at the end of a run of `duration` - negative when the run outlasts it.
function _seriesspare(change::SeriesLayerChange, duration)
    return _seriesreach(change) -
           (change.origin + duration)
end

# ---------------------------------------------------------------------------
# Priming: the values a run actually starts with
# ---------------------------------------------------------------------------
# A layer holds its current values in `matrix` and the rule that drives them in `change`, and a
# change only reaches `matrix` when `_applychange!` runs - which `update!` does at the *end* of a
# timestep, after that step's population dynamics. Left alone, then, a series-driven layer spends the
# whole of timestep one holding whatever its builder left behind: the spec's own value, or (for a
# read stack) slice one regardless of what epoch the run was given. The dynamics of the first step
# would run against an environment no part of the model asked for.
#
# Priming closes that: `build_ecosystem` evaluates each change at the run's start, so the
# environment a simulation begins with is the one its series and epoch describe. It is done there
# rather than in the `Ecosystem` constructor because that constructor is also reached mid-run - a
# `CachedEcosystem` step rebuilds an ecosystem at an arbitrary elapsed time - where resetting layers
# to the start would be plainly wrong.

function _primeseries!(layer::AbstractLayer, elapsed::Unitful.Time)
    _primechange!(layer.change, layer, elapsed)
    return layer
end

function _primeseries!(layer::LayerCollection, elapsed::Unitful.Time)
    foreach(l -> _primeseries!(l, elapsed), values(layer))
    return layer
end

# Dispatch is on the change's *mode*, as `_applychange!` does, because what it means to evaluate a
# change at the start of a run differs by mode rather than by shape. Absolute and relative changes
# are pure functions of elapsed time, so asking what they say at the start is exactly what the first
# timestep would have written - and idempotent, since neither compounds on what it wrote last.
_primechange!(::AbstractLayerChange{NoChange}, ::AbstractLayer, _) = nothing

function _primechange!(change::AbstractLayerChange{AbsoluteChange},
                       layer::AbstractLayer, elapsed::Unitful.Time)
    layer.matrix .= _changevalue(change, elapsed)
    return nothing
end

function _primechange!(change::AbstractLayerChange{RelativeChange},
                       layer::AbstractLayer, elapsed::Unitful.Time)
    layer.matrix .= change.baseline .+ _changevalue(change, elapsed)
    return nothing
end

# A rate is deliberately **not** primed, and this is the one case where priming would be wrong
# rather than merely unnecessary. A rate accumulates `value × timestep`; at the start of a run no
# time has passed, so nothing has accumulated. Writing anything here would add a phantom step's worth
# of drift before the first step - and there is no timestep at build to multiply by in any case.
_primechange!(::AbstractLayerChange{RateChange}, ::AbstractLayer, _) = nothing

_checkreach(::AbstractLayer, _) = nothing

function _checkreach(layer::AbstractLayer{Condition}, final::Unitful.Time)
    lo, hi = bounds(axisof(layer))
    (isnothing(lo) && isnothing(hi)) && return nothing
    reach = _reach(layer.change, layer, final)
    isnothing(reach) && return nothing
    _refusereach(layer, lo, reach, <, "below")
    _refusereach(layer, hi, reach, >, "above")
    return nothing
end

# Refuse a change whose values would take a layer past a declared bound. The `nothing` method is the
# no-bound case, which is not an error: an axis that declares no limit has nothing to exceed.
_refusereach(::AbstractLayer, ::Nothing, ::NamedTuple, _, _) = nothing

function _refusereach(layer::AbstractLayer, bound, reach::NamedTuple, outside,
                      which)
    limit = _inlayerunit(bound, layer)
    worst = outside(1, 0) ? reach.hi : reach.lo
    outside(worst, limit) || return nothing
    return error("this change would drive a $(nameof(axisof(layer))) condition $which its physical " *
                 "range before the run ends: the limit is $limit, and over the run the layer " *
                 "reaches $worst. A condition outside its range is impossible rather than merely " *
                 "exhausted - unlike a supply, which may legitimately run out - so the scenario is " *
                 "asking for something that cannot happen. Shorten the run, reduce the rate, or " *
                 "bound the change yourself. (Reported before the first timestep rather than at the " *
                 "step it would happen.)")
end

# The lowest and highest values a layer can hold at any point in `[0, final]`, or `nothing` where
# that cannot be computed exactly. A named tuple so the two ends cannot be swapped at the call site.
_reach(::AbstractLayerChange, layer::AbstractLayer, _) = nothing

function _reach(::NoLayerChange, layer::AbstractLayer, _)
    return (lo = minimum(layer.matrix), hi = maximum(layer.matrix))
end

# A constant rate is monotone, so the extremes are at the two ends of the run and nowhere between.
function _reach(change::SteadyLayerChange, layer::AbstractLayer,
                final::Unitful.Time)
    drift = change.value * final
    ends = (minimum(layer.matrix), maximum(layer.matrix),
            minimum(layer.matrix) + drift, maximum(layer.matrix) + drift)
    return (lo = minimum(ends), hi = maximum(ends))
end

# An absolute series takes its stored values and nothing else, so its reach is exactly their extremes
# - plus the layer's own, which stand before the first slice and again under `RevertToLayer`.
function _reach(change::SeriesLayerChange{AbsoluteChange}, layer::AbstractLayer,
                _)
    return (lo = min(minimum(change.slices), minimum(layer.matrix)),
            hi = max(maximum(change.slices), maximum(layer.matrix)))
end

# ---------------------------------------------------------------------------
# Bounds: a supply cannot be negative
# ---------------------------------------------------------------------------
# Not a policy bolted onto `Resource` - it restates what makes something a resource. A resource is
# *rival* and consumed against a demand, so a negative amount of it has no meaning. The floor
# therefore belongs to the layer and is enforced automatically, exactly as `NaN` already is
# (`_zerogaps`), and no change *operation* needs to exist to express it.
#
# The response deliberately differs by *when* the negative appears:
# * at **build** - a supply layer's own values, or the slices an absolute replacement will write -
# it is the user's data or spec being wrong, and they can fix it, so failing fast is a service;
# * at **run time** it is emergent (an increment has run the supply down past zero), and aborting a
# long simulation would be hostile when "you cannot have less than none of a consumable" is the
# right reading. So it warns once and clamps.
#
# Conditions are *not* covered here. Their bounds are axis-level rather than role-level and
# genuinely discriminate - `ClimateMoisture` and `SiteWaterBalance` are balances and `Altitude` goes
# below sea level - and both the bound shape (floor only, or a pair?) and the runtime policy (a
# temperature below absolute zero is impossible rather than exhausted) are still open.

# The values a change is holding to write on some later timestep. Only *stored* values can be checked
# ahead of time; a change that computes its values has none to inspect.
_storedvalues(::AbstractLayerChange) = ()

_storedvalues(change::SeriesLayerChange{AbsoluteChange}) = (change.slices,)

function _storedvalues(change::SumOfLayerChanges{AbsoluteChange})
    return mapreduce(_storedvalues, (a, b) -> (a..., b...), change.parts,
                     init = ())
end

# Refuse a change that would write negative values into a supply. Only an *absolute* replacement
# can be judged this way: its stored slices are the values the layer will take, whereas an offset or
# a rate is combined with whatever the layer holds at the time and cannot be checked in isolation.
_checkstoredbounds(::AbstractLayerChange, ::AbstractLayer) = nothing

function _checkstoredbounds(change::AbstractLayerChange,
                            layer::AbstractLayer{Resource})
    for values in _storedvalues(change)
        neg = count(<(zero(eltype(values))), values)
        iszero(neg) ||
            error("this change would set a supply negative: $neg of its $(length(values)) stored " *
                  "values are below zero (the smallest is $(minimum(values))). A supply is a " *
                  "*resource* - rival, and consumed against a demand - so a negative amount of it " *
                  "has no meaning. Fix the data or the spec: if the quantity genuinely takes both " *
                  "signs (a net flux, a balance) it is not a supply, and belongs on the regime side.")
    end
    return nothing
end

# Every supply layer reachable from a habitat's `supply`, checked at construction.
function _checksupplybounds(supply::AbstractLayer{Resource})
    neg = count(<(zero(eltype(supply.matrix))), supply.matrix)
    iszero(neg) ||
        error("this supply holds $neg negative value(s) (the smallest is " *
              "$(minimum(supply.matrix))). A supply is a *resource* - rival, and consumed against a " *
              "demand - so a negative amount of it has no meaning. Fix the data or the spec: if the " *
              "quantity genuinely takes both signs (a net flux, a balance) it is not a supply, and " *
              "belongs on the regime side.")
    return supply
end

function _checksupplybounds(supply::LayerCollection{Resource})
    foreach(_checksupplybounds, values(supply))
    return supply
end

# Hold a supply at its floor after a change has written to it. `count`/`max` over the matrix once
# per layer per timestep, not per cell per species - this is not the hot loop.
# `maxlog = 1` is the "warn once" of the agreed policy: one warning per run, not per step or per
# layer, because a supply held at zero warns on *every* subsequent step and would otherwise drown the
# output of a long run.
_enforcebounds!(::AbstractLayer) = nothing

# A **condition** out of bounds is an error, not a clamp, and the asymmetry with supplies is the
# point rather than an inconsistency. A supply hitting zero is *expected* - resources get consumed,
# and running out is a legitimate ecological outcome. A temperature below absolute zero is not an
# exhausted quantity but an **impossible** one: it says the scenario is wrong. Clamping it would let
# a physically nonsensical run continue and produce output that looks like a result.
# Most such failures are caught *before* the run by `check_bounds` - this is the backstop for the
# changes whose reach cannot be predicted.
function _enforcebounds!(layer::AbstractLayer{Condition})
    lo, hi = bounds(axisof(layer))
    _refuseoutofbounds(layer, lo, <, "below")
    _refuseoutofbounds(layer, hi, >, "above")
    return nothing
end

function _enforcebounds!(layer::AbstractLayer{Resource})
    floor = zero(eltype(layer.matrix))
    neg = count(<(floor), layer.matrix)
    iszero(neg) && return nothing
    @warn "A change has driven a $(nameof(axisof(layer))) supply below zero in $neg cell(s) " *
          "(the smallest is $(minimum(layer.matrix))); clamping to $floor. A supply is a resource, " *
          "so it cannot go negative - but running out of one is a legitimate outcome, which is why " *
          "this is a warning rather than an error. Reported once per run." maxlog=1
    layer.matrix .= max.(layer.matrix, floor)
    return nothing
end

# The bound is converted into the *layer's* unit rather than the values into canonical, so the
# comparison is one conversion instead of one per cell - and so an **affine** axis is handled by
# Unitful rather than by hand: absolute zero on a `°C` layer converts to `-273.15 °C`, which is what
# the values must actually be compared against.
_refuseoutofbounds(::AbstractLayer, ::Nothing, _, _) = nothing

function _refuseoutofbounds(layer::AbstractLayer, bound, outside, which)
    limit = _inlayerunit(bound, layer)
    bad = count(v -> outside(v, limit), layer.matrix)
    iszero(bad) && return nothing
    worst = outside(1, 0) ? maximum(layer.matrix) : minimum(layer.matrix)
    return error("a change has driven a $(nameof(axisof(layer))) condition $which its physical " *
                 "range in $bad cell(s): the limit is $limit and the worst value is $worst. Unlike " *
                 "a supply - which may legitimately run out, and is clamped - a condition outside " *
                 "its range is impossible rather than exhausted, so this is an error: the scenario " *
                 "asks for something that cannot happen. Shorten the run, reduce the rate, or bound " *
                 "the change yourself.")
end

# A dimensionless bound stays as it is; a dimensioned one is converted.
_inlayerunit(bound::Number, ::AbstractLayer) = bound

function _inlayerunit(bound::Unitful.Quantity, layer::AbstractLayer)
    return uconvert(unit(eltype(layer.matrix)), bound)
end

# Whether a change is holding values, to be written on some later timestep, that contain data gaps.
# Only *stored* values can carry a layer's gaps forward past a clean-up of its matrix; a change that
# computes its values has none to carry.
_hasgaps(::AbstractLayerChange) = false

_hasgaps(change::SeriesLayerChange) = any(isnan, change.slices)

_hasgaps(change::SumOfLayerChanges) = any(_hasgaps, change.parts)

# Rebuild `change` with its stored values passed through `clean!`. Cleaning a layer's matrix alone is
# not enough once a series holds the slices that will *become* that matrix: the first update would
# write the gaps straight back in, one step late and a long way from whatever removed them.
_cleanstored(change::AbstractLayerChange, clean!) = change

function _cleanstored(change::SeriesLayerChange{M, E, C},
                      clean!) where {M, E, C}
    slices = clean!(copy(change.slices))
    return SeriesLayerChange{M, E, C, typeof(slices),
                             typeof(change.baseline)}(slices,
                                                      change.times,
                                                      change.origin,
                                                      change.atend,
                                                      change.calendar,
                                                      change.baseline)
end

function _cleanstored(change::SumOfLayerChanges{M}, clean!) where {M}
    parts = map(p -> _cleanstored(p, clean!), change.parts)
    return SumOfLayerChanges{M, typeof(parts), typeof(change.baseline)}(parts,
                                                                        change.baseline)
end

# ---------------------------------------------------------------------------
# Declaring a change alongside a layer, at the `GridHabitat` boundary
# ---------------------------------------------------------------------------
# A change cannot be attached to a spec: it is validated against the layer's unit, and a relative
# change captures the layer's values, neither of which exists until the spec has been materialised
# on a grid. So a declaration is *carried* next to the spec by `Varying` and applied once the layer
# is built. Specs themselves stay naked - nothing upstream of `GridHabitat` knows that time
# variation exists, which is why the same spec can be handed to `StudyArea` (which cares only about
# the grid) and re-used wrapped at build time without the two drifting apart.

# Strip a `Varying`, giving the naked spec. Everything else passes through, so callers that do not
# care about change (the whole of `StudyArea`'s grid analysis) need no other guard.
_unwrapspec(x) = x

_unwrapspec(v::Varying) = v.spec

# The change a value declares, or `nothing` for a naked spec.
_changeof(::Any) = nothing

_changeof(v::Varying) = v.change

# Split a `regime`/`supply` argument into the naked spec(s) and the declared change(s). For a tuple
# both halves are tuples, aligned positionally, so a second element and its change stay together.
_splitvarying(x) = (spec = _unwrapspec(x), change = _changeof(x))

# A named tuple splits the same way and keeps its names on both halves - `NamedTuple <: Tuple` is
# `false`, so it needs saying explicitly or it would take the scalar path above.
function _splitvarying(t::Union{Tuple, NamedTuple})
    return (spec = map(_unwrapspec, t), change = map(_changeof, t))
end

# Apply the changes declared alongside a regime/supply to the layers actually built from them.
# Positional throughout: element `i` of a tuple regime becomes sub-layer `i` of the collection.
_applydeclared!(layer::AbstractLayer, ::Nothing) = layer

function _applydeclared!(coll::LayerCollection,
                         changes::Union{Tuple, NamedTuple})
    parts = values(coll)
    _samearity(parts, values(changes))
    foreach(_applydeclared!, parts, values(changes))
    return coll
end

function _applydeclared!(layer::AbstractLayer, change)
    # A layer read from a monthly stack already carries a `SeriesLayerChange`, which is what makes it vary.
    # A change declared on top of that is *added* to it rather than replacing it - a seasonal
    # pattern and a multi-year trend offset one another, they do not compete - which is why
    # composition had to exist before this could be anything but an error.
    layer.change isa NoLayerChange && return setchange!(layer, change)
    layer.change = _combineparts((layer.change,
                                  _attachchange(change, layer)), layer)
    return layer
end

# ---------------------------------------------------------------------------
# The simulation epoch
# ---------------------------------------------------------------------------
# A run's epoch is the real date its elapsed time zero corresponds to. It is resolved once, by
# `build_ecosystem`, from the series the environment already carries - never copied into a parallel
# structure on the habitat, which would leave two places to say the same thing. Once resolved it is
# applied *backwards*, by rewriting each series' `origin`: one number per series, after which the
# hot path indexes by elapsed time exactly as before and does no date arithmetic at all.
#
# It has to work this way round because the ecosystem does not exist when the layers are built.
# Series times are resolved during `GridHabitat`, long before `build_ecosystem` runs, so the
# epoch cannot simply be consulted where a series is created.

# The start dates of every dated series in a habitat, deduplicated and in a deterministic order.
# Determinism is load-bearing under MPI: every rank builds the layers redundantly and must land
# on the same epoch, so this walks the layer structure rather than anything rank-local or hashed.
# `habitat` is duck-typed on `.regime`/`.supply`, so it takes either a real habitat or the
# `(regime = ..., supply = ...)` pair the layer machinery passes around before one exists.
function _startdates(habitat::Union{AbstractHabitat, NamedTuple})
    return unique(filter(!isnothing,
                         map(_startdate,
                             vcat(_layercalendars(habitat.regime),
                                  _layercalendars(habitat.supply)))))
end

# Only a dated series knows a real date; the other calendars have none to contribute.
_startdate(::AbstractSeriesCalendar) = nothing

_startdate(calendar::DatedSeries) = calendar.start

# Every series calendar reachable from a layer, walking collections and summed changes alike.
function _layercalendars(layer::AbstractLayer)
    return _changecalendars(layer.change)
end

function _layercalendars(layer::LayerCollection)
    return mapreduce(_layercalendars, vcat, values(layer),
                     init = AbstractSeriesCalendar[])
end

# Every calendar a change carries, as a vector so several changes can be concatenated. Only a series
# has one; the empty default is what lets a caller ask any change without branching on its type.
_changecalendars(::AbstractLayerChange) = AbstractSeriesCalendar[]

function _changecalendars(change::SeriesLayerChange)
    return AbstractSeriesCalendar[change.calendar]
end

function _changecalendars(change::SumOfLayerChanges)
    return mapreduce(_changecalendars, vcat, change.parts,
                     init = AbstractSeriesCalendar[])
end

# Resolve the run's epoch: an explicit one always wins, otherwise adopt the habitat's if it is
# unambiguous and refuse to guess if it is not. The same "adopt if unambiguous, ask if not" rule
# `StudyArea` uses to settle a CRS, so the behaviour is one users have already met.
_resolveepoch(::Any, given::Dates.TimeType) = given

function _resolveepoch(habitat, ::Nothing)
    dates = _startdates(habitat)
    isempty(dates) && return nothing
    length(dates) == 1 && return only(dates)
    return error("this environment's series begin on $(length(dates)) different dates " *
                 "($(join(dates, ", "))), so there is no unambiguous date for elapsed time zero " *
                 "to mean. Pass `epoch = ` to `build_ecosystem` to say which date the run starts " *
                 "from; every series is then placed relative to it.")
end

# Re-point every series in a habitat against the resolved epoch. Rebuilds each change and reinstalls
# it: `SeriesLayerChange` is immutable but a layer is not, which is the same shape `_cleanstored`
# uses to rewrite a change in place.
function _repointseries!(layer::AbstractLayer, epoch)
    layer.change = _repointchange(layer.change, epoch)
    return layer
end

function _repointseries!(layer::LayerCollection, epoch)
    foreach(l -> _repointseries!(l, epoch), values(layer))
    return layer
end

# Re-express a change against the run's epoch. Only a series moves - its origin is computed from its
# calendar and the epoch - so everything else is returned untouched, which is what makes re-pointing
# idempotent and safe to apply to a whole environment.
_repointchange(change::AbstractLayerChange, epoch) = change

function _repointchange(change::SeriesLayerChange{M, E, C},
                        epoch) where {M, E, C}
    origin = _calendarorigin(change.calendar, change, epoch)
    origin == change.origin && return change
    return SeriesLayerChange{M, E, C, typeof(change.slices),
                             typeof(change.baseline)}(change.slices,
                                                      change.times, origin,
                                                      change.atend,
                                                      change.calendar,
                                                      change.baseline)
end

function _repointchange(change::SumOfLayerChanges{M}, epoch) where {M}
    parts = map(p -> _repointchange(p, epoch), change.parts)
    return SumOfLayerChanges{M, typeof(parts), typeof(change.baseline)}(parts,
                                                                        change.baseline)
end

# Where elapsed time zero falls in a series' own coordinate, given the run's epoch - the entirety of
# what an epoch does to a series.
# With no epoch, or for a series with no calendar identity for one to bind to, the series keeps the
# origin it was built with. An `UndatedSeries` is exactly the case `origin` is the knob for.
function _calendarorigin(::AbstractSeriesCalendar, change::SeriesLayerChange,
                         ::Nothing)
    return change.origin
end

function _calendarorigin(::UndatedSeries, change::SeriesLayerChange,
                         ::Dates.TimeType)
    return change.origin
end

# A dated series' coordinates are elapsed time from its own first slice, so the epoch's offset from
# that same date *is* the coordinate elapsed zero sits at.
# An epoch *before* the series begins is deliberately not an error. Nothing is invalid: the layer
# has perfectly good values of its own, and it is only the series that has nothing to say yet - which
# is precisely the case the layer stands for. "The record starts in 1970, run from 1950 on the spec's
# own climatology" is a thing to want, not a mistake. (It was an error while clamping to the first
# slice was the only alternative, which was silently wrong rather than merely unhelpful.)
function _calendarorigin(calendar::DatedSeries, change::SeriesLayerChange,
                         epoch::Dates.TimeType)
    return uconvert(s,
                    Dates.value(Dates.Millisecond(epoch - calendar.start)) *
                    Unitful.ms)
end

# A climatology is placed by *calendar month*, not by elapsed duration into the year, and the
# difference is a whole slice near the boundaries: the coordinates are `n * month_mean_duration`
# (30.44 d), so six of them reach 182.6 days while the real 1 July is day 181 - matching
# proportionally would start a July run on the June slice. Nothing checks whether a slice sits
# exactly at that month: `_seriesindex` takes the last slice at or before it and `atend` decides the
# rest, which is what makes a partial climatology (months 2:4, say) behave like any other series.
function _calendarorigin(::MonthOfYearSeries, change::SeriesLayerChange,
                         epoch::Dates.TimeType)
    return uconvert(s, float(Dates.month(epoch) * month_mean_duration))
end

# ---------------------------------------------------------------------------
# Layer change
# ---------------------------------------------------------------------------
# The per-layer change rule - what a layer does to itself each timestep. Drives any layer, regime or
# supply. The types live here (a layer has a `change` field, so they must precede it); the unit
# contract and the apply methods are in `LayerChange.jl`, included later.
#
# A change carried as a **function reference** could be neither dispatched on nor validated, and each
# such function would hard-code its own unit - `K` for temperature, `mm` for rainfall - which is
# exactly how a change's unit and its layer's unit drift apart. See `changeunit` in `LayerChange.jl`
# for the contract that replaces it.

# ---------------------------------------------------------------------------
# AbstractLayer - the materialised, hot-loop grid-layer family
# ---------------------------------------------------------------------------
# A materialised layer is `matrix + size + change`, tagged with a `Role` - `Condition` for a regime,
# `Resource` for a supply - and a `NicheAxis` saying what it measures. One family covers both sides
# of the environment, so a regime and a supply differ in their role rather than in their machinery.

# **`NicheAxis`, the root, is the default axis.** As the *ancestor* of every real axis it says "I
# could be anything", where a sibling of them would suggest it meant something specific. A case that
# genuinely means "I am something you do not represent" is a new axis declared with `@nicheaxis`.
#
# **The distinction it existed to draw survives, and is still load-bearing**: `nothing` means *no
# axis was named*, which is a different statement from *this axis is dimensionless* (`NoUnits`). While
# both were spelled `NoUnits` the two could not be told apart, and the dimensionless branch of
# `_tocanon` could not be made strict - `285.0K` on an axis-less layer is the **common** case and had
# to keep working, so `285.0K` on a genuinely dimensionless axis (`Isothermality`) had to be accepted
# too. An axis-less layer takes any unit (there is no axis to disagree with it); a
# declared-dimensionless one refuses.
# Two deliberate consequences elsewhere, now keyed on the root rather than on a sibling:
# `_checksupport` needs its own exact-`NicheAxis` method, or an axis-less tolerance is told it
# "declares `condition = nothing`" - true of the method table, wrong as advice; and `layerrate` stops
# dividing an axis-less layer by its period. The second reaches no shipped data: every row in all
# six catalogues names a real axis, none is blank.

# The role a layer is on, for `_sharedrole` (`collections.jl`), which both role-parameterised
# families - layers and species requirements - build their collections through.
_roleof(::AbstractLayer{R}) where {R} = R

# A layer's axis structure, straight off the type parameter - an axis for a leaf, a `Tuple` of
# them for a collection. See `AbstractLayer`'s docstring for why the two share one slot.
axisof(::AbstractLayer{R, A}) where {R, A} = A

# Named access to a collection's layers. Everything is forwarded to the backing, including `:nt`
# itself, so that a layer *named* `:nt` stays reachable; internal code uses `getfield(lc, :nt)`.
#
# **The bulk accessors are gone** - `regimes`/`supplies`, `regime_names`/`supply_names` and
# `named_regimes`/`named_supplies` were fifteen names across five families saying one thing, and a
# collection implements the **standard container interface** instead (`src/collections.jl`).
# Write `values(lc)`, `keys(lc)` and `NamedTuple(lc)`, which work identically on a bare layer -
# a leaf being a one-member container.

# The two layer roles as sides of a pairing check (see `_checkaligned`, `collections.jl`). A supply
# has no `kinds`: continuous-versus-categorical is a Condition distinction, and there is no
# categorical resource for it to tell apart.
_regimeside(r::AbstractLayer{Condition}) = _side(r, "environment regime", true)

_supplyside(s::AbstractLayer{Resource}) = _side(s, "environment supply", false)

# Wrap a plain 2-D array as a `(Y, X)` `DimArray` carrying **real coordinates** derived from the cell
# size: origin `(0, 0)`, step `size`, `Intervals(Start)`.
#
# **This is `_syntheticyx`'s construction, applied to a layer instead of a study area** - and it is
# what lets a layer's cell size be *derived from its own dims* rather than stored beside them
# (`[CELL-DERIVE]`). A grid with no coordinates has nothing for `getcellsizes` to read, which is the
# whole reason the `.size` field had to be supplied and could therefore lie.
#
# **Real coordinates, not a `NoLookup` wrapper.** The reproducibility argument for `NoLookup` - that
# two independently built pairs of the same size compare equal without a shared lookup object being
# threaded through - holds for a *data-driven* lookup and not for a **derived** one:
# `(0:(n - 1)) .* size` is reproducible from `(n, size)` alone, so two layers built independently
# from the same shape and cell size still compare equal.
#
# Only the *origin* is a convention; the spacing is real and follows from `size`.
function _sizedyx(M::AbstractMatrix, size)
    return DimArray(M,
                    map((n,
                         D) -> D(DimensionalData.Lookups.Sampled((0:(n - 1)) .*
                                                                 size,
                                                                 order = DimensionalData.Lookups.ForwardOrdered(),
                                                                 span = DimensionalData.Lookups.Regular(size),
                                                                 sampling = DimensionalData.Lookups.Intervals(DimensionalData.Lookups.Start()))),
                        Base.size(M), (Y, X)))
end

# The `(Y, X)` dims of a layer's own array - for a collection, every sub-layer must already
# agree (checked here, not assumed elsewhere) and the one shared value is returned. This is the
# single source of truth threaded through `GridHabitat`'s construction (active/supply reuse it
# rather than each independently rebuilding their own).
_yx(layer::AbstractLayer) = dims(layer.matrix, (Y, X))

function _yx(layer::LayerCollection)
    ls = values(layer)
    yx = _yx(first(ls))
    all(l -> _yx(l) == yx, Base.tail(ls)) ||
        error("Regime collection's sub-layers are on different grids.")
    return yx
end

# Whether two `(Y, X)` dims pairs may be combined into one `GridHabitat`. Always the same size, which
# is what catches an X/Y-order mismatch on a non-square grid; and where **both** sides carry real,
# non-`NoLookup` provenance, they must also match in value, which catches two independently sourced
# real-CRS arrays that do not in fact share a grid.
#
# A supply or mask deliberately built with no provenance of its own - a pre-built supply from a bare
# array - is a legitimate pattern rather than a bug, so `NoLookup` on either side falls back to the
# size check alone rather than rejecting the combination.
# **A layer must be able to say where its cells are.** A `NoLookup` grid carries pure indices, so
# `getcellsizes` cannot answer for it, `_uniformcellside` cannot set up dispersal, and there is
# nothing to derive `size` from - the field would have to be *supplied*, which is the second source
# of truth this rule exists to remove.
#
# **Refused at construction rather than tolerated.** Nothing in the package builds such a layer, since
# every `Matrix` constructor derives real coordinates, so this can only fire on one assembled by
# hand - and it is better to hear about that than to carry it. **The one sanctioned exception is
# `deprecations.jl`**, which reproduces the released coordinate-less budgets and reaches the *inner*
# constructor directly, bypassing this check on purpose.
# A layer's cell size, read from its own coordinates rather than supplied beside them.
#
# **This is what makes `size` a derived fact rather than a second source of truth.** Passed in and
# stored unchecked, it can disagree with the grid the values are actually on and nothing will say so.
# `_checkhascoords` runs first everywhere this is called, so the dims are always real.
# Reads the **declared** `Regular` span via `_axisstep`, not a difference of the first two
# coordinates - see the note there on why differencing drifts in the 13th digit.
function _derivecellsize(matrix)
    yx = dims(matrix, (Y, X))
    return _axisstep(yx[1])
end

# A layer must be able to say where its cells are, so `NoLookup` dims are refused at construction
# rather than failing later inside whatever needed the spacing.
function _checkhascoords(matrix)
    yx = dims(matrix, (Y, X))
    all(_isreallookup, yx) ||
        return error("a layer needs real `(Y, X)` coordinates, but this one was given `NoLookup` " *
                     "dims - bare cell indices, which cannot say how far apart the cells are. " *
                     "Build it on a `StudyArea`, or pass a matrix and a cell size and let the " *
                     "coordinates be derived.")
    # **And they must carry a unit.** `Y(1:3)` is a real lookup but a unitless one, so the cell
    # size derived from it would be a bare `1` - which the layer's `S <: Unitful.Quantity` then
    # rejects with a `MethodError` naming five type parameters. Caught here instead, where the
    # remedy can be stated. Consistent with what the check is for: a layer must say *where* its
    # cells are, and indices dressed as numbers do not.
    all(d -> eltype(parent(DimensionalData.lookup(d))) <: Unitful.Quantity,
        yx) ||
        return error("a layer's `(Y, X)` coordinates must carry a unit - these are bare numbers, " *
                     "so how far apart the cells are is unstated. Give the lookups a unit " *
                     "(`Y((0:3)km)`), or pass a matrix and a cell size.")
    return nothing
end

# Whether a dimension carries real coordinates rather than bare indices dressed as numbers.
function _isreallookup(d)
    return !(DimensionalData.lookup(d) isa DimensionalData.Lookups.NoLookup)
end

# Whether two `(Y, X)` dimension pairs describe the same grid. Shapes must match always; coordinate
# *values* are compared only when both sides carry real lookups, since a synthetic layer has none to
# compare and would otherwise be refused for saying nothing.
function _yxcompatible(yx1, yx2)
    length.(yx1) == length.(yx2) || return false
    (_isreallookup(yx1[1]) && _isreallookup(yx2[1])) || return true
    return yx1 == yx2
end

# ---------------------------------------------------------------------------
# Back-compat aliases (supply role)
# ---------------------------------------------------------------------------
# A supply is a `Resource`-role layer: the old supply structs are aliases over
# `ContinuousLayer{Resource, axis, V, Arr}`, all over `(Y, X)`. The axis records the
# resource measured; the (unused) `size` and the `change` rule are filled by the
# constructors in Energy.jl. Every supply names its axis - the free/axis-less family is gone. A supply that
# varies in time is not a separate type - it is one of these carrying a [`SeriesLayerChange`](@ref) change.

# ---------------------------------------------------------------------------
# Axis-driven canonicalisation + axis re-tagging (shared by the hand `*AE`
# constructors in `GridHabitat.jl` and the `build_*`/`materialise` path in `materialise.jl`)
# ---------------------------------------------------------------------------

# Canonicalise a regime *value* (a position) to its layer axis's unit, `canonicalunit(A))`: a unitful
# value converts (proper affine - canonical units are absolute, so no interval subtlety); a bare value
# attaches the unit. An axis with no canonical unit (`NicheAxis`/dimensionless, `NoUnits`/`nothing`)
# keeps the value's magnitude but still **absolutises** an affine unit (°C->K) - regimes are always in an
# absolute unit, which the downstream change machinery assumes. The single, axis-driven replacement
# for the old dimension-sniffing conversions.
function _canonical(x, axis::Type{<:NicheAxis})
    return _tocanon(x, canonicalunit(axis), axis)
end

# Three cases, as three methods rather than one branch, because they are three different
# statements that happen to share an implementation today - and keeping them apart is what lets any
# of them change without disturbing the others.
#
# 1. `nothing` - the axis has no condition unit: either nobody declared one, or it declared
# `condition = nothing` because it is a resource, not a condition. Nothing to check against, so
# the value passes. **This must stay permissive**: reference layers are exactly the layers with
# no canonical unit, and making them unbuildable is a proposal already withdrawn twice.
_tocanon(x, ::Nothing, axis) = _absolutise(float(x))

# 2. The axis declares itself **dimensionless** - meaning a bare **fraction**, so a value that is
# dimensionless but *scaled* (a percentage) is converted rather than passed through.
# **Why convert here and not leave `%` on the values.** `percent` and `NoUnits` share the
# dimension `NoDims` and differ only by a factor of 100, so `uconvert(NoUnits, 64.85%)` is exactly
# `0.6485` - a plain `Float64`, not a `Quantity`. Doing it once, here, is what stops a `%` unit
# riding along through every downstream multiplication and needing to be stripped at each use.
# A **dimensioned** value still passes through, and that hole stays open: `NicheAxis`, the
# default axis for any layer built without one, also declares `NoUnits`, and unit-bearing
# unclassified layers are ordinary (a `285.0K` matrix with no axis is the common case). So
# `NoUnits` still means both "genuinely dimensionless" and "no axis given", and only the first is
# tightened here. Separating them needs `canonicalunit(NicheAxis)` to become `nothing`,
# which is decided but not yet built.
_tocanon(x, ::typeof(NoUnits), axis) = _tocanondimensionless(float(x), axis)

# 3. A real unit: convert, which is also what rejects a dimension mismatch.
_tocanon(x, u, axis) = _tocanonu(float(x), u)

_tocanondimensionless(x::Real, axis) = x

function _tocanondimensionless(x::Unitful.AbstractQuantity, axis)
    dimension(x) == NoDims && return uconvert(NoUnits, x)
    return error("a regime value on axis $(nameof(axis)) must be dimensionless: it declares " *
                 "`condition = NoUnits`, so its values are bare fractions, and $x carries " *
                 "$(dimension(x)). Either strip the unit, or name the axis the value really " *
                 "belongs to.")
end

_tocanonu(x::Unitful.AbstractQuantity, u) = uconvert(u, x)

# A bare number carries no unit, so nothing can check it - stamping the canonical unit on would
# accept `12` as 12 K from someone who meant 12 °C, silently. Every other route into a layer is
# checked (`_checkaxisunit` covers every shipped catalogue row *and* every user `axis =` override,
# and the `AbstractQuantity` method above converts, erroring on a dimension mismatch), so this was
# the one remaining hole. Measured before closing it: **no test, canonical run, docs example or
# example script reaches this branch** - `core_test`, `extras_canonical`, `extras_docs` and
# `extras_examples` all passed unchanged with it erroring, so refusing costs nothing.
function _tocanonu(x::Real, u)
    return error("a regime value on an axis with canonical unit `$u` must carry a unit: got the " *
                 "bare number $x, which could be anything. Write `$(x)$u` if that is what you " *
                 "mean - a bare number cannot be checked against the axis.")
end

# The same quantity on its absolute scale, so an affine temperature becomes kelvin. Widths and bare
# numbers pass through: only a *position* on an affine scale needs converting, and treating an
# interval as one would shift it by the offset.
function _absolutise(x::Unitful.AbstractQuantity)
    return uconvert(Unitful.absoluteunit(unit(x)), x)
end

_absolutise(x::Real) = x

# The `Resource`-role mirror of `_canonical` above: a supply *value* - already a per-cell rate, so
# `cancel` has done the × area - converted to its axis's canonical resource unit,
# `canonicalunit(Resource, A)`. Lives here beside the Condition-role conversion because the two are
# one idea in two roles, and keeping them together is what stops either drifting.
#
# **This is where the unit guarantee lives.** `Supply{A}` leaves its value type free - deliberately,
# since `canonicalunit(Resource, A)` is the single statement of what the unit is - so nothing else
# checks it, and a per-resource type that pinned the value type would be a second statement of the
# same fact. Without this check a wrong **dimension** builds happily and then throws a
# `DimensionError` deep inside the threaded hot loop, nowhere near the call that made the mistake.
#
# Resource-hood is asked of `supplytype`, **not** of `canonicalunit(Resource, axis)`: the role form
# falls back to the *condition* unit for an axis that declares no resource, so
# `canonicalunit(Resource, WindSpeed)` is `m/s` rather than `nothing` and would report a dimension
# mismatch instead of the real mistake. `supplytype` is the authority the wind-speed case
# already keys on, and the macro's all-or-nothing rule keeps the two in step.
function _canonicalresource(x, axis::Type{<:NicheAxis})
    return uconvert(_resourceunit(typeof(x), axis), x)
end

# The array form resolves and checks the unit **once**, not once per cell: `_resourceunit` is the
# whole of the checking, so hoisting it out of the broadcast costs nothing in clarity.
function _canonicalresource(values::AbstractArray, axis::Type{<:NicheAxis})
    return uconvert.(_resourceunit(eltype(values), axis), values)
end

# The canonical resource unit for `axis`, checked against values of type `V` - the two failures a
# supply's values can have, told apart, since they are different mistakes with different fixes.
# `noun` exists so a **demand** gets the same two checks and says "demand" while doing it: both
# sides of a resource share one canonical unit (that is what makes them comparable at all), so they
# must share the check - but a message naming the wrong side sends the reader to the wrong call.
function _resourceunit(::Type{V}, axis::Type{<:NicheAxis},
                       noun::AbstractString = "supply") where {V}
    isnothing(supplytype(axis)) && _notaresource(axis, noun)
    u = canonicalunit(Resource, axis)
    dimension(V) == dimension(u) ||
        error("a $noun on axis $axis is measured in `$u`, but these values are " *
              "$(_unitdescription(V)). Either the values or the declared axis is wrong.")
    return u
end

# A bare number has no unit to name, and saying so beats printing an empty `NoUnits`.
_unitdescription(::Type{V}) where {V} = "`$(unit(V))` ($(dimension(V)))"

_unitdescription(::Type{<:Real}) = "bare numbers, which carry no unit at all"

# The niche axis a layer is tagged with. Only the concrete types carry it - `AbstractLayer{R}` is
# parameterised by role alone - so anything generic over layers has to ask for it here rather than
# match it in a signature. Used by the Plots recipe to name what it is drawing.
_layeraxis(::ContinuousLayer{R, A}) where {R, A} = A

_layeraxis(::CategoricalLayer{A}) where {A} = A

# Re-tag a materialised layer with its niche axis `A` - a phantom type parameter, so this shares the
# arrays. The low-level constructors build a root-axis layer that this narrows to the real axis.
# `S` is carried through explicitly rather than left free: re-tagging changes only the axis, so the
# cell-size type must survive it - and on a geographic layer that type is an angle, not a length.
function _reaxis(l::ContinuousLayer{Condition, A0, V, Arr, S},
                 ::Type{A}) where {A0, V, Arr, S, A <: NicheAxis}
    return ContinuousLayer{Condition, A, V, Arr, S}(l.matrix, l.size, l.change)
end

function _reaxis(l::CategoricalLayer{A0, V, Arr, S},
                 ::Type{A}) where {A0, V, Arr, S, A <: NicheAxis}
    return CategoricalLayer{A, V, Arr, S}(l.matrix, l.size, l.change)
end
