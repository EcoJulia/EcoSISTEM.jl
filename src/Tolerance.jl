# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a species can put up with: a built response distribution on a continuous axis, or a set of
# acceptable classes on a categorical one.

using Unitful

using Distributions

"""
    ContinuousTolerance{A, V <: Number} <: AbstractTolerance{A, V}

A species' preference along a **continuous** niche axis, where suitability falls off with distance
from an optimum rather than being a matter of membership. [`NicheTolerance`](@ref) is the concrete
type; the matching layer is a [`ContinuousRegime`](@ref).
"""
abstract type ContinuousTolerance{A, V <: Number} <: AbstractTolerance{A, V} end

"""
    AbstractCategoricalTolerance{A, V} <: AbstractTolerance{A, V}

A species' preference among **class labels** — the tolerances matched against a
[`CategoricalRegime`](@ref), where preference is membership rather than a distance.

Defined by what it **answers** rather than by what it stores: *the suitability weight species `sp`
gets in category `c`*. [`SimpleCategoricalTolerance`](@ref) answers that sparsely, from a set of
acceptable categories per species; a graded type could answer it densely, with a different weight for
every species in every category.

Categorical, so `iscontinuous` is `false`: these are class labels, and nothing between two of them is
meaningful. A species' response to *how much* of a class covers a cell is a continuous
[`NicheTolerance`](@ref) on [`SurfaceArea`](@ref) instead.
"""
abstract type AbstractCategoricalTolerance{A, V} <: AbstractTolerance{A, V} end

"""
    NicheTolerance{A <: NicheAxis, V, D} <: ContinuousTolerance{A, V}

One species' response along a continuous niche axis, as a **built** probability distribution: where
on the axis it does best, and how sharply that falls away.

The distributions are built once, at construction, in a **frame** — the unit their support is
measured in, which is also the tolerance's `eltype` and the unit the matching regime must be in. In
the hot loop the [`NicheSuitability`](@ref) nichefit fetches one and evaluates its density at a bare
number already in that frame, so there is no per-call conversion and no allocation.

# Fields

  - `dists`: one distribution per species. Reach them with [`getdist`](@ref).

# Type parameters

  - `A`: the niche axis.
  - `V`: the frame the distributions are built in.
  - `D`: the concrete distribution type, a `Distributions.ContinuousUnivariateDistribution` such as
    `Trapezoid{Float64}` or `Uniform{Float64}`.
"""
struct NicheTolerance{A <: NicheAxis, V, D} <: ContinuousTolerance{A, V}
    dists::Vector{D}
end

"""
    NicheTolerance(axis::Type{<:NicheAxis}, ::Type{D}, dist::Matrix; support = _defaultsupport(axis),
        offset = nothing, scale = nothing, probes = …)

Build a [`NicheTolerance`](@ref) from a **bare** parameter matrix, one species per row.

# Arguments

  - `axis`: the niche axis the tolerance is on.
  - `D`: the response distribution type to build.
  - `dist`: the parameters, one row per species, as bare numbers read in the `support` frame.
  - `support`: the **frame** the distributions are built in — the unit their support is measured in,
    the tolerance's `eltype`, and the unit the matching regime must be in. Defaults to the axis's
    canonical unit.
  - `offset`, `scale`: references in `support`, which place a *shape-only* distribution such as
    `Beta` or `LogNormal` onto the dimensioned axis through a `LocationScale`.
  - `probes`: how `D`'s parameter roles are determined; see `param_roles_resolved`.
"""
function NicheTolerance(axis::Type{A}, ::Type{D}, dist::Matrix;
                        support = _defaultsupport(axis), offset = nothing,
                        scale = nothing,
                        probes = _defaultprobes(D)) where {A <: NicheAxis, D}
    _checksupport(axis, support)
    return _buildniche(axis, D, dist, support, support, offset = offset,
                       scale = scale, probes = probes)
end

"""
    NicheTolerance(axis::Type{<:NicheAxis}, ::Type{D}, params::AbstractVector...; support = _defaultsupport(axis),
        offset = nothing, scale = nothing, probes = …)

Build a [`NicheTolerance`](@ref) from **one vector per parameter** of the distribution — the
programmatic counterpart of the matrix constructor:

```julia
NicheTolerance(Temperature, Normal, opts, vars)          # a mu and a sigma vector
NicheTolerance(Precipitation, Gamma, shape, scale_vec)
```

The vectors may carry units, and each is read according to its parameter's **role** — a location
converted as a position, a scale as an interval, a shape left dimensionless. So `support = K` and
`support = °C` build the *same* preference, expressed in each frame.

# Arguments

  - `axis`: the niche axis the tolerance is on.
  - `D`: the response distribution type to build.
  - `params`: one vector per parameter of `D`, each with one entry per species. All must share a
    unit, and all must be the same length.
  - `support`: the **frame** the distributions are built in, as for the matrix constructor.
  - `offset`, `scale`, `probes`: as for the matrix constructor.
"""
function NicheTolerance(axis::Type{A}, ::Type{D}, params::AbstractVector...;
                        support = _defaultsupport(axis), offset = nothing,
                        scale = nothing,
                        probes = _defaultprobes(D)) where {A <: NicheAxis, D}
    isempty(params) &&
        error("NicheTolerance needs at least one parameter vector for $D.")
    n = length(first(params))
    all(v -> length(v) == n, params) ||
        error("NicheTolerance parameter vectors must all have the same length (one entry per species).")
    _checksupport(axis, support)
    input_unit = _imputeinputunit(params, support)
    # A unitful vector becomes its magnitude in `input_unit` — its role, applied in `_buildniche`,
    # decides position against interval — while a bare vector passes straight through.
    cols = map(v -> unit(eltype(v)) === NoUnits ? float.(v) :
                    ustrip.(input_unit, v),
               params)
    return _buildniche(axis, D, reduce(hcat, cols), input_unit, support,
                       offset = offset, scale = scale, probes = probes)
end

"""
    TempTolerance

A temperature preference: a [`NicheTolerance`](@ref) on [`Temperature`](@ref), read in the kelvin
frame, with a `Trapezoid` response.

An **alias** rather than a type of its own, so it has the fields [`NicheTolerance`](@ref) has.
`TempTolerance(matrix)` is a shorthand for `NicheTolerance(Temperature, Trapezoid, matrix)`.
"""
const TempTolerance = NicheTolerance{Temperature,
                                     typeof(1.0 *
                                            canonicalunit(Temperature)),
                                     Trapezoid{Float64}}

# The alias shorthands, so that a temperature or rainfall tolerance can be built without naming its
# axis and response distribution.
TempTolerance(dist::Matrix) = NicheTolerance(Temperature, Trapezoid, dist)

"""
    RainTolerance

A rainfall preference: a [`NicheTolerance`](@ref) on [`Precipitation`](@ref), read in the `mm/day`
frame, with a `Uniform` response.

The frame is a **rate** rather than a depth: every shipped precipitation layer is read per unit time.

An **alias** rather than a type of its own, so it has the fields [`NicheTolerance`](@ref) has.
`RainTolerance(matrix)` is a shorthand for `NicheTolerance(Precipitation, Uniform, matrix)`.
"""
const RainTolerance = NicheTolerance{Precipitation,
                                     typeof(1.0 *
                                            canonicalunit(Precipitation)),
                                     Uniform{Float64}}

RainTolerance(dist::Matrix) = NicheTolerance(Precipitation, Uniform, dist)

"""
    SimpleCategoricalTolerance{A <: NicheAxis, V} <: AbstractCategoricalTolerance{A, V}
    SimpleCategoricalTolerance(vals; axis, penalty = 0.0)

The categories each species tolerates on niche axis `A` — one **set of acceptable
class values per species** — together with the weight a species gets **outside** its
set.

# Fields

  - `vals`: `vals[i]` lists the class values species `i` tolerates. A species with a
    single preferred class simply has a one-element set, so pass a plain vector of
    classes (one per species) and each is wrapped for you.
  - `penalty`: the suitability weight outside the set, in `[0, 1]`. Inside it the
    weight is always `1`.

**`penalty` is the whole soft-versus-hard distinction, and it is a number rather than a type.** `0.0`,
the default, is *hard* exclusion: the species cannot persist outside its classes at all, because a
zero suitability makes its death rate infinite. `0.5` is *soft* — it does worse there but survives.
Anything between is available and means what it says.

**`penalty` interacts with `params.survival`.** Suitability enters the demographics as
`suitability^±survival`, so at `survival = 0` every penalty is ignored and the tolerance does nothing
at all. That is a property of the model rather than of this type, and invisible from here, so it is
worth stating: a categorical tolerance needs a non-zero `survival` to have any effect.
"""
struct SimpleCategoricalTolerance{A <: NicheAxis, V} <:
       AbstractCategoricalTolerance{A, V}
    vals::Vector{Vector{V}}
    penalty::Float64

    # The sole constructor, so there is one spelling. `vals` is normalised per species by
    # `_categoryset`, which takes either a bare class or a collection of them, so the
    # one-preferred-class case and the set case are the same call.
    function SimpleCategoricalTolerance(vals::AbstractVector;
                                        axis::Type{A},
                                        penalty::Real = 0.0) where {A <:
                                                                    NicheAxis}
        0 ≤ penalty ≤ 1 ||
            error("a categorical tolerance's `penalty` is a suitability weight and must lie in " *
                  "[0, 1], but $penalty was given: 0 excludes a species from every category " *
                  "outside its set, 1 leaves it indifferent to them.")
        sets = map(_categoryset, vals)
        return new{A, eltype(eltype(sets))}(sets, Float64(penalty))
    end
end

# --- Reading a tolerance ----------------------------------------------------

"""
    getpref(tolerance::AbstractCategoricalTolerance, sp::Int64)

Return the class values species `sp` tolerates. Always a collection, with one element for a species
that has a single preferred class.
"""
function getpref(tolerance::AbstractCategoricalTolerance, sp::Int64)
    return _speciestolerance(tolerance, sp)
end

"""
    getpref(tolerance::NicheTolerance, sp::Int64)

Return the parameters of species `sp`'s response distribution, read off the distribution itself. Use
[`getdist`](@ref) to reach the distribution rather than its parameters.
"""
function getpref(tolerance::NicheTolerance, sp::Int64)
    return Distributions.params(_speciestolerance(tolerance, sp))
end

"""
    getdist(tolerance::NicheTolerance, sp::Int64)

Return the pre-built response distribution for species `sp`. This is what the hot loop reads: a plain
vector fetch, with no per-call construction and no allocation.
"""
function getdist(tolerance::NicheTolerance, sp::Int64)
    return tolerance.dists[sp]
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# The mirror of the layer displays in `Layer.jl`, member for member: same three slots, with
# **species where a layer has cells**. That symmetry is the point — a tolerance and the regime it is
# matched against should read alike, or the pairing the model rests on is invisible at the REPL.
#
# Without them the type parameter swamps the content: a 12-species two-axis tolerance collection
# prints over 2 000 characters, nearly all of it `Distributions.Normal{Float64}` repeated.
function Base.show(io::IO, t::NicheTolerance)
    return print(io,
                 "NicheTolerance($(nameof(axisof(t))), $(length(t.dists)) species, $(unit(eltype(t))))")
end

# `penalty` is shown because it is the whole of the soft/hard distinction, and it is a value rather
# than a type — so nothing else on the line reveals it.
function Base.show(io::IO, t::SimpleCategoricalTolerance)
    return print(io,
                 "SimpleCategoricalTolerance($(nameof(axisof(t))), $(length(t.vals)) species, ",
                 "penalty $(t.penalty))")
end

# ══ Functions ══════════════════════════════════════════════════════════════════════════════════

iscontinuous(trait::AbstractCategoricalTolerance) = false

iscontinuous(::NicheTolerance) = true

# One species' acceptable categories, from either spelling a caller may write: a bare class — the
# common single-preference case — or a collection of them. Dispatch rather than an `isa` branch,
# because the two spellings are distinguished by type alone.
_categoryset(v::AbstractVector) = collect(v)

_categoryset(v) = [v]

# Build a `NicheTolerance` whose distributions live in the `frame` unit — its `eltype`, and the unit
# the matching regime must be in — reading each row's bare parameters as being in `input_unit` and
# converting to `frame` **per role**: a location affinely, a scale as an interval, a shape left
# dimensionless. See `read_distribution`. The single place the storage frame is fixed.
function _buildniche(::Type{A}, ::Type{D}, mat::AbstractMatrix, input_unit,
                     frame;
                     offset, scale, probes) where {A <: NicheAxis, D}
    roles = param_roles_resolved(D, probes = probes)
    dists = [read_distribution(D, input_unit, mat[sp, :], canonical = frame,
                               offset = offset, scale = scale, roles = roles)
             for sp in Base.axes(mat, 1)]
    return NicheTolerance{A, typeof(1.0 * frame), eltype(dists)}(dists)
end

# The frame a tolerance is built in when the caller names none: `canonicalunit(axis)` for a declared
# axis, but **`NoUnits` for the bare root**, which declares `nothing`. Not because such a tolerance is
# dimensionless, but because no axis was named, so bare numbers mean "whatever the layer is in".
#
# **`::Type{NicheAxis}` exactly, not `::Type{<:NicheAxis}`** — this is the axis-less case, and the
# `<:` form would claim every axis in the package. Same rule as `_checksupport` below, and the one way
# either could go silently wrong while every same-axis test still passed.
_defaultsupport(axis::Type{<:NicheAxis}) = canonicalunit(axis)

_defaultsupport(::Type{NicheAxis}) = NoUnits

# The third state, the bare root, and a different statement from the other two: an axis-less tolerance
# is neither a missing declaration nor a refusal, but the ordinary case of building against a layer
# whose axis nobody named. Dispatch rather than a branch, because this is a static type; written as an
# `if` inside the general method it would sit alongside two runtime tests and read as a third
# exception.
#
# It deliberately **widens** what is accepted: any `support` is legal here, where the general method
# compares dimensions against the axis. Unit-bearing axis-less layers are ordinary — a `285.0K` matrix
# with no axis is the common case — and there is no axis to disagree with the frame, so there is
# nothing to check.
#
# **`::Type{NicheAxis}` exactly** — see `_defaultsupport` above. Widening this to `<:` would turn the
# support check off for *every* axis, which no test with a declared axis would notice.
_checksupport(::Type{NicheAxis}, support) = nothing

# Refuse a tolerance on an axis that has no canonical unit, naming the axis and which of the two
# reasons applies. Without it the failure is a bare `MethodError: no method matching
# dimension(::Nothing)`, which names neither the axis nor the tolerance and arrives from `Unitful`
# rather than from here.
#
# The two reasons need different remedies, and are told apart by **which method answers**: a
# deliberate `condition = nothing` — a supply-only axis such as `CarbonFlux` — has its own method,
# while an axis nobody has declared falls through to the root fallback. The first is a modelling
# mistake, matching species against something they consume; the second is a missing declaration.
function _checksupport(axis, support)
    cu = canonicalunit(axis)
    if isnothing(cu)
        stated = which(canonicalunit, Tuple{Type{axis}}) !==
                 which(canonicalunit, Tuple{Type{NicheAxis}})
        return error("cannot build a NicheTolerance on axis $axis: " *
                     (stated ?
                      "it declares `condition = nothing`, so it is not a condition at all — it is a resource species consume, not one they are matched against. Give the species a demand on it instead." :
                      "no canonical unit is declared for it, so there is no unit to build the tolerance in. Declare one with `@nicheaxis($axis <: …, condition = …)`."))
    end
    dimension(support) == dimension(cu) ||
        return error("NicheTolerance support unit $support and axis $axis's canonical unit " *
                     "$cu have different dimensions.")
    return nothing
end

# Impute the **input** unit of a set of parameter vectors — the unit their bare magnitudes are read
# in. That is their single shared dimensioned unit, ignoring dimensionless shape vectors, or the
# `support` frame where all of them are bare, so that bare vectors are read in the frame exactly as
# the matrix constructor reads them. Mixed units are an error.
function _imputeinputunit(params, support)
    us = unique(unit(eltype(v)) for v in params if unit(eltype(v)) !== NoUnits)
    isempty(us) && return support
    length(us) == 1 ||
        error("NicheTolerance parameter vectors carry differing units $(us); pass values in one unit.")
    return only(us)
end

# Several tolerances are a `SpeciesRequirementCollection` whose members happen to be tolerances,
# exactly as several regimes are a `LayerCollection{Condition}`; the role is read off the members. A
# collection implements the standard container interface (`src/collections.jl`), so write `values(x)`,
# `keys(x)` and `NamedTuple(x)`, which work identically on a leaf as a one-member container.

# The tolerance as a side of a pairing check — see `_checkaligned` (`collections.jl`). This is the
# **reference** side of the condition line, so it is what the regime and the nichefit are each
# compared against.
_toleranceside(t::AbstractTolerance) = _side(t, "species tolerance", true)

# Per-species niche optima — the distribution means — as a bare vector in the tolerance's own frame.
# What a caller reads when it needs the optima themselves rather than the distributions, as
# `continuous_evolve` does.
_nichemeans(tolerance::NicheTolerance) = Distributions.mean.(tolerance.dists)

# --- Building a tolerance from what a caller wrote ----------------------------

# A pre-built tolerance must already cover exactly `n` species — it was built elsewhere, so nothing
# here can fill it out, and a silent mismatch would pair species with the wrong niches.
function _checktolerancespecies(tolerance::AbstractTolerance, n::Integer)
    len = length(_nichedists(tolerance))
    len == n ||
        error("a pre-built tolerance covers $len species but `numspecies` is $n: it is used as " *
              "given, so it must already have one entry per species.")
    return nothing
end

# Whatever a tolerance stores per species — response distributions for a continuous one, class sets
# for a categorical one. Named for the shape rather than the field, so callers that only need "one
# entry per species" work on either without knowing which they hold.
_nichedists(t::NicheTolerance) = t.dists

_nichedists(t::SimpleCategoricalTolerance) = t.vals

# A single `(mean, width)` pair, filled out to `n` species as a Gaussian `NicheTolerance`.
function _nichetolerance(axis, tolerance, n::Integer)
    return NicheTolerance(axis, Normal,
                          _tofield(tolerance[1], n, "tolerance mean"),
                          _tofield(tolerance[2], n, "tolerance width"))
end
