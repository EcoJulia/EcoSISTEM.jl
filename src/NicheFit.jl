# SPDX-License-Identifier: LGPL-3.0-or-later
#
# How a tolerance is scored against the regime it is paired with, and how several such scores
# combine into one suitability.
#
# `Gauss`, `Trapeze` and `Unif` are also `AbstractNicheFit`s, but they are deprecated shims and
# stay in `deprecations.jl`, which this reorganisation does not touch.

using Unitful

using EcoSISTEM.Units

# **Fieldless and immutable, and both matter.** A nichefit carries no state - it is a callable marker
# built at the join - so a fieldless **immutable** is a singleton, `T() === T()`, and `isbits`, while
# a fieldless **mutable** is neither. Do not make one `mutable` without a field to justify it.
"""
    NicheSuitability{A, V} <: AbstractNicheFit{A, V}

The nichefit between a [`NicheTolerance`](@ref) continuous trait and its environment: the density of the
trait's response distribution evaluated at the current regime value, parameterised on any `V`.
Works for any `Distributions.ContinuousUnivariateDistribution` (e.g. [`Trapezoid`](@ref) or `Uniform`).
"""
struct NicheSuitability{A, V} <: AbstractNicheFit{A, V}
end

# == Functions ==================================================================================

# The `pdf` is a **density**, so it carries `1/x` and its stripped value depends on the frame.
# Multiplying by the axis's fixed physical `densitywidth`, expressed in that same frame, makes the
# result a dimensionless weight that is invariant to the axis's canonical unit - see
# `densitywidth`'s docstring. `_densityscale` is `1.0` for every shipped axis today (each width is
# one of its own current canonical units), so this changes no number and `reference.toml` must not
# move. An axis declaring no width is unscaled, exactly as before.
function (::NicheSuitability{A, V})(dist::ContinuousUnivariateDistribution,
                                    current) where {A, V}
    return pdf(dist, _toframe(V, current)) * _densityscale(A, V)
end

# == Functions ==================================================================================

# Convenience: the matching [`NicheSuitability`](@ref) nichefit for a `NicheTolerance`, taking its `V` (unit) from the
# trait's axis - so callers building an ecosystem by hand need not re-type the unit.
NicheSuitability(t::NicheTolerance) = NicheSuitability{axisof(t), eltype(t)}()

"""
    CategoricalSuitability{A, V} <: AbstractNicheFit{A, V}

The nichefit between a categorical trait on class labels of type `V` and its
environment: it reports the weight the species' tolerance gives the cell's class.

Carries no fields, like every other leaf nichefit - the whole response, both the
categories a species tolerates and the `penalty` it takes outside them, lives on the
species side in an [`AbstractCategoricalTolerance`](@ref). That is the same division
as the continuous branch, where [`NicheTolerance`](@ref) holds the response
distribution and [`NicheSuitability`](@ref) only evaluates it.

Soft and hard exclusion are therefore **not** different fits: they are different `penalty` values on
the tolerance.
"""
struct CategoricalSuitability{A, V} <: AbstractNicheFit{A, V}
end

# The deprecated `Gauss`/`Trapeze`/`Unif` nichefit shims (all `NicheSuitability` now) live in
# `src/deprecations.jl`.

# One argument, and deliberately: for a categorical fit the tolerance has already answered with
# the weight (`_categoryweight`), so there is nothing left for the fit to score. The continuous fits
# take two because a distribution genuinely has to be evaluated at the cell's value.
(::CategoricalSuitability)(weight::Real) = weight

# The two-argument form, for a caller invoking a nichefit by hand. It scores membership **hard** - 1
# inside, 0 outside - because the penalty belongs to the tolerance rather than to the fit.
function (::CategoricalSuitability{A, V})(niche,
                                          pref::AbstractVector{V}) where {A,
                                                                          V}
    return niche in pref ? 1.0 :
           0.0
end

function (::CategoricalSuitability{A, V})(niche, pref::V) where {A, V}
    return niche == pref ? 1.0 : 0.0
end

"""
    NoFitContinuous{A, V} <: AbstractNicheFit{A, V}

The absence of a nichefit between a continuous trait and its environment,
parameterised on its **axis structure** `A` and on any `V`. Returns the value 1.

"""
struct NoFitContinuous{A, V} <: AbstractNicheFit{A, V}
end

# The two-argument form `_suitability` actually calls (`nichefit(dist, current)`), matching every
# other continuous fit. The released three-argument form below is kept for callers of the v0.4.0
# `NoRelContinuous`, which was only ever invoked by hand - `_suitability` has never passed three.
(::NoFitContinuous)(_, _) = 1.0

function (::NoFitContinuous{A, V})(::V, ::V, ::V) where {A, V}
    return 1.0
end

"""
    NoFitCategorical{A, V} <: AbstractNicheFit{A, V}

The absence of a nichefit between a categorical trait and its environment,
parameterised on its **axis structure** `A` and on any `V`. Returns the value 1.

"""
struct NoFitCategorical{A, V} <: AbstractNicheFit{A, V}
end

# The one-argument form `_suitability` calls, matching `CategoricalSuitability` - it discards the
# tolerance's weight, which is exactly what "no fit" means. The two-argument form below is the
# released `NoRelDiscrete` spelling.
(::NoFitCategorical)(_) = 1.0

function (::NoFitCategorical{A, V})(niche::V, pref::V) where {A, V}
    return 1.0
end

"""
    CombiningFit{A, F, C <: NamedTuple} <: AbstractNicheFit{A, C}

Several nichefits, one per environmental layer, together with the function that combines their
per-layer suitabilities into the cell's one suitability. The fits are held by name in a
`NamedTuple`, and must line up name for name with the regime and tolerance collections they are
matched against; a `Tuple` is accepted and named by each member's own niche axis.

`combine` is handed the **whole `NamedTuple` of results**, not folded across them as a binary
operator, so a rule that treats layers differently is expressible and reads by name:

```julia
CombiningFit((summer = ..., winter = ..., moisture = ...)) do s
    return (s.summer + s.winter) * s.moisture
end
```

[`MultiplicativeFit`](@ref) and [`AdditiveFit`](@ref) are the two everyone actually uses, and are
aliases for `combine = prod` and `combine = sum`; over a tuple those are themselves unrolled
folds, so they cost exactly what the `*`/`+` of the fixed-arity fit types cost.
"""
struct CombiningFit{A, F, C <: NamedTuple} <: AbstractNicheFit{A, C}
    combine::F
    nt::C

    # The sole constructor - see `SpeciesRequirementCollection`. Pass a `NamedTuple` of fits to name them,
    # which is what lets `combine` read its argument by name.
    function CombiningFit(combine, nichefits::Union{Tuple, NamedTuple})
        nt = _asnamedtuple(nichefits)
        _checkbacking(AbstractNicheFit, values(nt), "nichefit")
        return new{_axisstructure(values(nt)), typeof(combine), typeof(nt)}(combine,
                                                                            nt)
    end
end

"""    MultiplicativeFit{A, C} - a [`CombiningFit`](@ref) whose per-layer suitabilities are multiplied. """
const MultiplicativeFit{A, C} = CombiningFit{A, typeof(prod), C}

# The aliases fix `combine` but leave the member params free, so no constructor is generated for
# them; forward to the sole `CombiningFit` constructor, in the same one spelling.
MultiplicativeFit(fits::Union{Tuple, NamedTuple}) = CombiningFit(prod, fits)

"""    AdditiveFit{A, C} - a [`CombiningFit`](@ref) whose per-layer suitabilities are added. """
const AdditiveFit{A, C} = CombiningFit{A, typeof(sum), C}

AdditiveFit(fits::Union{Tuple, NamedTuple}) = CombiningFit(sum, fits)

# --- Scoring a species against a cell ---------------------------------------
# The nichefit's own job. Every method dispatches on an `AbstractEcosystem`, or on a layer and a
# tolerance, all of which `Ecology.jl` declares first.

"""
    suitability(eco::AbstractEcosystem, pos::Int64, sp::Int64)
    suitability(eco::AbstractEcosystem, y::Int64, x::Int64, sp::Int64)

Return how well species `sp` suits a grid square - by linear position `pos`, or by cell coordinates
`(y, x)`.

Both forms exist because each has a caller that already holds what it takes: `populate_by_tolerance!`
works in linear indices, while the hot loop already has `(y, x)` for the supply. They are
distinguished by **arity**, so neither shadows the other.
"""
function suitability(eco::AbstractEcosystem, y::Int64, x::Int64, sp::Int64)
    return _suitability(eco.habitat.regime, eco.spplist.tolerance, eco.nichefit,
                        convert_coords(eco, (y, x)), sp)
end

function suitability(eco::AbstractEcosystem, pos::Int64, sp::Int64)
    regime = eco.habitat.regime
    tolerance = eco.spplist.tolerance
    nichefit = eco.nichefit
    return _suitability(regime, tolerance, nichefit, pos, sp)
end

"""
    nichefitcombine(nichefit::AbstractNicheFit)

Return the function that combines a nichefit's per-layer suitabilities into one number. It takes
the whole `NamedTuple` of per-layer results - see [`CombiningFit`](@ref) - so a single fit
"combines" by taking the only result, and `_suitability` needs no separate leaf path.
"""
nichefitcombine(nichefit::CombiningFit) = getfield(nichefit, :combine)

nichefitcombine(::AbstractNicheFit) = only

# The axis's density width expressed in the frame `V`, as a bare number - `1.0` when the axis
# declares none, so the multiplication is a no-op rather than a branch in the hot loop.
function _densityscale(::Type{A}, ::Type{V}) where {A, V}
    return _asscale(densitywidth(A), V)
end

# The width side is typed `::Unitful.Quantity`, not left free - with a free `w` the `Nothing`
# method and the `V <: Quantity` one are **ambiguous** for an axis with no width on a unitful
# frame, which is every legacy root-axis fit. Caught by `test_deprecations`/`test_NicheFit`.
_asscale(::Nothing, ::Type) = 1.0

function _asscale(w::Unitful.Quantity, ::Type{V}) where {V <: Quantity}
    return ustrip(unit(V), w)
end

_asscale(w::Unitful.Quantity, ::Type) = ustrip(w)

# **A DIMENSIONLESS axis's width is a bare number, not a `Quantity`** - `1.0NoUnits` evaluates to
# `1.0::Float64`, because Unitful collapses a quantity with no units to the number itself. All eight
# dimensionless axes (`SurfaceArea`, `TypologyAxis`, `Heterogeneity`, ...) declare their width that
# way, so without this method every continuous niche on one of them dies with
# `MethodError: no method matching _asscale(::Float64, ::Type{Float64})` in the hot loop.
# Nothing to strip: the width is already a bare number in the only frame such an axis has.
_asscale(w::Real, ::Type) = float(w)

# Convert `current` to the frame `V` expects before stripping. When `V` is a concrete `Quantity` type
# (the real case, via `_defaultsuitability`'s tolerance-derived `V`), this is a real `uconvert` - so a
# dimensionally-compatible-but-differently-scaled regime is corrected rather than silently wrong, and a
# genuinely incompatible dimension still throws, just via `uconvert` instead of dispatch. `V` can also be
# an abstract dimension (e.g. `Unitful.Temperature`) or a bare non-quantity type - neither names a single
# fixed target scale to convert to, so just strip as before.
_toframe(::Type{V}, current) where {V <: Quantity} = ustrip(unit(V), current)

_toframe(::Type{V}, current) where {V} = ustrip(current)

iscontinuous(nichefit::NicheSuitability) = true

iscontinuous(nichefit::CategoricalSuitability) = false

iscontinuous(nichefit::NoFitContinuous) = true

iscontinuous(nichefit::NoFitCategorical) = false

# Named access to the members, as for `LayerCollection`. Everything - `:combine` included - is
# forwarded to the backing, so no member name is reserved; `nichefitcombine` reaches the combining
# function with `getfield`.

# A nichefit's axis structure, straight off the type parameter - the third of the three families
# that carry one, and what lets `_checkmembers` compare all three sides rather than only two.
axisof(::AbstractNicheFit{A}) where {A} = A

# The nichefit as a side of a pairing check - see `_checkaligned` (`collections.jl`). It carries a
# real axis, so the condition line is checked three ways, tolerance against regime against nichefit,
# rather than two.
_nichefitside(f::AbstractNicheFit) = _side(f, "trait nichefit", true)

# One arity-general method for every collection: evaluate each layer's own suitability, then hand
# the whole `NamedTuple` of results to the nichefit's combining function. `_zipmap` unrolls at
# compile time and the names come from the type, so both are free - this is the hot path, per cell
# per species.
function _suitability(regime::LayerCollection{Condition},
                      tolerance::SpeciesRequirementCollection{Condition},
                      nichefit::NF,
                      pos::Int64,
                      sp::Int64) where {NF <: AbstractNicheFit}
    results = _zipmap(values(regime), values(tolerance),
                      values(nichefit)) do r, t, f
        return _suitability(r, t, f, pos, sp)
    end
    return nichefitcombine(nichefit)(NamedTuple{keys(nichefit)}(results))
end

function _suitability(regime::ContinuousRegime,
                      tolerance::NicheTolerance,
                      nichefit::NF,
                      pos::Int64,
                      sp::Int64) where {NF <: AbstractNicheFit}
    (y, x) = convert_coords(pos, Base.size(regime.matrix, 1))
    h = _cellvalue(regime, y, x)
    # Fetch the species' pre-built response distribution and evaluate it via `nichefit` - no per-call
    # construction or allocation (the distributions were built once when the `NicheTolerance` was made).
    # One method serves a static and a time-varying regime alike: a series writes the current slice
    # into the layer's own matrix, so by the time this reads it there is nothing left to choose.
    return nichefit(_speciestolerance(tolerance, sp), h)
end

function _suitability(regime::CategoricalRegime,
                      tolerance::AbstractCategoricalTolerance,
                      nichefit::NF,
                      pos::Int64,
                      sp::Int64) where {NF <: AbstractNicheFit}
    (y, x) = convert_coords(pos, Base.size(regime.matrix, 1))
    # The tolerance answers with the weight and the nichefit reads it, which is why one method serves
    # every categorical tolerance: a single preferred class is a one-element set, so there is no
    # second shape to dispatch on.
    return nichefit(_categoryweight(tolerance, sp, _cellvalue(regime, y, x)))
end

# ---------------------------------------------------------------------------
# The accessor interface
# ---------------------------------------------------------------------------
# **Two layers, because the hot path and a user want different things.** The private accessors
# return a bare value shaped for their one call site; the public ones return a `NamedTuple` keyed by
# the collection's own names. Neither is a compromise: `_suitability` genuinely wants a number, and a
# user genuinely wants to know which layer each value came from.
#
# The private names deliberately do **not** mirror the public ones. `_cellvalue` is singular where
# the public pair is `cellregime`/`cellsupply`, so it cannot be mistaken for a copy of them - and it
# is singular for a reason: a regime and a supply are both `AbstractLayer`s with a `.matrix`, so
# there is nothing to duplicate.

# The value a layer holds at one cell - one function for both roles.
_cellvalue(layer::AbstractLayer, y::Int64, x::Int64) = layer.matrix[y, x]

# What `_suitability` needs from a tolerance for one species: exactly one thing per kind, so there is
# one method each and no merge to agonise over.
function _speciestolerance(tolerance::NicheTolerance, sp::Int64)
    return tolerance.dists[sp]
end

function _speciestolerance(tolerance::SimpleCategoricalTolerance, sp::Int64)
    return tolerance.vals[sp]
end

# The suitability weight species `sp` gets in category `c` - the single question every categorical
# tolerance must answer, and the reason there is now one `_suitability` method rather than two.
# It lives on the tolerance, not the nichefit, because the general (per species per category) case
# is only expressible here: a nichefit field would be one value for the whole run.
# Allocation-free: `in` over a small `Vector` neither allocates nor is measurably slower than the
# scalar `==` it replaces.
function _categoryweight(tolerance::SimpleCategoricalTolerance, sp::Int64, c)
    return c in tolerance.vals[sp] ? 1.0 : tolerance.penalty
end

# The per-species demand for one resource, as the raw unitful rate.
_speciesdemand(demand::AbstractDemand, sp::Int64) = demand.resource[sp]

# --- Inferring the fit, and checking the pair share an axis --------------------

# A trait may only be matched to a layer on the same niche axis.
# **One policy, shared with `_checkmembers`**: `_axesagree` (`collections.jl`) is the single answer to
# "do these two axes agree?", and both routes ask it, so the question cannot have two answers
# depending on which entry point a caller came in by. The root is not special here either - see
# `_axesagree`.
function _checkaxis(t, h)
    at, ah = axisof(t), axisof(h)
    _axesagree(at, ah) ||
        error("A trait on the `$(nameof(at))` niche axis cannot be matched to a layer on the " *
              "`$(nameof(ah))` axis - paired layers must be on the *same* axis. Declare the same " *
              "axis on both sides (a grouping axis is fine, so long as both say it), or check the " *
              "order of the species niches against the environment layers.")
    return nothing
end

# Infer the trait-environment nichefit from the trait type, checking the niche axes.
# A `NicheTolerance`'s response distribution is built per-species in `_suitability`; its nichefit's `V`
# is the tolerance's own frame (what the distribution's parameters are actually expressed in) - not the
# regime's incidental unit, so a genuine tolerance/regime unit disagreement is caught by the pairing
# check at Ecosystem construction (`_checkaligned`, `collections.jl`) instead of the nichefit
# silently mirroring whatever the regime happens to be.
function _defaultsuitability(t::NicheTolerance, regime)
    _checkaxis(t, regime)
    return NicheSuitability{axisof(t), eltype(t)}()
end

# Same reasoning as the `NicheTolerance` method above: `V` comes from the tolerance, not the regime.
# One method where there were two, and it now checks the axis like its continuous sibling - which
# it could not before, because neither categorical tolerance carried one.
function _defaultsuitability(t::AbstractCategoricalTolerance, regime)
    _checkaxis(t, regime)
    return CategoricalSuitability{axisof(t), eltype(t)}()
end

# A collection of tolerances over a matching collection of regimes infers each sub-tolerance's
# nichefit and combines them multiplicatively. The derived fit takes the tolerance's own backing -
# positional stays positional, named keeps its names - so the three stay comparable.
function _defaultsuitability(t::SpeciesRequirementCollection{Condition},
                             h::LayerCollection{Condition})
    ts, hs = values(t), values(h)
    length(ts) == length(hs) ||
        error("A collection of $(length(ts)) tolerances cannot be matched to a collection of " *
              "$(length(hs)) regimes - pass one tolerance per environment layer.")
    return MultiplicativeFit(_likebacking(t,
                                          _zipmap(_defaultsuitability, ts,
                                                  hs)))
end

function _defaultsuitability(traits::AbstractTolerance, regime)
    return error("Cannot infer a trait nichefit for traits of type $(typeof(traits)); pass one explicitly with `nichefit = ...`.")
end
