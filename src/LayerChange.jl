# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A layer's change as it is actually held: a pure function of elapsed time, materialised
# from a spec and unit-checked against the layer it drives.

using Unitful
using Dates: Dates
using Unitful.DefaultSymbols
using DimensionalData
using DimensionalData.Lookups: NoLookup
using EcoSISTEM.Units

"""
    AbstractLayerChange{M <: AbstractChangeMode}

Abstract supertype of the materialised per-timestep change rules a layer can hold, parameterised
by the [`AbstractChangeMode`](@ref) that says how to read its values.

A change is applied as a pure function of `(layer, elapsed, timestep)` — never of the ecosystem —
which is what lets MPI apply it redundantly on every rank without diverging.
"""
abstract type AbstractLayerChange{M <: AbstractChangeMode} end

"""    NoLayerChange <: AbstractLayerChange{NoChange} — the layer never changes. """
struct NoLayerChange <: AbstractLayerChange{NoChange} end

"""
    SteadyLayerChange{V} <: AbstractLayerChange{RateChange}

A constant rate of change: the layer gains `value * timestep` every step. `value` is a scalar or a
per-cell matrix, already converted to the layer's `changeunit` when the change was attached.

The mode is fixed to [`RateChange`](@ref) because only a rate accumulates. A constant *value* —
absolute or relative — writes the same thing every step and is idempotent after the first, which
makes it a one-off operation on the ecosystem rather than a per-timestep change of the layer.
"""
struct SteadyLayerChange{V} <: AbstractLayerChange{RateChange}
    value::V
end

"""
    PatternedLayerChange{M <: AbstractChangeMode, F, V, B} <: AbstractLayerChange{M}

A change following an arbitrary `shape` of elapsed simulation time: the layer gets
`amplitude * shape(elapsed / timescale)`, read according to `M`. `shape` defaults to
[`sinusoidal`](@ref), one full cycle per `timescale`; a sigmoid, ramp or step function is equally
valid.

The phase handed to `shape` is deliberately **not** wrapped into `[0, 1)`. A sinusoid is periodic
by itself so wrapping would gain nothing, while a sigmoid is not periodic at all and wrapping would
destroy it. `timescale` is therefore the elapsed time mapping to one unit of `shape`'s argument — a
period for a sinusoid, a transition width for a sigmoid.

Under [`RelativeChange`](@ref) the pattern is added to `baseline`, captured from the layer when the
change was attached so that the result stays a pure function of elapsed time rather than compounding
on what the previous step wrote. The other modes leave `baseline` as `nothing`.
"""
struct PatternedLayerChange{M <: AbstractChangeMode, F, V, B} <:
       AbstractLayerChange{M}
    shape::F
    amplitude::V
    timescale::typeof(1.0s)
    baseline::B
end

"""
    SeriesLayerChange{M <: AbstractChangeMode, E <: AbstractSeriesEnd, C <: AbstractSeriesCalendar, D, B} <: AbstractLayerChange{M}

A stack of stored slices, indexed by elapsed simulation time: a slice is current from its own stored
time until the next slice's, so at `elapsed` the layer takes the last slice at or before
`origin + elapsed`, read according to `M`. Under
[`AbsoluteChange`](@ref) the layer *is* that slice (a read climate series); under
[`RelativeChange`](@ref) the slice is added to the layer's captured values.

`times` holds one time per slice in the same order, and `origin` is the point in that coordinate
which elapsed time is measured from — the first stored time unless the caller says otherwise, so
that a series starts at its own beginning whatever its axis is anchored to.

`calendar` records what the coordinates *mean* — see [`AbstractSeriesCalendar`](@ref). It is what
[`build_ecosystem`](@ref) reads when it resolves a run's epoch, and what it re-points the series
against: resolving an epoch rewrites `origin`, leaving the stored slices and the hot path untouched.

**Outside its own span a series contributes nothing**, and the layer is whatever it would be
without it: before the first slice always, and past the last under [`RevertToLayer`](@ref). `baseline`
is what that costs — the layer's values as they stood when the change was attached, kept for the
modes that need something to fall back *to*.

Indexing by *time* rather than by a step counter is the substance of this type: it makes the result
independent of the timestep (one twelve-month step and twelve one-month steps land on the same
slice), which a cursor advanced once per call can never be. `atend` decides what happens past the
last slice — see [`AbstractSeriesEnd`](@ref).
"""
struct SeriesLayerChange{M <: AbstractChangeMode, E <: AbstractSeriesEnd,
                         C <: AbstractSeriesCalendar, D, B} <:
       AbstractLayerChange{M}
    slices::D
    times::Vector{typeof(1.0s)}
    origin::typeof(1.0s)
    atend::E
    calendar::C
    baseline::B
end

"""
    SumOfLayerChanges{M <: AbstractChangeMode, P <: Tuple, B} <: AbstractLayerChange{M}

Several changes driving one layer, added together — a stored monthly series plus a multi-year
warming trend, say, where the trend offsets the whole seasonal pattern rather than contradicting it.

Composition is a sum of **values as functions of elapsed time**, evaluated once per step, never a
sequence of mutations. Chaining the applications instead would silently fail: an absolute write each
step erases whatever a rate had accumulated by then.

`M` follows from the parts. A part that is a *position* (a read series) makes the whole an
[`AbsoluteChange`](@ref) — there can be at most one, since two positions cannot be added — and
`parts` keeps it first so the fold adds intervals to it in an order affine units accept. With no
such part the whole is a [`RelativeChange`](@ref) over the layer's captured `baseline`. A
[`RateChange`](@ref) part contributes its integral, so it must have one in closed form: a steady
rate does (`value × elapsed`), a patterned one would need the integral of its own shape and is
refused rather than approximated.
"""
struct SumOfLayerChanges{M <: AbstractChangeMode, P <: Tuple, B} <:
       AbstractLayerChange{M}
    parts::P
    baseline::B
end

"""
    LegacyLoss{V} <: AbstractLayerChange{NoChange}

Carries the loss rate of the deprecated `HabitatLoss` change function. Its mode is
[`NoChange`](@ref): the rate is a plain per-time probability, in no way the layer's own unit.

**Transitional, and not a layer change at all.** Habitat loss mutates the *ecosystem* — supply
and abundances — from a layer's change slot, and draws randomly while doing it, so it satisfies
neither invariant the other changes do. It survives only to keep `HabitatLoss` callable, and is
superseded by an explicit cell-deactivating intervention.
"""
struct LegacyLoss{V} <: AbstractLayerChange{NoChange}
    rate::V
end

# `<: AbstractSpec` because a `Varying` *is* a build-time layer recipe — a spec plus the
# declaration of how it changes — and saying so is what lets `regime`/`supply` be typed
# `Union{AbstractSpec, Tuple, NamedTuple}` in the builders' own signatures rather than left to a
# runtime failure. Safe by construction: nothing in the package dispatches on `AbstractSpec`, which
# is a documentation root, so this can add no method ambiguity and change no existing call.
"""
    Varying(spec, change)

Declare that the layer built from `spec` carries `change`. `change` is a change recipe
([`ReplaceWith`](@ref), [`OffsetBy`](@ref), [`IncrementBy`](@ref)) or an already materialised
[`AbstractLayerChange`](@ref).

Pass it to `GridHabitat`'s `regime` or `supply` keyword. A `StudyArea` ignores the wrapper
entirely — a change has no meaning until there is a grid to apply it to, so the area is decided from
the spec alone.

```julia
GridHabitat(regime = Varying(SourceSpec(WorldClim{BioClim}, :bio1),
                                   IncrementBy(0.02K/yr)),
                  supply = GradientSpec(…), area = area)
```

Each layer names its own change, so a multi-variable regime wraps its *elements*, never the tuple:
`(Varying(temp, …), rain)`, not `Varying((temp, rain), …)`.
"""
struct Varying{S, C} <: AbstractSpec
    spec::S
    change::C

    function Varying(spec,
                     change::Union{AbstractChangeSpec, AbstractLayerChange})
        spec isa Varying &&
            error("`Varying` cannot wrap another `Varying`: a layer carries exactly one change. " *
                  "Combine the two declarations into one instead.")
        spec isa Union{Tuple, NamedTuple} &&
            error("`Varying` declares a change for one layer, so it cannot wrap a tuple of " *
                  "$(length(spec)) of them: a change is checked against a single layer's unit, and " *
                  "one declaration covering several would silently mean different things on each. " *
                  "Wrap the elements instead — `($(join(fill("…", length(spec)), ", ")))` with " *
                  "`Varying` on the ones that vary.")
        return new{typeof(spec), typeof(change)}(spec, change)
    end
end

# What a *constant* shape is: a single value, or one value per cell. Named rather than spelled
# `Union{…}` at each of its two `_attachshape` methods, because it is the definition that keeps them
# from colliding with the shape-typed ones — a `PatternedChange` is deliberately neither.
const _ConstantShape = Union{Number, AbstractArray}

# ══ Functions ══════════════════════════════════════════════════════════════════════════════════

# ---------------------------------------------------------------------------
# The unit contract
# ---------------------------------------------------------------------------
# What unit a change's values must be in, decided by the layer it is attached to and the change's
# mode — never hard-coded per change, which is how a change's unit and its layer's unit drift apart.

"""
    changeunit(mode::AbstractChangeMode, layer::AbstractLayer)

Return the unit a change's values must be in to be applied to `layer` under `mode`, or `nothing`
for a [`NoChange`](@ref) change, which carries no values in the layer's unit at all.

The three value-carrying modes are exactly the three answers to "what kind of quantity is this?":

| mode | the value is | unit on a `°C` layer |
|---|---|---|
| [`AbsoluteChange`](@ref) | a position — the layer's value *is* this | `°C` |
| [`RelativeChange`](@ref) | an interval from the layer's captured values | `K` |
| [`RateChange`](@ref) | an interval per unit time | `K/s` |

The absolutising in the latter two is load-bearing, not cosmetic: a temperature is a position but a
*change* in temperature is a width, and Unitful rejects adding two affine positions outright
(`10.0°C + 1.0°C` throws `AffineError`). Reading a relative change in `°C` would therefore not be
merely imprecise, it would not run.
"""
changeunit(::NoChange, ::AbstractLayer) = nothing
changeunit(::AbsoluteChange, layer::AbstractLayer) = unit(eltype(layer.matrix))
function changeunit(::RelativeChange, layer::AbstractLayer)
    return Unitful.absoluteunit(unit(eltype(layer.matrix)))
end
function changeunit(::RateChange, layer::AbstractLayer)
    return Unitful.absoluteunit(unit(eltype(layer.matrix))) / s
end

# Convert `x` (a scalar or a per-cell matrix) into the unit `mode` demands of `layer`, once, at
# attach. This is the check that `LayerUpdate` did not do: it validated the rate against a dimension
# passed in by hand and then discarded it, never consulting the layer or the change function at all.
# The error names both units and the layer's axis, because the bare `DimensionError` a raw `uconvert`
# throws says nothing about which of the two is wrong.
function _tochangeunit(mode::AbstractChangeMode, layer::AbstractLayer, x)
    u = changeunit(mode, layer)
    if dimension(eltype(x)) != dimension(u)
        error("a $(nameof(typeof(mode))) on a $(nameof(axisof(layer))) layer must be in " *
              "$u (the layer holds $(unit(eltype(layer.matrix)))), but got " *
              "$(unit(eltype(x))), which has dimension $(dimension(eltype(x))) " *
              "rather than $(dimension(u))")
    end
    return uconvert.(u, x)
end
# A `NoChange` change has no values to convert, so nothing to check against.
_tochangeunit(::NoChange, ::AbstractLayer, x) = x

# ---------------------------------------------------------------------------
# Change recipes — what a caller writes, before it meets a layer
# ---------------------------------------------------------------------------
# A recipe names a *shape* and a *mode* but no unit-checked values; `_attachchange` turns it into a
# materialised `AbstractLayerChange` against the layer it is destined for. Deliberately thin: this
# is only what the current attachment sites and the deprecation shims need, and the user-facing
# spec surface arrives with the attachment wrapper.

"""
    sinusoidal(phase)

Return one full sine cycle over a unit of dimensionless `phase`, the default shape of a
[`PatternedChange`](@ref) or [`PatternedLayerChange`](@ref).

`phase` is elapsed time divided by the change's `timescale`, so a whole number of timescales returns
to zero. Any function of a dimensionless phase may be used instead — a sigmoid, a ramp, a step.

# Arguments

  - `phase`: the dimensionless phase, elapsed time over `timescale`.
"""
sinusoidal(phase) = sinpi(2 * phase)   # `sinpi`, for accuracy at exact half-turns

# `Varying` prints as the two-argument call that builds it, and both halves print through their own
# compact form — which is the case the two-method `show` split exists for.
function Base.show(io::IO, v::Varying)
    return print(io, "Varying(", sprint(show, v.spec), ", ",
                 sprint(show, v.change), ")")
end

# Nesting a combination inside a combination would mean the same thing as one flat sum, so it is
# flattened here rather than left to produce two shapes for one meaning.
_flattenspecs(specs::Tuple) = mapreduce(_asspecs, (a, b) -> (a..., b...), specs)
_asspecs(spec::AbstractChangeSpec) = (spec,)
_asspecs(spec::CombinedChange) = spec.specs

# The mode each recipe fixes. Separate methods rather than a field, so the mode is available as a
# type and can become a change's type parameter.
_changemode(::ReplaceWith) = AbsoluteChange()
_changemode(::OffsetBy) = RelativeChange()
_changemode(::IncrementBy) = RateChange()

# The single funnel where a change is attached to a layer: build it, then check the values it will
# *store* against the layer's own bounds. Building is `_buildchange`, so that the check happens once,
# on the finished change, rather than once per nested part.
# Materialise a change recipe against the layer it will drive, converting and validating its values
# against [`changeunit`](@ref) exactly once. Returns an [`AbstractLayerChange`](@ref).
function _attachchange(spec, layer::AbstractLayer)
    change = _buildchange(spec, layer)
    _checkstoredbounds(change, layer)
    return change
end

# Turn a change recipe into the materialised change it describes, against the layer that will carry
# it. An already-materialised change passes through, which is what lets `Varying` accept either.
function _buildchange(spec::AbstractChangeSpec, layer::AbstractLayer)
    return _attachshape(spec.shape, _changemode(spec), layer)
end
_buildchange(change::AbstractLayerChange, ::AbstractLayer) = change
function _buildchange(spec::CombinedChange, layer::AbstractLayer)
    return _combineparts(map(s -> _buildchange(s, layer), spec.specs), layer)
end

# Assemble already-materialised changes into one that sums them. The mode follows from the parts —
# see `SumOfLayerChanges` — and every part is evaluated once here, at attach, so a part that cannot
# contribute a value function of elapsed time says so now rather than on some later timestep.
function _combineparts(parts::Tuple, layer::AbstractLayer)
    positions = filter(p -> _modeof(p) isa AbsoluteChange, parts)
    length(positions) <= 1 ||
        error("$(length(positions)) of these changes give the layer's value outright, and two " *
              "positions cannot be added — only intervals can be added to a position. Combine one " *
              "`ReplaceWith` with `OffsetBy`/`IncrementBy` changes.")
    # Positions first: the fold is left-to-right, and an affine unit accepts `10.0°C + 1.0K` but
    # not the same sum begun from the interval.
    ordered = (positions...,
               filter(p -> !(_modeof(p) isa AbsoluteChange),
                      parts)...)
    foreach(p -> _partvalue(p, 0.0s), ordered)
    mode = isempty(positions) ? RelativeChange() : AbsoluteChange()
    baseline = _patternbaseline(mode, layer)
    return SumOfLayerChanges{typeof(mode), typeof(ordered), typeof(baseline)}(ordered,
                                                                              baseline)
end

# A change's mode as a value, recovered from the type parameter it already carries.
_modeof(::AbstractLayerChange{M}) where {M} = M()

# What a part contributes to the sum at `elapsed`. For a position or an interval that is just its
# value; a rate contributes its *integral*, which is why only a steady one qualifies.
_partvalue(change::AbstractLayerChange, elapsed) = _changevalue(change, elapsed)
_partvalue(change::SteadyLayerChange, elapsed) = change.value .* elapsed
function _partvalue(change::AbstractLayerChange{RateChange}, _)
    return error("a $(nameof(typeof(change))) cannot be combined with other changes: a sum is over " *
                 "values as functions of elapsed time, so a rate contributes its integral. A " *
                 "steady rate has one in closed form (`value × elapsed`); this one would need the " *
                 "integral of its own shape, which is not generally available, and approximating " *
                 "it would make the result depend on the timestep. Use a steady rate.")
end
function _partvalue(change::AbstractLayerChange{NoChange}, _)
    return error("a $(nameof(typeof(change))) carries no values in the layer's unit, so there is " *
                 "nothing for it to contribute to a sum of changes.")
end

# A bare quantity or matrix under `IncrementBy` is a steady drift — the only mode in which a
# constant does something new each step.
function _attachshape(shape::_ConstantShape, mode::RateChange,
                      layer::AbstractLayer)
    value = _tochangeunit(mode, layer, shape)
    return SteadyLayerChange{typeof(value)}(value)
end
# …and under the other two it is a one-off, so it is refused here rather than silently installed as
# a change that rewrites the same values forever.
function _attachshape(shape::_ConstantShape, mode::AbstractChangeMode,
                      ::AbstractLayer)
    return error("`$(nameof(typeof(mode)))` with a constant ($shape) is not a layer change: " *
                 "writing the same values every timestep is idempotent after the first, so it " *
                 "is a one-off operation on the ecosystem. Use `IncrementBy` for a steady rate " *
                 "of change, or a shape that varies with time.")
end
# These two dispatch on the shape being a *constant*, not on `Any`. Against `Any` they were
# ambiguous with every shape-typed method under `RateChange` — one method narrowing the shape and
# the other the mode, neither more specific — so `IncrementBy(PatternedChange(…))`, an oscillating
# rate, always failed despite being documented as valid. Anything that is neither a constant nor a
# known shape falls through to here.
function _attachshape(shape, mode::AbstractChangeMode, ::AbstractLayer)
    return error("$(typeof(shape)) is not a change shape. Use a constant (with `IncrementBy`, " *
                 "for a steady rate), a `PatternedChange` for a shape of elapsed time, or a " *
                 "`SeriesChange` for a stack of stored slices.")
end
function _attachshape(shape::PatternedChange, mode::AbstractChangeMode,
                      layer::AbstractLayer)
    amplitude = _tochangeunit(mode, layer, shape.amplitude)
    baseline = _patternbaseline(mode, layer)
    return PatternedLayerChange{typeof(mode), typeof(shape.shape),
                                typeof(amplitude), typeof(baseline)}(shape.shape,
                                                                     amplitude,
                                                                     uconvert(s,
                                                                              float(shape.timescale)),
                                                                     baseline)
end

function _attachshape(shape::SeriesChange, mode::AbstractChangeMode,
                      layer::AbstractLayer)
    coords = _seriescoords(shape.source, shape.times, shape.calendar)
    _checkseriesend(shape.atend, coords.times)
    slices = _tochangeunit(mode, layer, shape.source)
    _checkseriesgrid(slices, layer)
    origin = _seriesorigin(coords.calendar, coords.times, shape.origin)
    baseline = _seriesbaseline(mode, layer)
    # The slices are stored stripped of their dims (`parent` is a no-op on a plain array): each is
    # broadcast into the layer's own matrix, which carries the grid's real lookups, and matching two
    # sets of lookups elementwise would reject a legitimately-independent series over an identical
    # grid. `_checkseriesgrid` has just established that they are the same shape, which is the part
    # that actually has to hold.
    return SeriesLayerChange{typeof(mode), typeof(shape.atend),
                             typeof(coords.calendar), typeof(parent(slices)),
                             typeof(baseline)}(parent(slices), coords.times,
                                               origin, shape.atend,
                                               coords.calendar, baseline)
end

# A series' elapsed-time coordinates and its calendar, resolved together because a dated source
# settles both at once: the same lookup supplies the dates that become the coordinates *and* the
# start date the calendar has to record. A named tuple so the two cannot be swapped at the call site.
function _seriescoords(source, times, calendar)
    found = isnothing(times) ? _sourcetimes(source) : _giventimes(times)
    resolved = isnothing(calendar) ? found.calendar : calendar
    checked = _checkseriestimes(found.times, source)
    _checkcalendar(resolved, checked)
    return (times = checked, calendar = resolved)
end

# The stored times, validated. Strictly increasing, because "the slice current at this time" is only
# well defined if no two slices claim the same instant.
function _checkseriestimes(given, source)
    ndims(source) == 3 ||
        error("a series' source must be a stack of slices — a 3-dimensional `(Y, X, Ti)` array — " *
              "but got a $(ndims(source))-dimensional one. A single slice is not a change.")
    length(given) == size(source, 3) ||
        error("a series needs one time per slice, but got $(length(given)) times for " *
              "$(size(source, 3)) slices.")
    resolved = collect(uconvert.(s, float.(given)))
    all(>(0.0s), diff(resolved)) ||
        error("a series' times must strictly increase, so that each slice covers a distinct " *
              "stretch of the simulation, but got $(given).")
    return resolved
end

# Times passed explicitly by the caller. Dates are still dates — giving a series real times is the
# documented way to say when an otherwise undated stack really begins — so they resolve exactly as a
# dated lookup would; anything else carries no calendar identity of its own.
_giventimes(times) = (times = times, calendar = UndatedSeries())
_giventimes(times::AbstractVector{<:Dates.TimeType}) = _datedtimes(times)

# A source's own slice times, taken from its `Ti` lookup.
function _sourcetimes(source::DimensionalData.AbstractDimArray)
    hasdim(source, Ti) ||
        error("a series' source has no `Ti` dimension to take slice times from — pass `times = ` " *
              "explicitly, or read the source with a time axis.")
    return _lookuptimes(DimensionalData.lookup(source, Ti))
end
# A bare 3-D array carries no coordinates at all, and monthly is what such a stack has always meant
# here. Note that only the **spacing** is monthly: the series is still `UndatedSeries`, so nothing
# phase-locks it to January, and a bare 10-slice stack is ten slices a month apart rather than
# "months 1 to 10 of the year".
function _sourcetimes(source::AbstractArray{<:Any, 3})
    return (times = axes(source, 3) .* month_mean_duration,
            calendar = UndatedSeries())
end

# Real dates become elapsed coordinates measured from the first slice, and that first date is what
# the calendar keeps: it is the series' entire calendar identity, and what an epoch is resolved
# against. Converting once, here, is what keeps date arithmetic out of the hot path altogether.
function _datedtimes(dates)
    start = first(dates)
    return (times = [uconvert(s,
                              Dates.value(Dates.Millisecond(d - start)) *
                              Unitful.ms) for d in dates],
            calendar = DatedSeries(start))
end

# Monthly is the spacing, never the meaning: an axis-less stack gets monthly coordinates for
# back-compatibility but stays `UndatedSeries`, so no epoch phase-locks it. A 10-slice synthetic
# stack read as "months 1–10 of the year" was the bug this distinction exists to prevent.
function _lookuptimes(lookup::NoLookup)
    return (times = eachindex(lookup) .* month_mean_duration,
            calendar = UndatedSeries())
end
# Elapsed-time coordinates are genuinely ambiguous — `1, 2, 3` months is equally "the first three
# months of the year" and "three months into my experiment" — so they infer the reading that phases
# nothing, and a climatology opts in with `calendar = MonthOfYearSeries()`.
function _lookuptimes(lookup::AbstractVector{<:Unitful.Time})
    return (times = collect(lookup), calendar = UndatedSeries())
end
# A calendar axis is now usable: its dates become elapsed coordinates and its start is kept, so a
# run's epoch can place it. Before an epoch existed this was rejected outright, since the simulation
# clock counted from construction and "18 months into the run" named no date.
_lookuptimes(lookup::AbstractVector{<:Dates.TimeType}) = _datedtimes(lookup)
function _lookuptimes(::AbstractVector{T}) where {T}
    return error("a series' `Ti` lookup holds $T, which is not a time. Pass `times = ` to give " *
                 "each slice an elapsed time explicitly.")
end

# The zero point of the lookup: the coordinate elapsed time zero corresponds to. `origin` is the
# knob for an undated series and refused for the other two, where the epoch fixes the phase and a
# second setting for the same thing could only contradict it.
# The `(UndatedSeries, ::Nothing)` method is not redundant — without it that call is genuinely
# ambiguous, the general method being more specific in its last argument and the undated one in its
# first, so neither wins. The same shape as `_stackcoords`'s ambiguity in `datasetread.jl`.
_seriesorigin(::AbstractSeriesCalendar, times, ::Nothing) = first(times)
_seriesorigin(::UndatedSeries, times, ::Nothing) = first(times)
_seriesorigin(::UndatedSeries, times, origin) = uconvert(s, float(origin))
function _seriesorigin(calendar::AbstractSeriesCalendar, times, origin)
    return error("`origin` is not accepted for a $(nameof(typeof(calendar))): its slices have a " *
                 "calendar identity of their own, so where elapsed time zero falls is fixed by the " *
                 "run's epoch (see `build_ecosystem`), and an `origin` beside it could only " *
                 "contradict it. `origin` is for an `UndatedSeries`, which has no other way to say " *
                 "where zero sits. To place an undated series at a real date, pass `times = ` as " *
                 "dates instead — that says what every slice is, not just where zero is.")
end

# A month-of-year series' coordinates have to *be* month numbers, because that is how the epoch
# reads them back. Checked where the calendar is declared rather than at epoch resolution, so a
# mislabelled series is refused at the line that mislabelled it.
_checkcalendar(::AbstractSeriesCalendar, times) = nothing
function _checkcalendar(::MonthOfYearSeries, times)
    months = _monthnumbers(times)
    (all(m -> 1 <= m <= 12, months) &&
     all(t -> isapprox(t / month_mean_duration,
                       round(ustrip(NoUnits, t / month_mean_duration)),
                       atol = 1.0e-6), times)) ||
        error("`calendar = MonthOfYearSeries()` says these slices are months of the year, so each " *
              "time must be a whole month number from 1 to 12 (`n * month_mean_duration`, as a " *
              "monthly climatology is read), but got $(times). Use `UndatedSeries()` for a series " *
              "whose coordinates are plain offsets.")
    return nothing
end

# Month numbers back out of month-of-year coordinates. Rounding rather than truncating: the
# coordinates are built as `n * month_mean_duration` and float division does not return `n` exactly.
function _monthnumbers(times)
    return [round(Int, ustrip(NoUnits, t / month_mean_duration)) for t in times]
end

# A repeating series needs a period, and its period is only defined if its slices are evenly spaced —
# an irregular axis (an ERA read, or explicit irregular `times`) has no turn length to derive, so it
# is refused rather than guessed at. The other two policies need nothing but an end.
_checkseriesend(::AbstractSeriesEnd, times) = nothing
function _checkseriesend(::RepeatAtEnd, times)
    length(times) > 1 ||
        error("a single-slice series cannot repeat: one slice has no spacing, so there is no " *
              "period to cycle over. Use `atend = HoldAtEnd()`.")
    gaps = diff(times)
    all(g -> isapprox(g, first(gaps), rtol = 1e-6), gaps) ||
        error("`atend = RepeatAtEnd()` needs evenly spaced slices to derive a period from, but this " *
              "series' gaps range from $(minimum(gaps)) to $(maximum(gaps)). Use `atend = HoldAtEnd()`, " *
              "or pass evenly spaced `times = `.")
    return nothing
end

# A series' slices must be on the layer's own grid. Checked at attach because a mismatch would
# otherwise surface as an opaque broadcast error on the first timestep, a long way from the line
# that caused it.
function _checkseriesgrid(slices, layer::AbstractLayer)
    size(slices)[1:2] == size(layer.matrix) ||
        error("a series' slices are $(size(slices)[1:2]) but the layer is " *
              "$(size(layer.matrix)) — the series is on a different grid.")
    return nothing
end

# The values a relative change is measured from — a snapshot taken at attach, so the change stays a
# pure function of elapsed time rather than compounding on what it wrote last step. The baseline
# belongs to the *mode*: only a relative change has one.
_patternbaseline(::RelativeChange, layer::AbstractLayer) = copy(layer.matrix)
_patternbaseline(::AbstractChangeMode, ::AbstractLayer) = nothing

# A *series* keeps the same snapshot for one further reason, so it needs it under one further mode.
# Outside its own span a series contributes nothing and the layer stands as it would without it —
# which under `AbsoluteChange` means there has to be something to stand *as*. Without it a lead-in
# would leave whatever the previous write happened to be, and `RevertToLayer` would have nothing to
# revert to. `RelativeChange` already keeps it, and `RateChange` needs none: contributing nothing is
# a zero rate, not a grid.
# Captured unconditionally rather than only when a lead-in or `RevertToLayer` is in play, because
# neither is knowable at attach — an epoch rewrites `origin` afterwards (`_repointseries!`) and can
# turn a series with no lead-in into one with a lead-in. It is one grid beside the N a series already
# stores, so the cost is a fraction of what is there anyway.
function _seriesbaseline(mode::AbstractChangeMode, layer::AbstractLayer)
    return _patternbaseline(mode,
                            layer)
end
_seriesbaseline(::AbsoluteChange, layer::AbstractLayer) = copy(layer.matrix)

# ---------------------------------------------------------------------------
# Applying a change
# ---------------------------------------------------------------------------
# A layer change is a pure function of `(layer, elapsed, timestep)` — never of the ecosystem. That
# is what lets the MPI update loop apply layer changes redundantly on every rank and stay in step:
# with no ecosystem and no randomness, every rank computes the same thing. `elapsed` is the clock
# *after* this step's advance (see `update!`), so a change sees the time it is changing to.

# Dispatch is on the *mode*, carried in the supertype, so each of the three ways of reading a value
# is written once regardless of which shape produced it. Each shape supplies `_changevalue` instead.
_applychange!(::AbstractLayerChange{NoChange}, ::AbstractLayer, _, _) = nothing

# Each write is followed by `_enforcebounds!`, uniformly rather than only for the rate mode. An
# absolute *series* is already refused at attach if its stored slices go negative, but a
# `PatternedLayerChange{AbsoluteChange}` computes its values and cannot be checked ahead of time — so
# the write sites are the only place that catches every case.
function _applychange!(change::AbstractLayerChange{AbsoluteChange},
                       layer::AbstractLayer, elapsed::Unitful.Time, _)
    layer.matrix .= _changevalue(change, elapsed)
    return _enforcebounds!(layer)
end
function _applychange!(change::AbstractLayerChange{RelativeChange},
                       layer::AbstractLayer, elapsed::Unitful.Time, _)
    layer.matrix .= change.baseline .+ _changevalue(change, elapsed)
    return _enforcebounds!(layer)
end
function _applychange!(change::AbstractLayerChange{RateChange},
                       layer::AbstractLayer, elapsed::Unitful.Time,
                       timestep::Unitful.Time)
    layer.matrix .+= _changevalue(change, elapsed) .* timestep
    return _enforcebounds!(layer)
end

# What a change is worth at `elapsed`, before the mode decides what to do with it.
_changevalue(change::SteadyLayerChange, _) = change.value
function _changevalue(change::PatternedLayerChange, elapsed::Unitful.Time)
    phase = uconvert(NoUnits, elapsed / change.timescale)
    return change.amplitude .* change.shape(phase)
end
function _changevalue(change::SeriesLayerChange, elapsed::Unitful.Time)
    _inspan(change, change.origin + elapsed) || return _outofspan(change)
    return view(change.slices, :, :, _seriesindex(change, elapsed))
end

# Whether a series has anything to say at `t`. Before its first slice it never does — a series that
# has not started yet is silent, which needs no policy because there is only one reading of it. Past
# its last slice it depends on the end policy, and only `RevertToLayer` steps aside.
function _inspan(change::SeriesLayerChange, t)
    _atorafterstart(change, t) || return false
    return _withinend(change.atend, change, t)
end

# Carries the same `_DRIFT` whisker as `_current`, and for the same reason: an accumulated elapsed
# time lands a hair short of a stored coordinate, and without the tolerance a series would count as
# not-yet-started on the very step that reaches it.
function _atorafterstart(change::SeriesLayerChange, t)
    times = change.times
    t >= first(times) && return true
    length(times) == 1 && return false
    return (first(times) - t) <= (times[2] - times[1]) * _DRIFT
end

# Whether a time is still within what the series will answer for. Only `RevertToLayer` has a limit —
# the others hold, repeat or error, all of which are answers — so the default is `true`.
_withinend(::AbstractSeriesEnd, ::SeriesLayerChange, _) = true
function _withinend(::RevertToLayer, change::SeriesLayerChange, t)
    return t <= _seriesreach(change)
end

# What a series is worth outside its own span: "nothing", expressed in each mode's own terms. That
# is why the three `_applychange!` methods need no change at all — writing the captured values,
# adding a zero offset and accumulating a zero rate each leave the layer exactly as it would be with
# no series attached.
_outofspan(change::SeriesLayerChange{AbsoluteChange}) = change.baseline
function _outofspan(change::SeriesLayerChange{RelativeChange})
    return zero(eltype(change.slices))
end
_outofspan(change::SeriesLayerChange{RateChange}) = zero(eltype(change.slices))
# Summed left to right over parts ordered position-first (see `_combineparts`), and broadcast
# throughout because a part contributes either one value or one per cell.
function _changevalue(change::SumOfLayerChanges, elapsed::Unitful.Time)
    values = map(p -> _partvalue(p, elapsed), change.parts)
    return foldl((a, b) -> a .+ b, Base.tail(values), init = first(values))
end

# Which stored slice is current at `elapsed`: the end-of-series policy first maps the requested time
# back into the stored range, then the last slice at or before it is taken.
function _seriesindex(series::SeriesLayerChange, elapsed::Unitful.Time)
    return _current(series.times,
                    _inrange(series.atend, series, series.origin + elapsed))
end

# **A coordinate is when its slice *becomes* current**, so slice `k` is the one in force over
# `[times[k], times[k+1])` and the index is the last coordinate at or before `t`.
#
# **Not a nearest-coordinate rule**, which would treat each coordinate as a slice *centre* and so put
# every transition at a midpoint, half a step early, leaving the first slice half-length.
# `_mkstackaxis` builds the coordinates as slice **identifiers**, `Ti((1:12) .* month_mean_duration)`,
# which is what makes "becomes current" the right reading and "centre" the wrong one.
#
# Nearest is not a careless choice, though, and its reason has to be honoured rather than dropped: an
# accumulated elapsed time does not land exactly on a stored coordinate, so a bare floor turns that
# drift into an off-by-one that only shows up sometimes. `_DRIFT` is the direct answer — a coordinate
# within a whisker *ahead* of `t` is one `t` has really reached. Nearest was avoiding having to write
# this tolerance, at the price of a half-step offset nobody had stated.
#
# The whisker is a fraction of the local gap, so it scales with the series and carries its units: for
# a monthly series it is ~2.5 s, against an accumulated float drift of well under a microsecond over a
# century of stepping, and against a smallest-plausible timestep of an hour. Comfortably larger than
# the error it absorbs and comfortably smaller than anything real.
const _DRIFT = 2^-20
# Which stored slice is current at elapsed time `t`: the last one at or before it, clamped to the
# ends. Clamped rather than refused because running before the first slice or after the last is the
# `atend` policy's business, and this only says which slice those policies are reasoning about.
function _current(times, t)
    i = searchsortedlast(times, t)
    i < firstindex(times) && return firstindex(times)
    i == lastindex(times) && return i
    gap = times[i + 1] - times[i]
    return (times[i + 1] - t) <= gap * _DRIFT ? i + 1 : i
end

# Map a requested time back into the stored range, according to the end-of-series policy.
# Only ever asked about a time the series *is* in span for: `_changevalue` has already stepped
# aside otherwise. In particular a time before the first slice never reaches here, so there is no
# undershoot to clamp — which was once justified by `origin` always defaulting to the first slice,
# and stopped being true the moment an epoch could rewrite `origin` to something earlier.
_inrange(::HoldAtEnd, series::SeriesLayerChange, t) = min(t, last(series.times))
function _inrange(::ErrorAtEnd, series::SeriesLayerChange, t)
    t <= _seriesreach(series) && return t
    return error("this series ran out: the simulation reached " *
                 "$(uconvert(month_mean_duration, t - series.origin)) of elapsed time but the series only " *
                 "covers $(uconvert(month_mean_duration, last(series.times) - series.origin)). Use " *
                 "`atend = HoldAtEnd()` to keep the last slice, `atend = RepeatAtEnd()` to cycle, or a longer " *
                 "series.")
end
# A true modulus, which keeps the overshoot: 13 months into a 12-slice series lands on slice 2, not
# slice 1. One full turn is the stored span **plus** the gap back round to the start, so a 12-slice
# monthly series repeats every 12 months rather than every 11. The spacing was checked to be even at
# attach, so taking the first gap as the step needs no scan here.
# In span, `RevertToLayer` is the identity — out of span it never gets here, because the series has
# already stepped aside and given the layer back.
_inrange(::RevertToLayer, ::SeriesLayerChange, t) = t
function _inrange(::RepeatAtEnd, series::SeriesLayerChange, t)
    times = series.times
    period = last(times) - first(times) + (times[2] - times[1])
    return first(times) + mod(t - first(times), period)
end

# How far a series reaches: its last slice, plus the **whole** gap that slice is in force for — since
# a coordinate is when a slice becomes current, the last one runs for a full step like every other.
#
# A **half** gap here — which is what a nearest-coordinate rule would want — would make a 12-slice
# monthly series cover 11.5 months rather than 12, so a full year of simulation against a twelve-month
# climatology would error a fortnight early under the default `ErrorAtEnd` and never fully reach the
# twelfth slice.
function _seriesreach(series::SeriesLayerChange)
    times = series.times
    length(times) == 1 && return only(times)
    return last(times) + (last(times) - times[end - 1])
end

# `LegacyLoss` is `{NoChange}` but is not a no-op, so it needs a method on its own concrete type —
# more specific than the blanket `{NoChange}` one above, which would otherwise silently swallow it.
function _applychange!(::LegacyLoss, ::AbstractLayer, _, _)
    return error("habitat loss cannot be applied as a layer change: it mutates the " *
                 "ecosystem's supply and abundances and draws at random, so it is neither a " *
                 "pure function of the layer nor reproducible across MPI ranks. Call the " *
                 "deprecated `HabitatLoss(eco, layer, timestep)` directly, or wait for the " *
                 "cell-deactivating intervention that supersedes it.")
end
