# SPDX-License-Identifier: LGPL-3.0-or-later
#
# How a layer changes in time, as a caller declares it — recipes, before anything is built.

using Unitful
using Dates: Dates

"""
    AbstractChangeSpec

Abstract supertype of the change *recipes* — a declared change that has not yet met a layer, and
so has not yet been unit-checked. Materialised into an [`AbstractLayerChange`](@ref) at attach.
"""
abstract type AbstractChangeSpec end

"""
    ReplaceWith(shape)

Declare that `shape`'s values are absolute: the layer's value *is* `shape`, in the layer's own unit.

# Arguments

  - `shape`: something that genuinely differs from step to step — a [`SeriesChange`](@ref), or a
    [`PatternedChange`](@ref) whose shape returns the full value. A bare constant is rejected,
    because writing the same value every step is a one-off operation on the ecosystem rather than a
    change of the layer.
"""
struct ReplaceWith{S} <: AbstractChangeSpec
    shape::S
end

"""
    OffsetBy(shape)

Declare that `shape`'s values are offsets from the layer's values as they stand when the change is
attached — an *interval*, so `K` rather than `°C` on a temperature layer.

# Arguments

  - `shape`: as [`ReplaceWith`](@ref), something that differs from step to step. A bare constant is
    rejected: a fixed offset applied every step is idempotent after the first, so it too is a
    one-off operation rather than a change.
"""
struct OffsetBy{S} <: AbstractChangeSpec
    shape::S
end

"""
    IncrementBy(shape)

Declare that `shape`'s values are rates: they are accumulated into the layer, in its unit per unit
time.

# Arguments

  - `shape`: a rate, in the layer's unit per unit time. This is the one recipe that accepts a bare
    constant, which gives a steady drift.
"""
struct IncrementBy{S} <: AbstractChangeSpec
    shape::S
end

"""
    PatternedChange(amplitude, timescale; shape = sinusoidal)

A shape of elapsed simulation time: `amplitude * shape(elapsed / timescale)`.

# Arguments

  - `amplitude`: what `shape`'s output is scaled by, and where the unit comes from.
  - `timescale`: the elapsed time that maps to one unit of `shape`'s argument. The phase is not
    wrapped, so this is a period for a periodic shape and a transition width for a monotone one.
  - `shape`: any function of a dimensionless phase — a sigmoid, a ramp, a step. Defaults to
    [`sinusoidal`](@ref), one full cycle per `timescale`.

Like [`SeriesChange`](@ref) this is a shape rather than a change in itself: the recipe wrapping it
says whether its values are absolute, an offset or a rate.
"""
struct PatternedChange{F, V, P}
    shape::F
    amplitude::V
    timescale::P

    # The sole constructor, so there is one way to build this and it is the documented one. It also
    # has to be the inner one: the fields are ordered `shape` first (it is the type's defining
    # parameter) while the call takes it last, as a keyword — leaving the generated positional
    # constructor in place would expose that mismatch as a second, silently different spelling.
    function PatternedChange(amplitude, timescale; shape = sinusoidal)
        return new{typeof(shape), typeof(amplitude), typeof(timescale)}(shape,
                                                                        amplitude,
                                                                        timescale)
    end
end

"""
    SeriesChange(source; times = nothing, origin = nothing, atend = ErrorAtEnd(),
                 calendar = nothing)

A stack of stored slices to be indexed by elapsed simulation time — a read climate series, or any
`(Y, X, Ti)` array of values.

# Arguments

  - `source`: the slices. Anything indexable by a time dimension.
  - `times`: one time per slice. Taken from `source`'s own `Ti` lookup by default, or assumed monthly
    for an array carrying no lookup at all.
  - `origin`: the point in the series' **own coordinate** that elapsed time zero corresponds to.
    Defaults to the first stored time, so a series starts at its own beginning whatever its axis is
    anchored to.
  - `atend`: what happens past the last slice — an [`AbstractSeriesEnd`](@ref), one of
    [`ErrorAtEnd`](@ref) (the default), [`HoldAtEnd`](@ref) or [`RepeatAtEnd`](@ref).
  - `calendar`: what the coordinates *mean* — an [`AbstractSeriesCalendar`](@ref). A source whose
    lookup holds real dates is inferred to be a [`DatedSeries`](@ref) and anything else an
    [`UndatedSeries`](@ref). A monthly climatology has to say `calendar = MonthOfYearSeries()` for a
    run's epoch to start it in the right month, because month-of-year coordinates cannot be told
    from plain offsets by their values alone.

**`origin = 0` is not the same as leaving `origin` out**, and the difference is a real one rather
than a technicality. Omitting it means *start at the first slice*; zero means *the coordinate zero*,
which on an ordinary monthly climatology — slices at 1 to 12 months — is a month **before** the first
slice, and a series before its first slice says nothing, so the layer holds its own values for that
month. Write the lead-in when you want one; leave `origin` out when you mean "from the beginning".

`origin` is accepted only for an [`UndatedSeries`](@ref), the one case where nothing else can set the
zero point: for the other two the epoch fixes the phase, and a second knob for the same thing could
only contradict it. To say that an undated series really begins at a particular date, give it real
`times` rather than an `origin` — that states what every slice is, rather than only where zero sits.

Like [`PatternedChange`](@ref) this is a shape, so the mode is named by the recipe wrapping it:
[`ReplaceWith`](@ref) for stored values the layer takes on, [`OffsetBy`](@ref) for a stored deviation
from the layer's own values, [`IncrementBy`](@ref) for a stored rate.
"""
struct SeriesChange{S, T, O, E <: AbstractSeriesEnd, C}
    source::S
    times::T
    origin::O
    atend::E
    calendar::C

    # The sole constructor: everything but `source` is optional and named, so a positional spelling
    # would be an unlabelled `(source, times, origin, atend, calendar)` that no caller should be
    # writing.
    function SeriesChange(source; times = nothing, origin = nothing,
                          atend::AbstractSeriesEnd = ErrorAtEnd(),
                          calendar::Union{Nothing, AbstractSeriesCalendar} = nothing)
        return new{typeof(source), typeof(times), typeof(origin),
                   typeof(atend), typeof(calendar)}(source,
                                                    times,
                                                    origin,
                                                    atend,
                                                    calendar)
    end
end

"""
    spec1 + spec2
    EcoSISTEM.CombinedChange(specs...)

Declare that a layer carries several changes at once, added together — a stored monthly series plus
a multi-year trend, say, where the trend offsets the whole seasonal pattern.

**Write it as a sum.** That is the spelling this exists for, and the only one most callers need:

```julia
ReplaceWith(SeriesChange(monthly)) + IncrementBy(0.02K / year)
```

`CombinedChange` is the type a sum builds, and is supported (though unexported, so it needs
qualifying) for the case a `+` chain cannot express — combining a collection whose length is not
known when the code is written:

```julia
EcoSISTEM.CombinedChange(specs...)     # rather than `reduce(+, specs)`
```

Either way the parts are summed as functions of elapsed time, not applied one after another; see
[`SumOfLayerChanges`](@ref) for which combinations are meaningful.
"""
struct CombinedChange{S <: Tuple} <: AbstractChangeSpec
    specs::S

    function CombinedChange(specs::AbstractChangeSpec...)
        flat = _flattenspecs(specs)
        length(flat) > 1 ||
            error("`CombinedChange` adds several changes together, so it needs more than one; " *
                  "a single change needs no combining.")
        return new{typeof(flat)}(flat)
    end
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# The same rule as the layer specs: the one-liner is the expression that builds it. These nest —
# a change spec sits inside a `Varying`, which sits inside a habitat's keyword — so the compact form
# is the one that matters, and it is why `CombinedChange` prints as the `+` that spells it.
function Base.show(io::IO, c::ReplaceWith)
    return print(io, "ReplaceWith(", sprint(show, c.shape), ")")
end

function Base.show(io::IO, c::OffsetBy)
    return print(io, "OffsetBy(", sprint(show, c.shape), ")")
end

function Base.show(io::IO, c::IncrementBy)
    return print(io, "IncrementBy(", sprint(show, c.shape), ")")
end

# `shape` is named only when it is not the default, which keeps the common case to the two arguments
# a caller actually wrote.
function Base.show(io::IO, c::PatternedChange)
    shape = c.shape === sinusoidal ? "" : ", shape = $(nameof(c.shape))"
    return print(io, "PatternedChange($(c.amplitude), $(c.timescale)$(shape))")
end

# The source is a stack of slices, so it is reported by size rather than printed; everything else
# appears only where it was set.
function Base.show(io::IO, c::SeriesChange)
    parts = String[]
    isnothing(c.times) || push!(parts, "times = $(length(c.times)) slices")
    isnothing(c.origin) || push!(parts, "origin = $(c.origin)")
    c.atend isa ErrorAtEnd ||
        push!(parts, "atend = $(nameof(typeof(c.atend)))()")
    isnothing(c.calendar) ||
        push!(parts, "calendar = $(nameof(typeof(c.calendar)))()")
    tail = isempty(parts) ? "" : ", " * join(parts, ", ")
    return print(io, "SeriesChange(<source>$(tail))")
end

# Printed as the `A + B` that builds it, because that spelling is the only one a caller writes.
function Base.show(io::IO, c::CombinedChange)
    return print(io, join(map(x -> sprint(show, x), c.specs), " + "))
end
