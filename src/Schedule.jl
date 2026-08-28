# SPDX-License-Identifier: LGPL-3.0-or-later
#
# *When* an intervention acts — always in elapsed simulation time, never in step counts, so the
# answer is the same at any timestep.

using Unitful

"""
    AbstractSchedule

**When** an [`Intervention`](@ref) fires — [`EveryStep`](@ref), [`AtTime`](@ref), [`AtTimes`](@ref),
[`BetweenTimes`](@ref) or [`NeverScheduled`](@ref).

A type rather than a predicate function, so that a schedule can be reported and checked rather than
merely called, and so each rule is a method instead of a branch retaken every step.
"""
abstract type AbstractSchedule end

"""    EveryStep() <: AbstractSchedule — fires on every timestep. """
struct EveryStep <: AbstractSchedule end

"""
    NeverScheduled() <: AbstractSchedule

Never fires. For disabling an intervention without removing it from a set, so a configuration can
keep its shape while one part of it is turned off.
"""
struct NeverScheduled <: AbstractSchedule end

"""
    AtTime(time::Unitful.Time)

Fires on the single step that *reaches* `time` — the first step whose elapsed time is at or past it.

**Reaches, not equals**: elapsed time accumulates as a float and a run's steps need not land on
`time` exactly, so an equality test would silently never fire. `_current` solves the same problem the
same way for series slices.

# Arguments

  - `time`: the elapsed simulation time to fire at.
"""
struct AtTime{T <: Unitful.Time} <: AbstractSchedule
    time::T
end

"""
    AtTimes(times)

Fires once for each of `times`, on the step that reaches each — [`AtTime`](@ref) repeated.

# Arguments

  - `times`: the elapsed simulation times to fire at.
"""
struct AtTimes{T} <: AbstractSchedule
    times::T
end

"""
    BetweenTimes(from::Unitful.Time, to::Unitful.Time)

Fires on every step whose elapsed time lies in `[from, to]`.

# Arguments

  - `from`, `to`: the inclusive elapsed-time bounds.
"""
struct BetweenTimes{F <: Unitful.Time, T <: Unitful.Time} <: AbstractSchedule
    from::F
    to::T
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# These are small, but they nest inside an `Intervention`, and a parametric struct's default `show`
# prints its full type signature — `AtTime{Quantity{Float64, 𝐓, Unitful.FreeUnits{(yr,), 𝐓,
# nothing}}}(5.0 yr)` for what a caller wrote as `AtTime(5.0year)`. The one-liner is that call.
#
# The fieldless schedules need nothing: `EveryStep()` and `NeverScheduled()` already print exactly
# as they are written.
Base.show(io::IO, s::AtTime) = print(io, "AtTime($(s.time))")
Base.show(io::IO, s::AtTimes) = print(io, "AtTimes($(s.times))")
function Base.show(io::IO, s::BetweenTimes)
    return print(io, "BetweenTimes($(s.from), $(s.to))")
end
