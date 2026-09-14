# SPDX-License-Identifier: LGPL-3.0-or-later
#
# *When* an intervention acts - always in elapsed simulation time, never in step counts, so the
# answer is the same at any timestep.

using Unitful
using Dates: Dates

"""
    AbstractSchedule

**When** an [`Intervention`](@ref) fires - [`EveryStep`](@ref), [`AtTime`](@ref), [`AtTimes`](@ref),
[`BetweenTimes`](@ref), [`EveryInterval`](@ref), [`AtDates`](@ref), [`EveryYear`](@ref) or
[`NeverScheduled`](@ref).

A type rather than a predicate function, so that a schedule can be reported and checked rather than
merely called, and so each rule is a method instead of a branch retaken every step.
"""
abstract type AbstractSchedule end

"""    EveryStep() <: AbstractSchedule - fires after every timestep, and not at the start. """
struct EveryStep <: AbstractSchedule end

"""
    NeverScheduled() <: AbstractSchedule

Never fires. For disabling an intervention without removing it from a set, so a configuration can
keep its shape while one part of it is turned off.
"""
struct NeverScheduled <: AbstractSchedule end

"""
    AtTime(time::Unitful.Time)

Fires on the single step that *reaches* `time` - the first step whose elapsed time is at or past it.
A `time` at or before the start of the run acts on the starting state, before the first step.

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

Fires once for each of `times`, on the step that reaches each - [`AtTime`](@ref) repeated.

# Arguments

  - `times`: the elapsed simulation times to fire at.
"""
struct AtTimes{T} <: AbstractSchedule
    times::T
end

"""
    BetweenTimes(from::Unitful.Time, to::Unitful.Time)

Fires on every step whose elapsed time lies in `[from, to]`. It acts over steps taken, so not at the
start of the run even when `from` is zero.

# Arguments

  - `from`, `to`: the inclusive elapsed-time bounds.
"""
struct BetweenTimes{F <: Unitful.Time, T <: Unitful.Time} <: AbstractSchedule
    from::F
    to::T
end

"""
    EveryInterval(interval::Unitful.Time)

Fires at every whole multiple of `interval` - `0`, `interval`, `2interval` and on - each on the step
that reaches it, exactly as [`AtTimes`](@ref) would over that unending list. The multiples are
counted from the start of the run, so the answer does not depend on the timestep; the multiple at zero
acts on the starting state, before the first step.

# Arguments

  - `interval`: the elapsed simulation time between firings, which must be positive.
"""
struct EveryInterval{T <: Unitful.Time} <: AbstractSchedule
    interval::T

    function EveryInterval(interval::Unitful.Time)
        interval > zero(interval) ||
            throw(ArgumentError("an `EveryInterval` needs a positive interval, but got " *
                                "$interval: every multiple of it would fall on the same instant."))
        return new{typeof(interval)}(interval)
    end
end

"""
    AtDates(dates)

Fires once for each of `dates`, on the step that reaches each - [`AtTimes`](@ref) given real dates
rather than elapsed times. A date becomes an elapsed time through the run's epoch and calendar (see
[`build_ecosystem`](@ref)), so a run with no epoch refuses one before its first step. A date at or
before the epoch acts on the starting state, before the first step, as an elapsed time at or before
zero does. A date the run reaches must fall at the end of a step, or the step be no longer than a
day: one falling inside a longer step would be acted on late, so the run is refused before its first
step.

# Arguments

  - `dates`: the `Date`s or `DateTime`s to fire at.
"""
struct AtDates{D <: AbstractVector{<:Dates.TimeType}} <: AbstractSchedule
    dates::D
end

"""
    EveryYear(; month = 1, day = 1)

Fires on the step that reaches `day` of `month` in every year of the run, from the epoch on - an
anniversary on the epoch acts on the starting state, and one before it does not fire. Its dates are placed through the run's epoch and
calendar as [`AtDates`](@ref)'s are, so a run with no epoch refuses it before its first step, and
so is each anniversary checked against the steps as a date is: under `ExactDates` 1 January is
365 or 366 days after the one before, so no step of `month_mean_duration` or `year` ends on every
one, while under `MeanMonths` every twelfth `month_mean_duration` does.

# Arguments

  - `month`: the month of the year, from 1 to 12.
  - `day`: the day of that month, which must be one every year has - so 29 February is refused.
"""
struct EveryYear <: AbstractSchedule
    month::Int
    day::Int

    function EveryYear(; month::Integer = 1, day::Integer = 1)
        1 <= month <= 12 ||
            throw(ArgumentError("`month` must be from 1 to 12, but got $month."))
        days = Dates.daysinmonth(2001, month)
        1 <= day <= days ||
            throw(ArgumentError("`day` must be one every year has, but " *
                                "$(Dates.monthname(month)) has $days days in a common year " *
                                "and got $day."))
        return new(month, day)
    end
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# These are small, but they nest inside an `Intervention`, and a parametric struct's default `show`
# prints its full type signature - `AtTime{Quantity{Float64, 𝐓, Unitful.FreeUnits{(yr,), 𝐓,
# nothing}}}(5.0 yr)` for what a caller wrote as `AtTime(5.0year)`. The one-liner is that call.
#
# The fieldless schedules need nothing: `EveryStep()` and `NeverScheduled()` already print exactly
# as they are written.
Base.show(io::IO, s::AtTime) = print(io, "AtTime($(s.time))")
Base.show(io::IO, s::AtTimes) = print(io, "AtTimes($(s.times))")
function Base.show(io::IO, s::BetweenTimes)
    return print(io, "BetweenTimes($(s.from), $(s.to))")
end
Base.show(io::IO, s::EveryInterval) = print(io, "EveryInterval($(s.interval))")
Base.show(io::IO, s::AtDates) = print(io, "AtDates(", s.dates, ")")
function Base.show(io::IO, s::EveryYear)
    return print(io, "EveryYear(month = $(s.month), day = $(s.day))")
end
