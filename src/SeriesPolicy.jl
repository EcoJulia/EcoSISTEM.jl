# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a stored series does past its last slice, and what its time coordinates mean.

using Unitful
using Dates: Dates

"""
    AbstractSeriesEnd

What a [`SeriesLayerChange`](@ref) does once elapsed simulation time runs past its last stored
slice - [`ErrorAtEnd`](@ref), [`HoldAtEnd`](@ref), [`RepeatAtEnd`](@ref) or [`RevertToLayer`](@ref).
Written as the `atend` keyword of [`SeriesChange`](@ref), and a type rather than a name drawn from a
list, so an unrecognised policy is refused by the signature where it is written and the behaviours
are separate methods rather than a branch retaken every step.

There is no matching `atstart`, because before its first slice a series has only one sensible
reading: it has not started, so it says nothing and the layer stands. A policy with one option is not
a policy.
"""
abstract type AbstractSeriesEnd end

"""    ErrorAtEnd <: AbstractSeriesEnd - running past the last slice is an error. The default. """
struct ErrorAtEnd <: AbstractSeriesEnd end

"""    HoldAtEnd <: AbstractSeriesEnd - the last slice persists for the rest of the simulation. """
struct HoldAtEnd <: AbstractSeriesEnd end

"""    RepeatAtEnd <: AbstractSeriesEnd - the series cycles, by a true modulus of its own period. """
struct RepeatAtEnd <: AbstractSeriesEnd end

"""
    RevertToLayer <: AbstractSeriesEnd

Past its last slice the series stops driving the layer, which returns to the values it held before
the series was attached - its own spec's, for a [`ReplaceWith`](@ref) series.

This is the not-yet-started rule of [`AbstractSeriesEnd`](@ref) applied at the far end, where the
alternatives are genuine. It is also the only end policy for which the layer's own values are still
needed: the other three keep the series in charge.
"""
struct RevertToLayer <: AbstractSeriesEnd end

"""
    AbstractSeriesCalendar

What kind of time coordinate a [`SeriesLayerChange`](@ref)'s slices carry, and so what a simulation
epoch is able to do with them - [`DatedSeries`](@ref), [`MonthOfYearSeries`](@ref) or
[`UndatedSeries`](@ref). Written as the `calendar` keyword of [`SeriesChange`](@ref).

**A series has to say which kind it is; it cannot be inferred from the values.** Coordinates of one,
two and three months read equally well as *the first three months of the year* and as *one, two and
three months into my experiment*, and phase-locking the second to January would be silently wrong.
Only the two unambiguous cases are inferred: real dates in the lookup give a `DatedSeries`, and no
time lookup at all gives an `UndatedSeries`. A climatology must say `calendar = MonthOfYearSeries()`.

The set is extensible: a series subdivided some other way (day-of-year, say) is a new subtype
supplying `_calendarorigin`, not a change to these.
"""
abstract type AbstractSeriesCalendar end

"""
    DatedSeries(start::Dates.TimeType)

Slices carry real dates, `start` being the first slice's. The epoch places the series directly: the
slice current at elapsed zero is the one covering the epoch's own date.

Inferred from a source whose `Ti` lookup holds `Dates` values, so it is rarely written by hand.
"""
struct DatedSeries <: AbstractSeriesCalendar
    start::Dates.TimeType
end

"""
    MonthOfYearSeries()

Slices are months of the year with the year unknown - a monthly climatology. The epoch selects the
slice matching its **own calendar month**, so a run beginning in July starts on the July slice
rather than on slice one.

Matching is by **month number**, never by elapsed duration into the year. The slices are spaced by
`month_mean_duration` (30.44 d), so six of them reach day 182.6 while the real 1 July is day 181 -
proportional matching would land a July epoch on the June slice.
"""
struct MonthOfYearSeries <: AbstractSeriesCalendar end

"""
    UndatedSeries()

No calendar identity at all - a synthetic stack, or any series whose coordinates are simply offsets.
There is nothing for an epoch to bind to, so elapsed zero is the first slice unless the series'
`origin` says otherwise. This is the only calendar for which `origin` is meaningful, and the default
for any source that does not carry real dates.
"""
struct UndatedSeries <: AbstractSeriesCalendar end
