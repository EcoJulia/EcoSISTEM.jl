# SPDX-License-Identifier: LGPL-3.0-or-later

module Units

import Unitful
using Unitful: @unit, uconvert

# Month **numbers** - `January == 1` to `December == 12`, a 1-based ordinal matching the `Ti` axis the
# readers build. Re-exported from `Dates` so that naming a month needs no `using Dates` at the call
# site, which matters because `Dates` and this module both export `day`, `week` and `year`: a blanket
# `using Dates` alongside this one makes all three ambiguous, and someone who only wants to say
# "January" should not have to pay that.
#
# **Re-exported, never redefined.** A `const January = 1` of our own would work, since Julia resolves
# a name exported by two modules when the values are identical, but it would duplicate a definition
# `Dates` owns and the duplicate would be free to drift. Re-exporting is the same binding, so
# `using EcoSISTEM.Units, Dates` stays unambiguous for these names by construction.
#
# An ordinal is also calendar-agnostic: every CF calendar numbers its months 1 to 12.
using Dates: Dates
using Dates: January, February, March, April, May, June, July, August,
             September, October, November, December
const day = Unitful.d
const week = Unitful.wk
const year = Unitful.yr
# --- Calendar-month durations -------------------------------------------------------------------
#
# How long a named calendar month lasts - the bridge between the duration world above and the calendar
# world, and what is needed to turn a month's accumulated total into an honest rate.
#
# Not general-purpose duration units: use `day`, `week` or `year` for those. These answer "how long
# did *this slice* accumulate over?", and dividing by the wrong one is a silent error of 1 to 7%.
#
# Named `_duration` and never `_length`, because in this package `Unitful.Length` is the 𝐋 dimension -
# cell sizes, in metres - so a `_length` holding a time would read as a distance.
#
# Exact `Rational`s, for the same reason `arcminute` and `arcsecond` below are: `28.25` in floating
# point would not let the twelve sum to exactly one Julian year, which the tests assert.
@unit january_duration "january_duration" JanuaryDuration 31day false
@unit march_duration "march_duration" MarchDuration 31day false
@unit april_duration "april_duration" AprilDuration 30day false
@unit may_duration "may_duration" MayDuration 31day false
@unit june_duration "june_duration" JuneDuration 30day false
@unit july_duration "july_duration" JulyDuration 31day false
@unit august_duration "august_duration" AugustDuration 31day false
@unit september_duration "september_duration" SeptemberDuration 30day false
@unit october_duration "october_duration" OctoberDuration 31day false
@unit november_duration "november_duration" NovemberDuration 30day false
@unit december_duration "december_duration" DecemberDuration 31day false

# February has three, and deliberately **no bare `february_duration`**: it is the one month whose
# duration is not a single number, so a bare name would be the only one that silently meant a mean.
# Which to use is decided by what the data knows - a real date picks one of the first two, and a
# climatology, where the year is genuinely unknown, the third.
@unit february_common_year_duration "february_common_year_duration" FebruaryCommonYearDuration 28day false
@unit february_leap_year_duration "february_leap_year_duration" FebruaryLeapYearDuration 29day false
@unit february_mean_duration "february_mean_duration" FebruaryMeanDuration (113//4)day false

# The mean month, for a layer whose month is not merely unrecorded but **unknowable** - `bio13` is
# "precipitation of the wettest month", and which month that is varies by cell. Using it is a declared
# approximation, within 7% of any real month, rather than a fudge: the alternative is to pretend to
# know something the data does not record.
#
# The name says what it is. A plain `month` would be a fixed 30.4375 days - the length of no real
# calendar month - and a name that reads as though it were exact is precisely what this scheme exists
# to avoid.
@unit month_mean_duration "month_mean_duration" MonthMeanDuration (487//16)day false

# The mean quarter, for `bio16` to `bio19`. There is deliberately no `q1_duration` to `q4_duration`,
# and the catalogue's own wording is why: these layers are the wettest, driest, warmest and coldest
# quarter, each defined as *any three consecutive months* - a rolling window with twelve possible
# positions, and which one it is varies by cell. A bioclim quarter is therefore not a calendar quarter
# at all, and is unknowable in exactly the way `bio13`'s month is, so only a mean is honest. A genuine
# calendar quarter would just be a sum of the month durations above and needs no unit of its own.
#
# Exactly `3 * month_mean_duration`, so four of them are exactly one Julian year. Both asserted.
@unit quarter_mean_duration "quarter_mean_duration" QuarterMeanDuration (1461//16)day false

# Unitful ships `°` but no subdivisions of it, and every geographic grid this package reads is laid
# out on a whole-arcsecond lattice: WorldClim's 10 arcmin, EarthEnv's and CHELSA's 30 arcsec. Without
# these, such a cell size can only be shown as `0.16666666666666666°` and can only be reasoned about
# as a bare number. Defined as exact `Rational` fractions of a degree so `60arcminute == 1°` holds
# exactly, which `1//60` in floating point would not.
@unit arcminute "′" Arcminute (1 // 60) * Unitful.° false
@unit arcsecond "″" Arcsecond (1 // 3600) * Unitful.° false

# The twelve calendar-month durations in order, so that a slice index can be turned into the interval
# it accumulated over. February is the **mean**, because this indexes a climatology where the year is
# genuinely unknown; a real-dated series picks `february_common_year_duration` or
# `february_leap_year_duration` from its own date instead, which is why this takes a month number and
# cannot be given a year.
const _MONTH_DURATIONS = (january_duration, february_mean_duration,
                          march_duration, april_duration, may_duration,
                          june_duration, july_duration, august_duration,
                          september_duration, october_duration,
                          november_duration, december_duration)

"""
    month_duration(n::Integer)

Return the duration of a calendar month named by its number, using February's **mean** length, since
a bare month number carries no year.

# Arguments

  - `n`: the month, from `January` (1) to `December` (12).
"""
function month_duration(n::Integer)
    1 <= n <= 12 ||
        error("month_duration expects a month number from 1 (January) to 12 (December); got $n.")
    return _MONTH_DURATIONS[n]
end

"""
    month_duration(date::Union{Dates.Date, Dates.DateTime})

Return the duration of the calendar month a date falls in, using February's **real** length - 29 days
in a leap year and 28 otherwise.

This is the form to use whenever the year is known, and the only one that can be right for a dated
series: dividing February 2016's monthly total by February's mean 28.25 days rather than its actual 29
overstates that month's rate by 2.7%, every leap year. The `Integer` method cannot do better, because
a bare month number carries no year, which is why both forms exist.

Typed on the two **Gregorian** date types rather than on `Dates.TimeType`, because every duration it
can return is a Gregorian month length. A date on another calendar is refused rather than answered.

# Arguments

  - `date`: any `Date` or `DateTime`.
"""
function month_duration(date::Union{Dates.Date, Dates.DateTime})
    Dates.month(date) == 2 || return _MONTH_DURATIONS[Dates.month(date)]
    return Dates.isleapyear(date) ? february_leap_year_duration :
           february_common_year_duration
end

# A `Dates.TimeType` that is not one of the two Gregorian types - a `CFTime` calendar, or a bare
# `Dates.Time` with no month at all - is refused, and refusing is the point.
#
# Every alternative would be silently wrong. `Dates.isleapyear` is generic over `TimeType` and applies
# the **Gregorian** rule, which `CFTime` does not override, so a `DateTimeNoLeap` February 2016 would
# come back as 29 days when its calendar says 28, and every `DateTime360Day` month would take a
# Gregorian 30 or 31 when its own is 30. That is not only a February problem: on a 360-day calendar
# `_MONTH_DURATIONS` is wrong for all twelve.
#
# `CFTime` does publish the right answer, `daysinmonth(t)`, and using it would make one implementation
# correct on every calendar. It is deliberately not done here: these functions return a **unit**, a
# `Unitful.FreeUnits`, and there is no unit for "30 days on a 360-day calendar", so a general version
# would have to return a quantity instead and change what every caller receives. That is an API
# decision belonging with the question of whether to support non-Gregorian calendars at all.
function month_duration(date::Dates.TimeType)
    return error("month_duration only knows Gregorian months, so it takes a `Date` or a " *
                 "`DateTime`; got $(typeof(date)). Every duration it can return is a Gregorian " *
                 "month length, and applying one to another calendar is silently wrong rather " *
                 "than approximate - a `DateTimeNoLeap` February would be given 29 days, and " *
                 "every `DateTime360Day` month 30 or 31 instead of 30.")
end

# Anything else is refused rather than coerced. A caller holding real dates that flattened them to
# month numbers somewhere upstream would otherwise get the mean February silently - the exact failure
# the two methods above exist to prevent.
function month_duration(x)
    return error("month_duration takes a month number (1-12, February's mean length) or a date " *
                 "(February's real length, from its year); got $(typeof(x)). If this came from a " *
                 "dated series, pass the dates rather than month numbers - otherwise February is " *
                 "silently divided by its mean 28.25 days instead of its actual 28 or 29.")
end
public month_duration

const localunits = Unitful.basefactors

# Unitful requires a module defining its own units to register them at load, and the conversion
# factors to be merged into its table - a package's units do not survive precompilation otherwise.
function __init__()
    merge!(Unitful.basefactors, localunits)
    return Unitful.register(Units)
end

# The month a monthly slice is being asked for, as a 1-based month number.
#
# A month is an ordinal, so the natural way to name one is `January` to `December`, re-exported from
# `Dates` by this module - `plot(wc, January)` needs no `using Dates`. A duration is accepted too,
# `plot(wc, 1month_mean_duration)`, for a caller who already has one in hand.
#
# **1-based, matching the axis the reader actually builds**: `_mkstackaxis` (`src/datasetread.jl`)
# emits `Ti((1:12) .* month_mean_duration)`, so `1month_mean_duration` is January.
_monthindex(month_number::Integer) = _checkmonth(month_number)
# `uconvert(NoUnits, ...)` rather than a bare division: a caller may hand over a duration in an equal
# but differently named unit, and a bare division leaves that uncancelled, which `Int` then refuses.
function _monthindex(time::Unitful.Time)
    months = uconvert(Unitful.NoUnits, time / month_mean_duration)
    isinteger(months) ||
        error("a monthly layer is named by a whole month; got `$time`, which is $months months.")
    return _checkmonth(Int(months))
end

# Shared bounds check. Spelled with the month names themselves, which say what `1` and `12` mean.
function _checkmonth(month_number::Integer)
    January <= month_number <= December ||
        error("a monthly layer is named by a month from `January` (1) to `December` (12), " *
              "matching the `Ti` axis the reader builds; got `$month_number`.")
    return Int(month_number)
end

# `percent` is Unitful's own and is re-exported, not redefined: `LandmassesAbove(5percent)`
# needs it in scope, and `Unitful.DefaultSymbols` does not carry it.
using Unitful: percent

export day,
       week,
       year,
       percent,
       arcminute,
       arcsecond,
       january_duration,
       february_common_year_duration,
       february_leap_year_duration,
       february_mean_duration,
       march_duration,
       april_duration,
       may_duration,
       june_duration,
       july_duration,
       august_duration,
       september_duration,
       october_duration,
       november_duration,
       december_duration,
       month_mean_duration,
       quarter_mean_duration,
# month numbers (re-exported from Dates; `January == 1`), for naming a slice of a monthly layer
       January,
       February,
       March,
       April,
       May,
       June,
       July,
       August,
       September,
       October,
       November,
       December

end
