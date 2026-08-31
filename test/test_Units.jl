# SPDX-License-Identifier: LGPL-3.0-or-later

module TestUnits

using EcoSISTEM, EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Test
# Named imports, never a blanket `using Dates`: `Dates` and `EcoSISTEM.Units` both export `day`,
# `week` and `year`, which a blanket import would make ambiguous. The module name itself is bound so
# `Dates.Time`/`Dates.TimeType` can be reached without that.
using Dates: Dates, Date, DateTime

@testset "Checking Units" begin
    @test 12month_mean_duration == 1year
end

# Unitful ships `°` but no subdivisions of it, so these exist to describe the whole-arcsecond lattice
# every geographic grid this package reads is laid out on. Defined as exact `Rational` fractions of a
# degree, so the identities below hold with `==` rather than `≈` - which is the point of having them.
@testset "arcminutes and arcseconds subdivide a degree exactly" begin
    @test 60arcminute == 1°
    @test 3600arcsecond == 1°
    @test 60arcsecond == 1arcminute
    # The cell sizes actually in play: WorldClim's 10 arcmin and EarthEnv's/CHELSA's 30 arcsec, whose
    # degree forms (0.1666... and 0.008333...) are exactly what these units exist to avoid printing.
    @test uconvert(arcminute, (1 / 6)°) == 10arcminute
    @test uconvert(arcsecond, (1 / 120)°) == 30arcsecond
    @test uconvert(°, 30arcsecond) == (1 // 120)°
end

# These are how long a *named calendar month* lasts, not general-purpose durations - the bridge
# between the fixed-unit world Unitful requires and the calendar, and what turns a month's
# accumulated total into an honest rate. Getting the wrong one is a silent 1-7% error, so the
# identities below are asserted with `==`: they are defined as exact `Rational`s precisely so that
# they can be.
@testset "calendar-month durations" begin
    # The twelve must sum to exactly one Julian year, with February at its 28.25 d mean. This is
    # the invariant that makes the whole scheme self-consistent: it is what lets a monthly
    # climatology be re-expressed month by month without gaining or losing time over the year.
    months = [january_duration, february_mean_duration, march_duration,
        april_duration, may_duration, june_duration, july_duration,
        august_duration, september_duration, october_duration,
        november_duration, december_duration]
    @test sum(1 .* months) == 1year
    @test uconvert(day, sum(1 .* months)) == (1461 // 4)day     # exact, not 365.25 as a float

    # The three Februaries, and the reason there is no bare `february_duration`: its duration is not
    # one number, so any single name would silently mean a mean.
    @test 1february_common_year_duration == 28day
    @test 1february_leap_year_duration == 29day
    @test uconvert(day, 1february_mean_duration) == (113 // 4)day
    @test uconvert(day, 1february_mean_duration) ==
          (1february_common_year_duration * 3 + 1february_leap_year_duration) /
          4

    # `month_duration` is where the three Februaries earn their keep, and the two forms exist
    # because a bare month number *cannot* choose between them.
    md = EcoSISTEM.Units.month_duration
    @test md(2) == february_mean_duration                    # no year: the mean is all there is
    @test md(Date(2015, 2, 14)) == february_common_year_duration
    @test md(Date(2016, 2, 14)) == february_leap_year_duration
    @test md(DateTime(2016, 2, 1)) == february_leap_year_duration
    # Every other month is unaffected - its length is one number whether or not the year is known.
    for m in (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
        @test md(m) == md(Date(2016, m, 15)) == md(Date(2015, m, 15))
    end
    # The error this prevents: dividing a leap February's monthly total by the mean 28.25 days
    # rather than its real 29 overstates that month's rate by 2.7%, every four years.
    @test uconvert(day, 1md(Date(2016, 2, 1))) / uconvert(day, 1md(2)) ≈
          29 / 28.25
    # ...and anything that is neither is refused rather than coerced, so a caller that has flattened
    # real dates to month numbers upstream finds out instead of silently getting the mean.
    @test_throws ErrorException md("February")
    @test_throws ErrorException md(2.5)
    @test_throws ErrorException md(0)
    @test_throws ErrorException md(13)

    # A `TimeType` that is not Gregorian is refused too, and this is the guard that makes widening
    # `_istimeaxis` to `Dates.TimeType` safe: every duration this returns is a *Gregorian* month
    # length, so answering for another calendar would be silently wrong rather than approximate.
    # `Dates.Time` stands in for the general case here - it is a `TimeType` with no month at all, and
    # needs no extra dependency, where a real `CFTime` calendar would.
    @test Dates.Time(12, 30) isa Dates.TimeType
    @test_throws ErrorException md(Dates.Time(12, 30))
    # ...while the two Gregorian types keep working, which is what the guard must not break.
    @test md(Date(2016, 2, 1)) == md(DateTime(2016, 2, 1)) ==
          february_leap_year_duration

    # The eleven whose duration *is* one number.
    @test 1january_duration == 1march_duration == 31day
    @test 1april_duration == 1june_duration == 30day
    @test 1september_duration == 1november_duration == 30day
    @test 1may_duration == 1july_duration == 1august_duration == 31day
    @test 1october_duration == 1december_duration == 31day

    # The value is pinned directly rather than against a `month` unit, because there is none: as
    # with `quarter` below, `Units.month` does not exist and the equality is unwritable.
    # **Worth knowing, because it is how a real assertion becomes a vacuous one**: a rename sweep
    # over an earlier `1month_mean_duration == 1month` rewrote it to
    # `1month_mean_duration == 1month_mean_duration`, a **tautology that can never fail** - and
    # silently, where the `quarter` case at least errored.
    @test uconvert(day, 1month_mean_duration) == (487 // 16)day
    @test 12month_mean_duration == 1year

    # A mean month is not a real one - the error this whole scheme exists to remove. February is
    # the worst case at 7.2%.
    @test 1february_mean_duration < 1month_mean_duration < 1january_duration
    @test !(1month_mean_duration == 1february_mean_duration)
    @test isapprox(1february_mean_duration / (1month_mean_duration), 0.928,
                   atol = 0.001)

    # The mean quarter, for `bio16`-`bio19`. There is deliberately no `q1_duration`...`q4_duration`:
    # a bioclim "quarter" is *any 3 consecutive months* - a rolling window whose position varies by
    # cell - so it is unknowable in the way `bio13`'s month is, and only a mean is honest.
    # Pinned against `3month_mean_duration` rather than against a `quarter` unit, because there is
    # none: the identity below fixes the same value without naming a unit that does not exist.
    @test 1quarter_mean_duration == 3month_mean_duration
    @test uconvert(day, 1quarter_mean_duration) == (1461 // 16)day
    @test 4quarter_mean_duration == 1year
    # Real calendar quarters are 90.25, 91, 92 and 92 d, so the mean is the length of none of them -
    # which is exactly what the name admits.
    @test 1quarter_mean_duration >
          1january_duration + 1february_mean_duration + 1march_duration
    @test 1quarter_mean_duration <
          1july_duration + 1august_duration + 1september_duration

    # They parse from a string, which is what lets a catalogue cell name one.
    @test uparse("january_duration",
                 unit_context = [Unitful, EcoSISTEM.Units]) ==
          EcoSISTEM.Units.january_duration
end

# **`_monthindex` lives here now, not in a plot extension.** It was hoisted into `Units` when the
# monthly-climate recipes were split across two extensions and both needed it - and `Units` is where
# it belongs anyway, beside `month_mean_duration` and the `January`...`December` re-exports.
# **1-based, matching the axis the readers actually build.** The predecessor was 0-based and
# labelled every month one late; only a hand-built fixture no reader produces kept that hidden.
@testset "_monthindex" begin
    @test EcoSISTEM.Units._monthindex(January) == 1
    @test EcoSISTEM.Units._monthindex(December) == 12
    # A duration is accepted too, including one in an equal-but-differently-named unit - which is
    # what the `uconvert` in `_monthindex` is for.
    @test EcoSISTEM.Units._monthindex(1month_mean_duration) == 1
    # Out of range and fractional are refused, named in months rather than the old `error("NO")`.
    @test_throws ErrorException EcoSISTEM.Units._monthindex(0)
    @test_throws ErrorException EcoSISTEM.Units._monthindex(13)
    @test_throws ErrorException EcoSISTEM.Units._monthindex(1.5month_mean_duration)
end

end
