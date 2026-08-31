# SPDX-License-Identifier: LGPL-3.0-or-later

module TestLayerChange

using EcoSISTEM
# `[C7-VIS]` C: these are `public` rather than exported - a spec is what a user writes,
# and these are what it materialises into.
using EcoSISTEM: AbsoluteChange, NoChange, PatternedLayerChange, RateChange,
                 RelativeChange, SeriesLayerChange, SteadyLayerChange,
                 SumOfLayerChanges
# Abstract types are `public`, not exported, so `using EcoSISTEM` does not bring them in.
using EcoSISTEM: AbstractLayer
using Distributions
using Test
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using RasterDataSources
using DimensionalData
using DimensionalData: DimArray, X, Y, Ti
using DimensionalData.Lookups: NoLookup
import Dates

include("TestCases.jl")

# A bare continuous regime on a declared axis, for exercising the unit contract directly.
function _testregime(value, axis)
    layer = ContinuousRegime(fill(value, 5, 5), 1.0km, NoLayerChange())
    return EcoSISTEM._reaxis(layer, axis)
end

# A species list aligned to a **one-layer** habitat whose regime is a temperature: a `Temperature`
# tolerance and a single demand on `axis`.
#
# **Needed because neither shared fixture fits, which had never been noticed.** `Test1Ecosystem`
# carries `Int64` niche codes (it is built from `simplenichehabitat`) and `TestMultiEcosystem` carries
# *two* demands, so against a one-layer `K` habitat the first fails `_checkaligned` on the tolerance's
# type and the second on the layer count. Both assertions below sat behind this file's `AxisArray`
# fixture and so had never run.
function _alignedspecies(demand)
    n = 3
    abun = rand(Multinomial(100, n))
    tolerance = NicheTolerance(Temperature, Normal, fill(1.0K, n),
                               fill(0.1K, n))
    movement = BirthOnlyMovement(GaussianKernel.(fill(1.0km, n), 10.0e-4))
    # `longevity`/`survival`/`boost` are plain `Float64`s; only the two rates carry a unit.
    param = EqualPop(0.1 / month_mean_duration, 0.1 / month_mean_duration, 1.0,
                     0.1, 0.1)
    return SpeciesList(n, tolerance, abun, demand, movement, param,
                       fill(true, n))
end

@testset "Condition update" begin
    eco = Test1Ecosystem()
    @test_nowarn EcoSISTEM.regimeupdate!(eco, 1month_mean_duration)
    @test_nowarn EcoSISTEM.supplyupdate!(eco, 1month_mean_duration)

    eco = TestMultiEcosystem()
    @test_nowarn EcoSISTEM.regimeupdate!(eco, 1month_mean_duration)
    @test_nowarn EcoSISTEM.supplyupdate!(eco, 1month_mean_duration)

    # The two ecosystems above are built from `NicheSpec`/`UniformSpec`, so their layers carry a
    # `SteadyLayerChange`. This third one is the only place either updater is driven on a layer that
    # actually holds a **series** - which is the whole of what it is here for. What the series itself
    # does with elapsed time is `_seriesindex`'s business and is covered at length below, so none of
    # it is re-asserted through the ecosystem.
    # `Varying(spec, ReplaceWith(SeriesChange(stack)))` is how a spec declares a series, so this
    # needs no raster and no deprecated builder. A supply spec is stated **per area** and
    # multiplied by the cell area (0.25 km^2 here), so 40 kJ/km^2/day lands as the 10 kJ/day per cell
    # the fixture it replaces built by hand; the series slices are already per-cell.
    area = StudyArea(extent = (5.0km, 5.0km), cellsize = 0.5km,
                     verbosity = :silent)
    tempstack = cat((fill(i * 1.0K, 10, 10) for i in 1:12)..., dims = 3)
    solarstack = cat((fill(i * 10.0kJ / day, 10, 10) for i in 1:3)..., dims = 3)
    varying = GridHabitat(regime = Varying(UniformSpec(1.0K,
                                                       axis = Temperature),
                                           ReplaceWith(SeriesChange(tempstack))),
                          supply = Varying(UniformSpec(40.0kJ / km^2 /
                                                       day,
                                                       axis = SolarRadiation),
                                           ReplaceWith(SeriesChange(solarstack))),
                          area = area)
    # The guard that this fixture is still doing its job: both halves must really be series-backed,
    # or the two updates below would be exercising the steady path the cases above already cover.
    @test varying.regime.change isa SeriesLayerChange
    @test varying.supply.change isa SeriesLayerChange
    # The species list is incidental to what is under test - a series-backed layer advancing by
    # elapsed time - so it only has to align with this habitat: one `SolarRadiation` demand.
    eco = Ecosystem(_alignedspecies(Demand{SolarRadiation}(fill(2.0kJ / day, 3))),
                    varying,
                    NicheSuitability{EcoSISTEM.axisof(varying.regime),
                                     eltype(varying.regime)}())
    @test_nowarn EcoSISTEM.regimeupdate!(eco, 1month_mean_duration)
    @test_nowarn EcoSISTEM.supplyupdate!(eco, 1month_mean_duration)
end

@testset "The unit contract" begin
    # A change's unit comes from the layer it is attached to, never from the change itself. An
    # absolute change is in the layer's own unit; an increment is in that unit per unit time.
    temp = _testregime(274.0K, Temperature)
    @test isnothing(EcoSISTEM.changeunit(NoChange(), temp))
    @test EcoSISTEM.changeunit(AbsoluteChange(), temp) == K
    @test EcoSISTEM.changeunit(RelativeChange(), temp) == K
    @test EcoSISTEM.changeunit(RateChange(), temp) == K / s

    # Precipitation is the case this replacement exists to fix. Its canonical unit is a *rate*
    # (`mm/d`), so a change function hard-coding `mm` - a depth - adds a depth to a rate and throws
    # a bare `DimensionError` even at rate zero.
    rain = _testregime(3.0mm / day, Precipitation)
    @test EcoSISTEM.changeunit(AbsoluteChange(), rain) == mm / day
    @test EcoSISTEM.changeunit(RateChange(), rain) == mm / day / s

    # A whole year of drift accumulates exactly the declared amount.
    EcoSISTEM.setchange!(rain, IncrementBy(1.0mm / day / year))
    EcoSISTEM._layerupdate!(rain, 1.0year, 1.0year)
    @test all(≈(4.0mm / day), rain.matrix)

    # ...and a rate that is a depth per time, not a rate per time, is rejected with a message naming
    # both units rather than a bare DimensionError.
    err = try
        EcoSISTEM.setchange!(rain, IncrementBy(1.0mm / year))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("Precipitation", err.msg)
    @test occursin("mm", err.msg)

    # The affine case. A relative or rate change is an *interval*, so it is in K even where the
    # layer holds positions in °C - the layer's unit is absolutised first. Getting this wrong is not
    # merely imprecise: Unitful rejects `10.0°C + 1.0°C` outright with an `AffineError`, so a
    # relative change read in °C would not run at all.
    celsius = _testregime(10.0u"°C", Temperature)
    @test EcoSISTEM.changeunit(AbsoluteChange(), celsius) == u"°C"
    @test EcoSISTEM.changeunit(RelativeChange(), celsius) == K
    @test EcoSISTEM.changeunit(RateChange(), celsius) == K / s
    @test_nowarn EcoSISTEM.setchange!(celsius, IncrementBy(1.0K / year))

    # ...and it applies without throwing, which is the regression this actually guards.
    EcoSISTEM.setchange!(celsius, OffsetBy(PatternedChange(2.0K, 1.0year)))
    @test_nowarn EcoSISTEM._layerupdate!(celsius, 0.25year, 0.25year)
    @test all(≈(12.0u"°C"), celsius.matrix)
end

@testset "A constant is a change only as a rate" begin
    # Writing the same values every step is idempotent after the first, so it is a one-off
    # operation on the ecosystem rather than a change of the layer. Only a *rate* accumulates.
    temp = _testregime(274.0K, Temperature)
    @test_nowarn EcoSISTEM.setchange!(temp, IncrementBy(2.0K / year))
    @test temp.change isa SteadyLayerChange
    for spec in (ReplaceWith(300.0K), OffsetBy(2.0K))
        err = try
            EcoSISTEM.setchange!(temp, spec)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("one-off operation", err.msg)
    end
end

@testset "Applying a change" begin
    # A steady rate accumulates rate × timestep.
    temp = _testregime(274.0K, Temperature)
    EcoSISTEM.setchange!(temp, IncrementBy(2.0K / year))
    EcoSISTEM._layerupdate!(temp, 1.0year, 1.0year)
    @test all(≈(276.0K), temp.matrix)

    # A pattern is a function of elapsed time, not of the layer's own values - so two layers at the
    # same elapsed time agree however they got there. (The change function this replaces fed the
    # matrix back into `sin`, making it a path-dependent walk.)
    onestep = _testregime(274.0K, Temperature)
    EcoSISTEM.setchange!(onestep, OffsetBy(PatternedChange(5.0K, 1.0year)))
    EcoSISTEM._layerupdate!(onestep, 0.25year, 0.25year)

    manysteps = _testregime(274.0K, Temperature)
    EcoSISTEM.setchange!(manysteps, OffsetBy(PatternedChange(5.0K, 1.0year)))
    for i in 1:5
        EcoSISTEM._layerupdate!(manysteps, i * 0.05year, 0.05year)
    end
    @test all(onestep.matrix .≈ manysteps.matrix)
    # A quarter of the way round the default sinusoid is its peak.
    @test all(≈(279.0K), onestep.matrix)

    # `AbsoluteChange` - the layer's value *is* the pattern, with no baseline. Reachable with a
    # shape that returns the whole value rather than a deviation from one.
    absolute = _testregime(274.0K, Temperature)
    EcoSISTEM.setchange!(absolute,
                         ReplaceWith(PatternedChange(1.0K, 1.0year,
                                                     shape = x -> 280 +
                                                                  10sinpi(2x))))
    EcoSISTEM._layerupdate!(absolute, 0.25year, 0.25year)
    @test all(≈(290.0K), absolute.matrix)
    EcoSISTEM._layerupdate!(absolute, 0.75year, 0.5year)
    @test all(≈(270.0K), absolute.matrix)

    # An arbitrary shape need not be periodic: the phase handed to it is deliberately unwrapped, so
    # a sigmoid stays monotone across many timescales instead of resetting each one.
    sigmoid = _testregime(274.0K, Temperature)
    EcoSISTEM.setchange!(sigmoid,
                         OffsetBy(PatternedChange(10.0K, 1.0year,
                                                  shape = x -> 1 /
                                                               (1 + exp(-x)))))
    values = map(1:6) do i
        EcoSISTEM._layerupdate!(sigmoid, i * 1.0year, 1.0year)
        return first(sigmoid.matrix)
    end
    @test issorted(values)
    @test values[end] > values[1]

    # A collection applies each sub-layer's own change.
    warming = _testregime(274.0K, Temperature)
    static = _testregime(274.0K, Temperature)
    EcoSISTEM.setchange!(warming, IncrementBy(2.0K / year))
    EcoSISTEM._layerupdate!(LayerCollection((; warming, static)), 1.0year,
                            1.0year)
    @test all(≈(276.0K), warming.matrix)
    @test all(≈(274.0K), static.matrix)
end

# A stack of 12 monthly slices, slice `i` holding `base + i * step` everywhere.
function _teststack(base, step, dim = (5, 5), n = 12)
    return cat((fill(base + i * step, dim) for i in 1:n)..., dims = 3)
end

@testset "A shape is valid under every mode" begin
    # Regression: these two dispatched on the shape being *anything*, so under `RateChange` they
    # were ambiguous with every shape-typed method - one narrowing the shape, the other the mode,
    # neither more specific. `IncrementBy(PatternedChange(...))`, an oscillating rate, therefore
    # always failed despite being documented as valid.
    layer = _testregime(274.0K, Temperature)
    @test EcoSISTEM._attachchange(IncrementBy(PatternedChange(1.0K / year,
                                                              1.0year)),
                                  layer) isa
          PatternedLayerChange{RateChange}
    @test EcoSISTEM._attachchange(OffsetBy(PatternedChange(1.0K, 1.0year)),
                                  layer) isa
          PatternedLayerChange{RelativeChange}
    @test EcoSISTEM._attachchange(IncrementBy(SeriesChange(_teststack(0.0K /
                                                                      year,
                                                                      1.0K /
                                                                      year))),
                                  layer) isa SeriesLayerChange{RateChange}

    # A constant still means what it meant: a drift under `IncrementBy`, a one-off under the rest.
    @test EcoSISTEM._attachchange(IncrementBy(0.5K / year), layer) isa
          SteadyLayerChange
    @test_throws ErrorException EcoSISTEM._attachchange(OffsetBy(0.5K), layer)
    # ...and something that is neither a constant nor a known shape says so, rather than reporting a
    # missing method on an internal.
    err = try
        EcoSISTEM._attachchange(IncrementBy("tomorrow"), layer)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("not a change shape", err.msg)
end

# **Regression tests for the half-month phase bug - fixed 2026-08-05.** These were written
# deliberately failing, against the unfixed behaviour, and passed **unchanged** when the fix landed.
# That is why their expectations must not be edited: they are the statement of what correct means,
# written before the code that satisfies them.
#
# **The bug they pin.** A monthly series applied every slice half a month early, and slice 1 for half
# as long as the others. Coordinates are built 1-based (`_mkstackaxis` -> `Ti((1:12) .*
# month_mean_duration)`) as slice *identifiers*, `origin` defaults to `first(times)` = 1 month, and
# the lookup took the **nearest** coordinate - which treats each one as a slice *centre*, so
# transitions landed on midpoints instead of month boundaries. The fix (`_current`, `LayerChange.jl`)
# reads a coordinate as *when its slice becomes current* and floors, with a drift tolerance standing
# in for what nearest was implicitly buying.
#
# **Why it survived so long, and why these tests are sub-monthly**: at `timestep =
# 1month_mean_duration` every step lands exactly on a stored coordinate, so the right slice was picked
# every time and the bug was invisible. Every timestep in this repo is `1month_mean_duration`, so
# nothing exercised it. **A fix verified only at monthly steps would have proved nothing** - keep
# these at sub-monthly elapsed times.
@testset "a monthly series applies each slice for its whole month" begin
    stack = _teststack(273.0K, 1.0K)
    layer = _testregime(274.0K, Temperature)
    series = EcoSISTEM._attachchange(ReplaceWith(SeriesChange(stack)), layer)

    # Slice `k` represents month `k`, so it should be current for the whole of that month:
    # elapsed ∈ [(k-1)*month, k*month).
    monthof(k) = min(12,
                     floor(Int, ustrip(uconvert(month_mean_duration, k))) + 1)

    # These three passed even against the bug - the first half of each month landed correctly.
    @test EcoSISTEM._seriesindex(series, 0.0day) == 1
    @test EcoSISTEM._seriesindex(series, 10.0day) == 1
    @test EcoSISTEM._seriesindex(series, 45.0day) == 2

    # These are the ones that failed. 20 and 30 days elapsed are still inside the *first* month (a
    # month being 30.4375 d), so January's slice must still be current - the midpoint rule had already
    # switched to February's at 15.2 d.
    @test EcoSISTEM._seriesindex(series, 20.0day) == 1
    @test EcoSISTEM._seriesindex(series, 30.0day) == 1
    # ...and 60 days is still inside month 2 (which ends at 60.875 d), 90 inside month 3 (91.3 d).
    @test EcoSISTEM._seriesindex(series, 60.0day) == 2
    @test EcoSISTEM._seriesindex(series, 90.0day) == 3

    # And the half-length first slice, stated directly: slice 1 must cover a *whole* month, like
    # every other slice. Under the midpoint rule it was abandoned after 15.2 days.
    @test EcoSISTEM._seriesindex(series, 30.4day) == 1

    # The general property, so a fix cannot satisfy the cases above by accident: across a whole year
    # of daily steps, the current slice must always be the month elapsed time actually falls in.
    # `HoldAtEnd` only so the sweep can reach day 364 - see the coverage bug below for why it
    # otherwise cannot.
    held = EcoSISTEM._attachchange(ReplaceWith(SeriesChange(stack,
                                                            atend = HoldAtEnd())),
                                   layer)
    wrong = [d
             for d in 0:364
             if EcoSISTEM._seriesindex(held, d * 1.0day) != monthof(d * 1.0day)]
    @test isempty(wrong)
end

# **The same root cause, but this one was loud** - also fixed 2026-08-05, and also written before
# the fix. `_seriesreach` added only *half* a gap past the last coordinate, to match the
# nearest-coordinate rule, so a 12-slice monthly series covered 11.5 months: a full year of simulation
# against a twelve-month climatology errored with "this series ran out" a fortnight early, and the
# twelfth slice was never fully reached under the default `ErrorAtEnd`.
#
# Kept separate from the phase test above because this one *throws* rather than quietly returning
# the wrong slice, so it is the symptom a user hits first - and a fix to the lookup alone would have
# left it standing.
@testset "a 12-slice monthly series covers a full twelve months" begin
    stack = _teststack(273.0K, 1.0K)
    layer = _testregime(274.0K, Temperature)
    series = EcoSISTEM._attachchange(ReplaceWith(SeriesChange(stack)), layer)

    # `_seriesreach` is an **absolute** time (last slice + half a gap = 12.5 months), not a
    # duration - subtract `origin` to get the elapsed time actually covered. Comparing the absolute
    # value against `12month_mean_duration` passes vacuously, which is exactly the mistake this comment exists to
    # stop the next person repeating.
    covered = EcoSISTEM._seriesreach(series) - series.origin
    @test covered >= uconvert(s, 12.0month_mean_duration)          # was 11.5 months

    # And the observable consequence: a full year of simulation must not run off the end of a
    # twelve-month climatology. Caught rather than left to throw, so this registers as a clean
    # failure rather than an error.
    reached = try
        EcoSISTEM._seriesindex(series, 11.99month_mean_duration)
        true
    catch
        false
    end
    @test reached
end

@testset "A SeriesLayerChange is indexed by elapsed time" begin
    stack = _teststack(273.0K, 1.0K)

    # The wrap fix. The `StackCycle` this replaces reset to slice 1 and discarded the overshoot,
    # so 13 months into a 12-slice series landed on slice 1; a true modulus lands on slice 2.
    layer = _testregime(274.0K, Temperature)
    repeating = EcoSISTEM._attachchange(ReplaceWith(SeriesChange(stack,
                                                                 atend = RepeatAtEnd())),
                                        layer)
    @test EcoSISTEM._seriesindex(repeating, 13.0month_mean_duration) == 2
    @test EcoSISTEM._seriesindex(repeating, 12.0month_mean_duration) == 1
    @test EcoSISTEM._seriesindex(repeating, 24.0month_mean_duration) == 1
    @test EcoSISTEM._seriesindex(repeating, 0.0month_mean_duration) == 1

    # Indexed by elapsed time, not by call count: twelve one-month steps and one twelve-month
    # step reach the same slice. A cursor advanced once per call cannot do this, which is the whole
    # reason the series replaced one.
    accumulated = sum(fill(uconvert(s, 1.0month_mean_duration), 12))
    @test EcoSISTEM._seriesindex(repeating, accumulated) ==
          EcoSISTEM._seriesindex(repeating, 12.0month_mean_duration)

    # Applying it writes the slice, in the layer's own unit.
    EcoSISTEM._applychange!(repeating, layer, 2.0month_mean_duration,
                            1.0month_mean_duration)
    @test all(≈(276.0K), layer.matrix)
end

@testset "A SeriesLayerChange honours its end-of-series policy" begin
    stack = _teststack(273.0K, 1.0K)
    attach(policy) = EcoSISTEM._attachchange(ReplaceWith(SeriesChange(stack,
                                                                      atend = policy)),
                                             _testregime(274.0K,
                                                         Temperature))

    # Inside the series all three agree.
    for policy in (ErrorAtEnd(), HoldAtEnd(), RepeatAtEnd())
        @test EcoSISTEM._seriesindex(attach(policy), 11.0month_mean_duration) ==
              12
    end
    # Past the end they diverge: the default refuses, `HoldAtEnd` keeps the last slice,
    # `RepeatAtEnd` wraps.
    err = try
        EcoSISTEM._seriesindex(attach(ErrorAtEnd()), 13.0month_mean_duration)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("ran out", err.msg)
    @test EcoSISTEM._seriesindex(attach(HoldAtEnd()),
                                 13.0month_mean_duration) == 12
    @test EcoSISTEM._seriesindex(attach(HoldAtEnd()),
                                 120.0month_mean_duration) == 12
    @test EcoSISTEM._seriesindex(attach(RepeatAtEnd()),
                                 13.0month_mean_duration) == 2

    # An end-of-series policy *is* a type, so anything that is not one is refused by the signature,
    # where it is written - rather than being looked up from a list of accepted names and reported
    # only once the series runs out.
    @test_throws TypeError SeriesChange(stack, atend = :wrap)
end

@testset "A SeriesLayerChange validates what it cannot index" begin
    stack = _teststack(273.0K, 1.0K)
    layer = _testregime(274.0K, Temperature)

    # A calendar axis is now read rather than refused: its dates become elapsed coordinates
    # measured from the first slice, and that first date is kept as the series' calendar identity,
    # which is what a run's epoch is resolved against.
    calendar = DimArray(stack,
                        (Y(NoLookup()), X(NoLookup()),
                         Ti(Dates.DateTime(2000, 1, 1):Dates.Month(1):Dates.DateTime(2000,
                                                                                     12,
                                                                                     1))))
    dated = EcoSISTEM._attachchange(ReplaceWith(SeriesChange(calendar)), layer)
    @test dated isa SeriesLayerChange
    @test dated.calendar == DatedSeries(Dates.DateTime(2000, 1, 1))
    # Elapsed coordinates run from zero (the first slice) to the real gap to the last - 2000 was a
    # leap year, so January to December is 335 days, not eleven mean months.
    @test dated.times[1] == 0.0s
    @test dated.times[end] ≈ uconvert(s, 335.0 * u"d")
    @test dated.origin == 0.0s
    # ...and giving each slice an elapsed time explicitly still overrides the axis, leaving a series
    # with no calendar identity at all.
    plain = EcoSISTEM._attachchange(ReplaceWith(SeriesChange(calendar,
                                                             times = (1:12) .*
                                                                     month_mean_duration)),
                                    layer)
    @test plain isa SeriesLayerChange
    @test plain.calendar == UndatedSeries()

    # An unevenly spaced series has no period to cycle over (an ERA read is exactly this), so
    # `:repeat` is refused rather than guessed at - `:hold` remains available.
    uneven = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 20] .* month_mean_duration
    @test_throws ErrorException EcoSISTEM._attachchange(ReplaceWith(SeriesChange(stack,
                                                                                 times = uneven,
                                                                                 atend = RepeatAtEnd())),
                                                        layer)
    @test EcoSISTEM._attachchange(ReplaceWith(SeriesChange(stack,
                                                           times = uneven,
                                                           atend = HoldAtEnd())),
                                  layer) isa SeriesLayerChange

    # A series on another grid, or in another dimension, is caught at attach rather than in a
    # broadcast on the first timestep.
    @test_throws ErrorException EcoSISTEM._attachchange(ReplaceWith(SeriesChange(_teststack(273.0K,
                                                                                            1.0K,
                                                                                            (2,
                                                                                             2)))),
                                                        layer)
    @test_throws ErrorException EcoSISTEM._attachchange(ReplaceWith(SeriesChange(_teststack(0.0mm,
                                                                                            1.0mm))),
                                                        layer)
    # One time per slice, strictly increasing.
    @test_throws ErrorException EcoSISTEM._attachchange(ReplaceWith(SeriesChange(stack,
                                                                                 times = (1:5) .*
                                                                                         month_mean_duration)),
                                                        layer)
end

@testset "A SeriesLayerChange reads its own time axis, and its mode" begin
    # A real `Ti` lookup in months is used as given...
    stack = _teststack(273.0K, 1.0K)
    layer = _testregime(274.0K, Temperature)
    real = DimArray(stack,
                    (Y(NoLookup()), X(NoLookup()),
                     Ti((1:12) .* month_mean_duration)))
    onaxis = EcoSISTEM._attachchange(ReplaceWith(SeriesChange(real,
                                                              atend = RepeatAtEnd())),
                                     layer)
    @test EcoSISTEM._seriesindex(onaxis, 13.0month_mean_duration) == 2
    # ...and elapsed time is measured from the series' own start, so a 1-based axis and a 0-based one
    # both begin at their first slice rather than one being an off-by-one from the other.
    zerobased = DimArray(stack,
                         (Y(NoLookup()), X(NoLookup()),
                          Ti((0:11) .* month_mean_duration)))
    @test EcoSISTEM._seriesindex(EcoSISTEM._attachchange(ReplaceWith(SeriesChange(zerobased,
                                                                                  atend = RepeatAtEnd())),
                                                         layer),
                                 13.0month_mean_duration) == 2

    # Under `OffsetBy` the slices are an interval added to the layer's captured values, so the same
    # stack means something different - and, on an affine axis, must be read in `K` not `°C`.
    offset = _testregime(274.0K, Temperature)
    deviation = EcoSISTEM._attachchange(OffsetBy(SeriesChange(_teststack(0.0K,
                                                                         1.0K),
                                                              atend = RepeatAtEnd())),
                                        offset)
    EcoSISTEM._applychange!(deviation, offset, 2.0month_mean_duration,
                            1.0month_mean_duration)
    @test all(≈(277.0K), offset.matrix)
end

@testset "Changes compose additively" begin
    # The motivating case, and the forward requirement Phase 2 was rejected on: a monthly pattern
    # varying over the year *plus* a steady multi-year trend offsetting the whole pattern. They are
    # additive, not contradictory.
    layer = _testregime(274.0K, Temperature)
    stack = _teststack(273.0K, 1.0K)          # slice i is 273 K + i K
    combined = ReplaceWith(SeriesChange(stack, atend = RepeatAtEnd())) +
               IncrementBy(1.2K / year)
    EcoSISTEM.setchange!(layer, combined)
    @test layer.change isa SumOfLayerChanges{AbsoluteChange}

    # At 2 months: slice 3 (276 K) plus two months of trend.
    EcoSISTEM._layerupdate!(layer, 2.0month_mean_duration,
                            1.0month_mean_duration)
    @test layer.matrix[1, 1] ≈ 276.0K + 1.2K / year * 2.0month_mean_duration

    # The failure this design exists to prevent: applied one after another, the absolute write
    # each step would erase whatever the trend had accumulated, so the trend would never grow. Here
    # it does, because the sum is over values as functions of elapsed time and is written once.
    EcoSISTEM._layerupdate!(layer, 14.0month_mean_duration,
                            1.0month_mean_duration)
    @test layer.matrix[1, 1] ≈ 276.0K + 1.2K / year * 14.0month_mean_duration
    # ...and the series has wrapped back to the same slice, so the whole difference is the trend.
    @test EcoSISTEM._seriesindex(first(layer.change.parts),
                                 14.0month_mean_duration) ==
          EcoSISTEM._seriesindex(first(layer.change.parts),
                                 2.0month_mean_duration)

    # With no absolute part the sum is relative - read from the layer's captured values.
    relative = _testregime(274.0K, Temperature)
    EcoSISTEM.setchange!(relative,
                         OffsetBy(PatternedChange(5.0K, 1.0year)) +
                         IncrementBy(1.2K / year))
    @test relative.change isa SumOfLayerChanges{RelativeChange}
    EcoSISTEM._layerupdate!(relative, 0.25year, 0.25year)
    @test all(≈(274.0K + 5.0K + 1.2K / year * 0.25year), relative.matrix)

    # Combining flattens, so one sum of three means what it says.
    @test length((ReplaceWith(SeriesChange(stack)) +
                  OffsetBy(PatternedChange(1.0K,
                                           1.0year)) +
                  IncrementBy(1.0K / year)).specs) == 3

    # `+` is the spelling; the type it builds is supported but *unexported*, for the one case a `+`
    # chain cannot express - a collection whose length is not known when the code is written.
    parts = [ReplaceWith(SeriesChange(stack)), IncrementBy(1.0K / year)]
    @test EcoSISTEM.CombinedChange(parts...).specs ==
          (parts[1] + parts[2]).specs
    @test !Base.isexported(EcoSISTEM, :CombinedChange)
    @test Base.ispublic(EcoSISTEM, :CombinedChange)
    # ...and a sum of one is not a sum at all.
    @test_throws ErrorException EcoSISTEM.CombinedChange(parts[1])
end

@testset "Composition refuses what it cannot sum" begin
    layer = _testregime(274.0K, Temperature)
    stack = _teststack(273.0K, 1.0K)

    # Two positions cannot be added - only intervals can be added to a position.
    err = try
        EcoSISTEM.setchange!(layer,
                             ReplaceWith(SeriesChange(stack)) +
                             ReplaceWith(SeriesChange(stack)))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("two positions cannot be added", err.msg)

    # A rate contributes its *integral*, so it needs one in closed form. A steady rate has one;
    # an oscillating rate would need the integral of its own shape, and approximating it would make
    # the result depend on the timestep - so it is refused rather than approximated.
    err = try
        EcoSISTEM.setchange!(layer,
                             ReplaceWith(SeriesChange(stack)) +
                             IncrementBy(PatternedChange(1.0K / year,
                                                         1.0year)))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("integral", err.msg)
    # ...and it is refused at attach, not on some later timestep.
    @test layer.change isa NoLayerChange
end

@testset "A series' stored gaps are cleaned with the layer's" begin
    # Regression. A supply's data gaps become "no resource" when it is placed in a habitat, but a
    # series holds the slices it will write on *later* steps. Cleaning only the matrix would let the
    # first update put the `NaN`s straight back - one step late, and where they reach the resource
    # arithmetic and `populate!`'s supply-weighted placement.
    stack = cat(fill(10.0kJ / day, 3, 3), fill(20.0kJ / day, 3, 3), dims = 3)
    stack[1, 1, :] .= NaN * kJ / day
    supply = EcoSISTEM._setseries!(Supply{SolarRadiation}(stack[:, :, 1]),
                                   stack)
    @test any(isnan, supply.matrix)
    @test EcoSISTEM._hasgaps(supply.change)

    cleaned = EcoSISTEM._zerogaps(supply)
    @test !any(isnan, cleaned.matrix)
    @test !EcoSISTEM._hasgaps(cleaned.change)
    @test iszero(cleaned.matrix[1, 1])

    # ...and it stays clean once the series actually writes a later slice.
    EcoSISTEM._layerupdate!(cleaned, 1.0month_mean_duration,
                            1.0month_mean_duration)
    @test !any(isnan, cleaned.matrix)
    @test iszero(cleaned.matrix[1, 1])
    @test cleaned.matrix[2, 2] == 20.0kJ / day

    # The original is untouched - cleaning rebuilds rather than reaching back into the caller's
    # supply, and that has to hold for the stored slices too.
    @test any(isnan, supply.matrix)
    @test EcoSISTEM._hasgaps(supply.change)
end

@testset "Varying declares a change alongside a spec" begin
    spec = UniformSpec(274.0K, axis = Temperature)
    change = IncrementBy(2.0K / year)

    v = Varying(spec, change)
    @test v.spec === spec
    @test v.change === change

    # A layer carries exactly one change, so the wrapper cannot nest...
    @test_throws ErrorException Varying(v, change)
    # ...and it declares a change for one layer, so it cannot wrap a tuple of them. A change is
    # checked against a single layer's unit, so one declaration spanning several would mean
    # different things on each.
    err = try
        Varying((spec, spec), change)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("one layer", err.msg)

    # The second argument must be a change, not an arbitrary value.
    @test_throws MethodError Varying(spec, 3)

    # Splitting names its halves, and for a tuple both are tuples aligned positionally, so an
    # element and its change stay together.
    split = EcoSISTEM._splitvarying((v, spec))
    @test split.spec == (spec, spec)
    @test split.change == (change, nothing)
    # ...and anything unwrapped passes straight through.
    @test EcoSISTEM._splitvarying(spec) == (spec = spec, change = nothing)
end

@testset "A declared change is added to a series, not put in its place" begin
    # The promise Phase 2 made when it rejected this outright. A layer read from a monthly stack
    # already carries a `SeriesLayerChange`; a change declared on top is *added* to it, because a seasonal
    # pattern and a multi-year trend offset one another rather than competing. Replacing the series
    # would freeze the layer on one slice while still looking like a working seasonal layer.
    stack = cat((fill(i * 1.0K, 5, 5) for i in 1:12)..., dims = 3)
    layer = EcoSISTEM._setseries!(ContinuousRegime(fill(1.0K, 5, 5), 1.0km,
                                                   NoLayerChange()),
                                  stack)
    EcoSISTEM._applydeclared!(layer, IncrementBy(2.0K / year))
    @test layer.change isa SumOfLayerChanges{AbsoluteChange}
    @test any(p -> p isa SeriesLayerChange, layer.change.parts)
    @test any(p -> p isa SteadyLayerChange, layer.change.parts)

    # ...and both halves are live: the slice still advances with elapsed time, and the trend
    # accumulates on top of it rather than being erased by each absolute write.
    EcoSISTEM._layerupdate!(layer, 2.0month_mean_duration,
                            1.0month_mean_duration)
    @test layer.matrix[1, 1] ≈ 3.0K + 2.0K / year * 2.0month_mean_duration
    EcoSISTEM._layerupdate!(layer, 5.0month_mean_duration,
                            3.0month_mean_duration)
    @test layer.matrix[1, 1] ≈ 6.0K + 2.0K / year * 5.0month_mean_duration

    # A layer with no change of its own takes one directly.
    flat = ContinuousRegime(fill(274.0K, 5, 5), 1.0km, NoLayerChange())
    EcoSISTEM._applydeclared!(flat, IncrementBy(2.0K / year))
    @test flat.change isa SteadyLayerChange
    # ...and declaring nothing leaves it alone.
    EcoSISTEM._applydeclared!(flat, nothing)
    @test flat.change isa SteadyLayerChange
end

@testset "Condition loss" begin
    eco = Test1Ecosystem()
    # A regime carrying a HabitatLoss change whose rate destroys every active cell
    # over one timestep (rate * 1month_mean_duration == 1 -> loss probability 1).
    change = EcoSISTEM.LegacyLoss(1.0 / month_mean_duration)
    losshab = EcoSISTEM.ContinuousRegime(fill(1.0K, 10, 10), 1.0km, change)
    @test (@test_deprecated EcoSISTEM.HabitatLoss(eco, losshab,
                                                  1month_mean_duration)) === eco
    @test all(iszero, eco.habitat.supply.matrix)
    @test all(iszero, eco.abundances.matrix)

    # It is not a layer change: driving it from the update loop would need the ecosystem it
    # mutates, and would draw at random on every rank independently.
    @test_throws ErrorException EcoSISTEM._layerupdate!(losshab,
                                                        1.0month_mean_duration,
                                                        1month_mean_duration)
end

# ---------------------------------------------------------------------------
# The simulation epoch
# ---------------------------------------------------------------------------
# A 12-slice climatology whose slice `k` holds `279 + k` K, so the *value* names the slice it came
# from and an assertion can say which slice was selected rather than only that something changed.
function _namedclimatology(dim = (5, 5))
    stack = cat((fill((279.0 + i) * K, dim) for i in 1:12)..., dims = 3)
    return DimArray(stack,
                    (Y(NoLookup()), X(NoLookup()),
                     Ti((1:12) .* month_mean_duration)))
end

# Which slice a layer is showing, read back out of the value `_namedclimatology` encoded.
_shownslice(layer) = round(Int, ustrip(K, layer.matrix[1, 1]) - 279)

# A layer carrying that climatology, declared as months of the year.
function _climatologylayer(; atend = RepeatAtEnd(),
                           calendar = MonthOfYearSeries())
    layer = _testregime(280.0K, Temperature)
    EcoSISTEM.setchange!(layer,
                         ReplaceWith(SeriesChange(_namedclimatology(),
                                                  calendar = calendar,
                                                  atend = atend)))
    return layer
end

@testset "A series says what its coordinates mean" begin
    layer = _testregime(274.0K, Temperature)

    # What this distinction prevents: a lookup-less stack read as *months of the year*, which makes
    # a synthetic 10-slice stack silently "months 1 to 10". Monthly is its spacing - that is what the
    # v0.4.0 stack walk assumes - but it has no calendar identity, so nothing may phase-lock it.
    bare = EcoSISTEM._attachchange(ReplaceWith(SeriesChange(_teststack(273.0K,
                                                                       1.0K,
                                                                       (5, 5),
                                                                       10))),
                                   layer)
    @test bare.calendar == UndatedSeries()
    @test bare.times == collect((1:10) .* uconvert(s, 1.0month_mean_duration))

    # Elapsed-time coordinates are genuinely ambiguous, so they infer the reading that phases
    # nothing; a climatology opts in explicitly.
    @test EcoSISTEM._attachchange(ReplaceWith(SeriesChange(_namedclimatology())),
                                  layer).calendar == UndatedSeries()
    @test EcoSISTEM._attachchange(ReplaceWith(SeriesChange(_namedclimatology(),
                                                           calendar = MonthOfYearSeries())),
                                  layer).calendar == MonthOfYearSeries()

    # ...and the claim is checked where it is made: coordinates that are not whole month numbers
    # from 1 to 12 are not months of the year, whatever the caller says.
    @test_throws ErrorException EcoSISTEM._attachchange(ReplaceWith(SeriesChange(_namedclimatology(),
                                                                                 times = (1:12) .*
                                                                                         year,
                                                                                 calendar = MonthOfYearSeries())),
                                                        layer)
end

@testset "`origin` is accepted only where an epoch is powerless" begin
    layer = _testregime(274.0K, Temperature)
    # An undated series has no other way to say where elapsed zero sits, so `origin` is its knob.
    shifted = EcoSISTEM._attachchange(ReplaceWith(SeriesChange(_namedclimatology(),
                                                               origin = 3.0month_mean_duration)),
                                      layer)
    @test shifted.origin == uconvert(s, 3.0month_mean_duration)

    # For the other two the epoch fixes the phase, and a second knob for the same thing could only
    # contradict it - so it is refused at the line that wrote it.
    for calendar in (MonthOfYearSeries(), DatedSeries(Dates.DateTime(2000, 1, 1)))
        err = try
            EcoSISTEM._attachchange(ReplaceWith(SeriesChange(_namedclimatology(),
                                                             origin = 3.0month_mean_duration,
                                                             calendar = calendar)),
                                    layer)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("epoch", err.msg)
    end
end

@testset "An epoch phases a climatology by calendar month" begin
    # The point of the whole phase: without an epoch every run starts in January.
    layer = _climatologylayer()
    EcoSISTEM._repointseries!(layer, nothing)
    EcoSISTEM._layerupdate!(layer, 0.0s, 1.0day)
    @test _shownslice(layer) == 1

    # With one, a run starts on the slice for its own month - all twelve of them.
    for month in 1:12
        layer = _climatologylayer()
        EcoSISTEM._repointseries!(layer, Dates.DateTime(2015, month, 1))
        EcoSISTEM._layerupdate!(layer, 0.0s, 1.0day)
        @test _shownslice(layer) == month
    end

    # Matching is by month *number*, not by elapsed duration into the year, and this is the case
    # that separates them: the slices are spaced by `month_mean_duration` (30.44 d) so six of them
    # reach 182.6 days, while the real 1 July is day 181 - matching proportionally would select the
    # June slice. Pinned because it is the exact off-by-one the layer-units work found elsewhere.
    @test 6 * 30.4375 > 181
    layer = _climatologylayer()
    EcoSISTEM._repointseries!(layer, Dates.DateTime(2015, 7, 1))
    EcoSISTEM._layerupdate!(layer, 0.0s, 1.0day)
    @test _shownslice(layer) == 7

    # A July run wraps round to January after six months, rather than running off the end.
    EcoSISTEM._layerupdate!(layer, 6.0month_mean_duration, 1.0day)
    @test _shownslice(layer) == 1

    # An undated series is left exactly where it was, epoch or no epoch - nothing to bind to.
    undated = _climatologylayer(calendar = UndatedSeries())
    origin = undated.change.origin
    EcoSISTEM._repointseries!(undated, Dates.DateTime(2015, 7, 1))
    @test undated.change.origin == origin
    EcoSISTEM._layerupdate!(undated, 0.0s, 1.0day)
    @test _shownslice(undated) == 1
end

@testset "An epoch places a dated series at its own date" begin
    dates = Dates.DateTime(2010, 1, 1):Dates.Month(1):Dates.DateTime(2010, 12,
                                                                     1)
    stack = DimArray(cat((fill((279.0 + i) * K, 5, 5) for i in 1:12)...,
                         dims = 3),
                     (Y(NoLookup()), X(NoLookup()), Ti(dates)))
    build() = begin
        layer = _testregime(280.0K, Temperature)
        EcoSISTEM.setchange!(layer,
                             ReplaceWith(SeriesChange(stack,
                                                      atend = HoldAtEnd())))
        layer
    end

    # Elapsed zero is the slice covering the epoch's own date, so an April start shows April.
    layer = build()
    EcoSISTEM._repointseries!(layer, Dates.DateTime(2010, 4, 1))
    EcoSISTEM._layerupdate!(layer, 0.0s, 1.0day)
    @test _shownslice(layer) == 4

    # An epoch *before* the series begins is **not** an error - see the dedicated testset below.
    # It gives the run a lead-in during which the layer's own values stand, which is a thing to want
    # rather than a mistake. Here it is enough that it is accepted and dates the origin negatively.
    early = build()
    EcoSISTEM._repointseries!(early, Dates.DateTime(2009, 6, 1))
    @test early.change.origin < 0.0s

    # ...and a *cyclic* series has no beginning to precede at all, so any date phases it.
    @test EcoSISTEM._repointseries!(_climatologylayer(),
                                    Dates.DateTime(1850, 6, 1)) isa
          EcoSISTEM.AbstractLayer
end

@testset "An epoch is adopted if unambiguous and asked for if not" begin
    onedate(y) = begin
        dates = Dates.DateTime(y, 1, 1):Dates.Month(1):Dates.DateTime(y, 12, 1)
        stack = DimArray(_teststack(273.0K, 1.0K),
                         (Y(NoLookup()), X(NoLookup()), Ti(dates)))
        layer = _testregime(280.0K, Temperature)
        EcoSISTEM.setchange!(layer,
                             ReplaceWith(SeriesChange(stack,
                                                      atend = HoldAtEnd())))
        layer
    end
    habitat(regime, supply) = (regime = regime, supply = supply)

    # No dated series at all: no epoch, and everything behaves as a run that never mentions dates.
    plain = habitat(_climatologylayer(), _climatologylayer())
    @test isnothing(EcoSISTEM._resolveepoch(plain, nothing))

    # Exactly one real start date: adopt it.
    single = habitat(onedate(2001), _climatologylayer())
    @test EcoSISTEM._resolveepoch(single, nothing) == Dates.DateTime(2001, 1, 1)

    # Two that disagree: refuse to guess, and name the candidates.
    both = habitat(onedate(2001), onedate(2005))
    err = try
        EcoSISTEM._resolveepoch(both, nothing)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("2001", err.msg) && occursin("2005", err.msg)

    # An explicit epoch always wins, including over a disagreement.
    @test EcoSISTEM._resolveepoch(both, Dates.DateTime(2003, 1, 1)) ==
          Dates.DateTime(2003, 1, 1)
    @test EcoSISTEM._resolveepoch(single, Dates.DateTime(2003, 1, 1)) ==
          Dates.DateTime(2003, 1, 1)
end

@testset "A run is checked against the series driving it up front" begin
    # `ErrorAtEnd` would fail anyway, but at the step it happens - so it is reported before the
    # first one instead.
    err = try
        EcoSISTEM._checkcoverage(_climatologylayer(atend = ErrorAtEnd()),
                                 5.0year)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("does not cover the run", err.msg)

    # `HoldAtEnd` never fails, which is why it earns a warning instead - with the proportion in
    # it, since that is the number that says whether it was intended.
    held = _climatologylayer(atend = HoldAtEnd())
    @test_logs (:warn, r"80.0% of the simulation") EcoSISTEM._checkcoverage(held,
                                                                            5.0year)

    # A cycling series always covers the run, and a run inside the series has nothing to report.
    @test_logs EcoSISTEM._checkcoverage(_climatologylayer(atend = RepeatAtEnd()),
                                        5.0year)
    @test_logs EcoSISTEM._checkcoverage(_climatologylayer(atend = HoldAtEnd()),
                                        6.0month_mean_duration)
end

@testset "A run starts in the state its series describes" begin
    # The defect this closes: a change only reaches a layer's `matrix` when `_applychange!` runs,
    # and `update!` does that at the *end* of a timestep - so an unprimed series layer spent the
    # whole of step one holding whatever its builder left behind, and step one's population dynamics
    # ran against an environment nothing had asked for.
    layer = _climatologylayer()
    EcoSISTEM._repointseries!(layer, Dates.DateTime(2015, 7, 1))
    @test _shownslice(EcoSISTEM._primeseries!(layer, 0.0s)) == 7

    # Without an epoch it is the first slice, which is what elapsed zero selects.
    plain = _climatologylayer()
    EcoSISTEM._repointseries!(plain, nothing)
    @test _shownslice(EcoSISTEM._primeseries!(plain, 0.0s)) == 1

    # Priming is idempotent - it evaluates a pure function of elapsed time rather than compounding.
    @test _shownslice(EcoSISTEM._primeseries!(plain, 0.0s)) == 1

    # ...but a **rate** must not be primed: it accumulates `value × timestep`, and at the start of a
    # run nothing has accumulated. Priming would add a phantom step of drift before the first step.
    rate = _testregime(274.0K, Temperature)
    EcoSISTEM.setchange!(rate, IncrementBy(2.0K / year))
    EcoSISTEM._primeseries!(rate, 0.0s)
    @test all(≈(274.0K), rate.matrix)
    EcoSISTEM._layerupdate!(rate, 1.0year, 1.0year)
    @test all(≈(276.0K), rate.matrix)

    # A relative change *is* primed, against the baseline it captured at attach.
    offset = _testregime(274.0K, Temperature)
    EcoSISTEM.setchange!(offset, OffsetBy(PatternedChange(4.0K, 1.0year)))
    EcoSISTEM._primeseries!(offset, 0.25year)
    @test all(≈(278.0K), offset.matrix)
end

@testset "A supply cannot go negative" begin
    # Not a policy bolted onto `Resource`: it restates what makes something a resource. A resource
    # is rival and consumed against a demand, so a negative amount of it has no meaning - which is
    # why the rule needs no per-axis opt-in and admits no exceptions.
    negative = Supply{SolarRadiation}(fill(-1.0kJ / day, 3, 3))
    err = try
        EcoSISTEM._checksupplybounds(negative)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("negative", err.msg)
    # ...and the message points at the resolution rather than just refusing: a quantity that
    # genuinely takes both signs is not a supply.
    @test occursin("regime side", err.msg)
    @test EcoSISTEM._checksupplybounds(Supply{SolarRadiation}(fill(1.0kJ / day,
                                                                   3, 3))) isa
          AbstractLayer

    # A collection is checked member by member, so one bad supply among several is still caught.
    @test_throws ErrorException EcoSISTEM._checksupplybounds(LayerCollection((Supply{SolarRadiation}(fill(1.0kJ /
                                                                                                          day,
                                                                                                          3,
                                                                                                          3)),
                                                                              Supply{SolarRadiation}(fill(-1.0kJ /
                                                                                                          day,
                                                                                                          3,
                                                                                                          3)))))

    # An absolute replacement is refused at *attach*, because its stored slices are exactly the
    # values the layer will take.
    stack(v) = DimArray(cat((fill(v, 3, 3) for _ in 1:4)..., dims = 3),
                        (Y(NoLookup()), X(NoLookup()),
                         Ti((1:4) .* month_mean_duration)))
    supply() = Supply{SolarRadiation}(fill(10.0kJ / day, 3, 3))
    @test_throws ErrorException EcoSISTEM.setchange!(supply(),
                                                     ReplaceWith(SeriesChange(stack(-5.0kJ /
                                                                                    day))))
    @test_nowarn EcoSISTEM.setchange!(supply(),
                                      ReplaceWith(SeriesChange(stack(5.0kJ /
                                                                     day))))

    # ...but a **regime** may go as negative as it likes: the bound is role-level, not universal.
    regime = _testregime(10.0K, Temperature)
    negativeslices = DimArray(cat((fill(-5.0K, 5, 5) for _ in 1:4)...,
                                  dims = 3),
                              (Y(NoLookup()), X(NoLookup()),
                               Ti((1:4) .* month_mean_duration)))
    @test_nowarn EcoSISTEM.setchange!(regime,
                                      ReplaceWith(SeriesChange(negativeslices)))

    # At run time the response is deliberately different: a supply driven below zero by an
    # increment is *emergent*, and aborting a long simulation when "you cannot have less than none of
    # a consumable" is the right reading would be hostile. So it warns and clamps.
    running = Supply{SolarRadiation}(fill(10.0kJ / day, 3, 3))
    EcoSISTEM.setchange!(running,
                         IncrementBy(-6.0kJ / day / month_mean_duration))
    EcoSISTEM._layerupdate!(running, 1.0month_mean_duration,
                            1.0month_mean_duration)
    @test all(≈(4.0kJ / day), running.matrix)          # still positive, untouched

    # Warned **once**, not once per step - a supply held at zero would otherwise warn on every
    # subsequent step and drown a long run's output. A single `(:warn, ...)` pattern asserts exactly
    # one matching record, so a second warning fails this.
    # Both steps must sit inside **one** `@test_logs`: `maxlog` counts per *logger*, and
    # `@test_logs` installs a fresh one each time - so splitting them resets the counter and would
    # test nothing. (Found the hard way; the production behaviour was right all along.)
    @test_logs (:warn, r"below zero") begin
        EcoSISTEM._layerupdate!(running, 2.0month_mean_duration,
                                1.0month_mean_duration)
        EcoSISTEM._layerupdate!(running, 3.0month_mean_duration,
                                1.0month_mean_duration)
    end
    @test all(iszero, running.matrix)
end

@testset "A condition cannot leave its physical range" begin
    # The asymmetry with supplies is the point, not an inconsistency: a supply hitting zero is
    # *expected* (resources get consumed) and is clamped; a temperature below absolute zero is
    # **impossible**, so it errors. Clamping it would let a physically nonsensical run continue and
    # produce output that looks like a result.
    area = StudyArea(extent = (3.0km, 5.0km), cellsize = 1.0km,
                     verbosity = :silent)
    eco(rate) = begin
        env = GridHabitat(regime = Varying(UniformSpec(285.0K,
                                                       axis = Temperature),
                                           IncrementBy(rate)),
                          supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                               axis = SolarRadiation),
                          area = area)
        spp = build_species(3, tolerance = (285.0K, 5.0K),
                            toleranceaxis = Temperature,
                            demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                            abundance = 300, seed = 1)
        build_ecosystem(spp, env, seed = 1)
    end

    # Caught BEFORE the first timestep: a constant rate over a known duration has a predictable
    # end value, so "shorten your run" beats a simulation that dies part-way through.
    err = try
        simulate!(eco(-1.0K / year), 300.0year, 1.0month_mean_duration)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("before the run ends", err.msg)

    # A run that stays in range is untouched - 285 K falling 1 K/yr for 100 years reaches 185 K.
    @test_nowarn simulate!(eco(-1.0K / year), 100.0year, 1.0month_mean_duration)
    # ...and `Temperature` has no ceiling, so warming as far as you like is fine.
    @test_nowarn simulate!(eco(1.0K / year), 300.0year, 1.0month_mean_duration)

    # The write-site backstop, for a change whose reach cannot be predicted: a `PatternedChange`
    # takes an arbitrary function of elapsed time, so its amplitude does not bound it in general and
    # the pre-flight check deliberately stays silent rather than guessing.
    # This also exercises the **affine** path: the layer holds °C, so absolute zero must be
    # compared as -273.15 °C - converting the bound, not the values.
    celsius = _testregime(5.0u"°C", Temperature)
    EcoSISTEM.setchange!(celsius, OffsetBy(PatternedChange(500.0K, 1.0year)))
    err2 = try
        EcoSISTEM._layerupdate!(celsius, 0.75year, 0.1year)   # the trough of the sinusoid
        nothing
    catch e
        e
    end
    @test err2 isa ErrorException
    @test occursin("physical range", err2.msg)
    @test occursin("-273.15 °C", err2.msg)

    # A **balance** axis is unbounded and must not be refused - this is the case the catalogue's
    # `Category = balance` marker exists to protect.
    balance = _testregime(5.0mm / day, ClimateMoisture)
    EcoSISTEM.setchange!(balance, IncrementBy(-1.0mm / day / year))
    # One call applies **one timestep**, not `elapsed`-worth of drift - the elapsed argument only
    # tells the change what its value is *now*, which for a steady rate is constant.
    @test_nowarn EcoSISTEM._layerupdate!(balance, 100.0year, 100.0year)
    @test all(<(0.0mm / day), balance.matrix)      # genuinely negative, and allowed
end

@testset "Outside its span a series contributes nothing" begin
    # Slices at months 12-23, so an `origin` of zero gives a year of lead-in before the series has
    # anything to say. 500 K is unmistakably the layer's own.
    stack = DimArray(cat((fill((279.0 + i) * K, 5, 5) for i in 1:12)...,
                         dims = 3),
                     (Y(NoLookup()), X(NoLookup()),
                      Ti((12:23) .* month_mean_duration)))
    build(recipe) = begin
        layer = _testregime(500.0K, Temperature)
        EcoSISTEM.setchange!(layer, recipe)
        EcoSISTEM._primeseries!(layer, 0.0s)
        layer
    end
    shown(l) = l.matrix[1, 1]
    at(l, m) = (
                EcoSISTEM._layerupdate!(l, m * month_mean_duration,
                                        1.0month_mean_duration);
                shown(l))

    # Before its first slice the layer stands - including at build, which priming must not
    # overwrite with a slice the series has not reached.
    lead = build(ReplaceWith(SeriesChange(stack,
                                          origin = 0.0month_mean_duration,
                                          atend = HoldAtEnd())))
    @test shown(lead) == 500.0K
    @test at(lead, 0) == 500.0K
    @test at(lead, 11) == 500.0K
    @test at(lead, 12) == 280.0K          # the first slice, exactly when it becomes current
    @test at(lead, 13) == 281.0K

    # `RevertToLayer` is the same rule at the far end: the last slice serves its full month, then
    # the layer is given back.
    revert = build(ReplaceWith(SeriesChange(stack, atend = RevertToLayer())))
    @test shown(revert) == 280.0K         # default origin, so no lead-in
    @test at(revert, 11) == 291.0K
    @test at(revert, 12) == 291.0K        # still in force for its own month
    @test at(revert, 13) == 500.0K        # ...and then the layer is back

    # The default is untouched: with `origin` defaulting to the first slice, elapsed time can
    # never fall short of it, so no run that does not ask for a lead-in gets one.
    plain = build(ReplaceWith(SeriesChange(stack, atend = HoldAtEnd())))
    @test shown(plain) == 280.0K
    @test at(plain, 12) == 291.0K

    # ...but omitting `origin` and passing zero are **not** the same thing, and this is the pairing
    # that says so: omitting means "the first slice", zero means "coordinate zero", which is before
    # it. On an ordinary 1-12 climatology that difference is a whole month of lead-in.
    ordinary = DimArray(cat((fill((279.0 + i) * K, 5, 5) for i in 1:12)...,
                            dims = 3),
                        (Y(NoLookup()), X(NoLookup()),
                         Ti((1:12) .* month_mean_duration)))
    @test shown(build(ReplaceWith(SeriesChange(ordinary, atend = HoldAtEnd())))) ==
          280.0K
    @test shown(build(ReplaceWith(SeriesChange(ordinary,
                                               origin = 0.0month_mean_duration,
                                               atend = HoldAtEnd())))) ==
          500.0K

    # "Contributes nothing" is expressed in each mode's own terms, which is what lets the three
    # `_applychange!` methods stay as they are: a captured value, a zero offset, a zero rate.
    offsets = DimArray(cat((fill(i * 1.0K, 5, 5) for i in 1:12)..., dims = 3),
                       (Y(NoLookup()), X(NoLookup()),
                        Ti((12:23) .* month_mean_duration)))
    off = build(OffsetBy(SeriesChange(offsets,
                                      origin = 0.0month_mean_duration,
                                      atend = RevertToLayer())))
    @test at(off, 0) == 500.0K            # the baseline, offset by nothing
    @test at(off, 12) == 501.0K

    rates = DimArray(cat((fill(12.0K / year, 5, 5) for i in 1:12)..., dims = 3),
                     (Y(NoLookup()), X(NoLookup()),
                      Ti((12:23) .* month_mean_duration)))
    inc = build(IncrementBy(SeriesChange(rates,
                                         origin = 0.0month_mean_duration,
                                         atend = RevertToLayer())))
    @test at(inc, 0) == 500.0K
    @test at(inc, 6) == 500.0K            # a rate out of span accumulates *nothing*
    @test at(inc, 12) ≈ 501.0K
end

@testset "An epoch before a dated series is not an error" begin
    # Nothing is invalid: the layer has values of its own and only the *series* has nothing to say
    # yet, which is exactly the case the layer stands for. "The record starts in 1970, run from 1950
    # on the spec's own climatology" is a thing to want.
    dates = Dates.DateTime(1970, 1, 1):Dates.Month(1):Dates.DateTime(1970, 12,
                                                                     1)
    stack = DimArray(cat((fill((279.0 + i) * K, 5, 5) for i in 1:12)...,
                         dims = 3),
                     (Y(NoLookup()), X(NoLookup()), Ti(dates)))
    layer = _testregime(500.0K, Temperature)
    EcoSISTEM.setchange!(layer,
                         ReplaceWith(SeriesChange(stack, atend = HoldAtEnd())))
    @test_nowarn EcoSISTEM._repointseries!(layer, Dates.DateTime(1950, 1, 1))
    EcoSISTEM._primeseries!(layer, 0.0s)
    @test layer.matrix[1, 1] == 500.0K
    EcoSISTEM._layerupdate!(layer, 20.0year, 1.0day)
    @test layer.matrix[1, 1] == 280.0K
end

@testset "Coverage is checked against where the run really ends" begin
    # `simulate!` takes `length(0s:timestep:duration)` steps - a range including both ends - so a
    # twelve-month run in one-month steps advances the clock *thirteen* times. Checking against
    # `duration` would have passed runs that then failed mid-flight, which is precisely what this
    # check exists to pre-empt, so it is checked against the elapsed time the run actually reaches.
    eco = Test1Ecosystem()
    @test EcoSISTEM._finalelapsed(eco, 12.0month_mean_duration,
                                  1.0month_mean_duration) ≈
          uconvert(s, 13.0month_mean_duration)

    # ...and it counts from where the clock already is, since `simulate!` does not reset it.
    EcoSISTEM._advanceclock!(eco, 5.0month_mean_duration)
    @test EcoSISTEM._finalelapsed(eco, 12.0month_mean_duration,
                                  1.0month_mean_duration) ≈
          uconvert(s, 18.0month_mean_duration)

    # That whole-step difference is enough to change the verdict on its own: a twelve-slice
    # climatology covers twelve months but not the thirteenth step such a run really takes.
    layer = _climatologylayer(atend = ErrorAtEnd())
    @test_nowarn EcoSISTEM._checkcoverage(layer, 12.0month_mean_duration)
    @test_throws ErrorException EcoSISTEM._checkcoverage(layer,
                                                         13.0month_mean_duration)
end

end
