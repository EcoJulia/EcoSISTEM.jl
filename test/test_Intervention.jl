# SPDX-License-Identifier: LGPL-3.0-or-later

module TestIntervention

using EcoSISTEM
# These are `public` rather than exported - a spec is what a user writes,
# and these are what it materialises into.
using EcoSISTEM: SteadyLayerChange, materialise
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Test
using Random
using DataFrames: DataFrame
using DimensionalData: DimensionalData, DimArray, X, Y
using Rasters
import Dates

include("rasterfixtures.jl")

# A small synthetic ecosystem on a 10 × 10 grid. A **fresh** `StudyArea` each time: an area is
# reusable by design, and sharing one between ecosystems is exactly the bug this file first exposed
# (see the regression test at the end). Any other keyword - an `epoch`, a `calendar` - goes to
# `build_ecosystem`.
function _eco(; abundance = 3000, seed = 1, kw...)
    area = StudyArea(extent = (10.0km, 10.0km), cellsize = 1.0km,
                     verbosity = :silent)
    env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                      supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                           axis = SolarRadiation),
                      area = area)
    spp = build_species(SpeciesPool(3),
                        tolerance = ToleranceSpec(Temperature, μ = 285.0K,
                                                  σ = 5.0K),
                        demand = DemandSpec(SolarRadiation, 1.0e9kJ / day),
                        seed = seed, verbosity = :silent)
    return build_ecosystem(BySupply(total = abundance), spp, env; seed = seed,
                           verbosity = :silent, kw...)
end
_active(eco) = count(parent(eco.habitat.active))

# An ecosystem with nothing placed on a projected grid, 8 rows by 9 columns of 2.5 km over central
# Scotland, so that a place can be named by its latitude and longitude and a shape laid over it.
function _placedeco()
    regime = _reg(_bngraster(EcoSISTEM.SyntheticData, fill(298.0K, 8, 9),
                             north = (640000.0:2500.0:657500.0) .* m),
                  axis = Temperature)
    env = GridHabitat(regime = regime,
                      supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                           axis = SolarRadiation),
                      area = StudyArea(regime = regime, verbosity = :silent))
    spp = build_species(SpeciesPool(2),
                        tolerance = ToleranceSpec(Temperature, μ = 298.0K,
                                                  σ = 5.0K),
                        demand = DemandSpec(SolarRadiation, 1.0e9kJ / day),
                        seed = 1, verbosity = :silent)
    return build_ecosystem(NoAbundance(), spp, env, seed = 1,
                           verbosity = :silent)
end

@testset "Schedules fire when they say" begin
    fires(sch, elapsed) = EcoSISTEM._fires(sch, elapsed, 1.0month_mean_duration)
    # Every step, but the start of a run is not a step: nothing has been taken yet to act over.
    @test !fires(EveryStep(), 0.0s)
    @test fires(EveryStep(), 100.0year)
    @test !fires(NeverScheduled(), 0.0s)
    @test !fires(NeverScheduled(), 100.0year)

    # A one-off schedule fires on the step that *reaches* its instant - the half-open window
    # `(elapsed - timestep, elapsed]`. Not an equality test: elapsed time accumulates as a float
    # and a run's steps need not land on the instant exactly, so `==` would silently never fire.
    at = AtTime(5.0month_mean_duration)
    @test !fires(at, 4.0month_mean_duration)
    @test fires(at, 5.0month_mean_duration)          # lands exactly on it
    @test !fires(at, 6.0month_mean_duration)
    # ...and it fires exactly once even when no step lands on it: 5 months falls inside the step that
    # covers months 4.5->5.5.
    @test EcoSISTEM._fires(at, 5.5month_mean_duration, 1.0month_mean_duration)
    @test !EcoSISTEM._fires(at, 4.5month_mean_duration, 1.0month_mean_duration)

    # ...and `AtTime(0)` fires at the start of the run, on the starting state, and not again on the
    # first step, whose window `(0, timestep]` excludes zero.
    start = AtTime(0.0s)
    @test fires(start, 0.0s)
    @test !EcoSISTEM._fires(start, 1.0month_mean_duration,
                            1.0month_mean_duration)
    @test fires(AtTime(-1.0month_mean_duration), 0.0s)

    many = AtTimes([2.0month_mean_duration, 7.0month_mean_duration])
    @test fires(many, 2.0month_mean_duration)
    @test fires(many, 7.0month_mean_duration)
    @test !fires(many, 4.0month_mean_duration)

    span = BetweenTimes(2.0month_mean_duration, 4.0month_mean_duration)
    @test !fires(span, 1.0month_mean_duration)
    @test fires(span, 2.0month_mean_duration)
    @test fires(span, 3.0month_mean_duration)
    @test fires(span, 4.0month_mean_duration)
    @test !fires(span, 5.0month_mean_duration)
    # A span acts over steps taken, so not at the start even when it opens there.
    @test !fires(BetweenTimes(0.0s, 4.0month_mean_duration), 0.0s)

    # Every whole multiple of an interval, on the step that reaches each - the answer `AtTimes`
    # gives over the explicit list, step for step, including where the steps do not land on the
    # multiples and where the interval and timestep are floats in different units.
    function firingsteps(schedule, timestep, nsteps)
        elapsed = 0.0s
        steps = Int[]
        for step in 1:nsteps
            elapsed += uconvert(s, float(timestep))
            EcoSISTEM._fires(schedule, elapsed, timestep) && push!(steps, step)
        end
        return steps
    end
    # The multiple at zero is the start of the run, which no step reaches.
    @test EcoSISTEM._fires(EveryInterval(3.0day), 0.0s, 1.0day)
    @test firingsteps(EveryInterval(3.0day), 1.0day, 10) == [3, 6, 9]
    for (interval, timestep, nsteps) in ((3.0day, 1.0day, 30),
        (10.0day, 3.0day, 40),
        (1.0year, 1.0month_mean_duration, 120),
        (12.0month_mean_duration,
        1.0month_mean_duration, 120))
        listed = AtTimes([k * interval
                          for k in 0:ceil(Int,
                                          ustrip(NoUnits,
                                                 nsteps * timestep / interval))])
        @test firingsteps(EveryInterval(interval), timestep, nsteps) ==
              firingsteps(listed, timestep, nsteps)
    end
    @test length(firingsteps(EveryInterval(1.0year), 1.0month_mean_duration,
                             121)) == 10
    @test_throws ArgumentError EveryInterval(0.0s)
    @test sprint(show, EveryInterval(3.0day)) == "EveryInterval(3.0 d)"
end

@testset "Date schedules fire on the step that reaches the date" begin
    epoch = Dates.DateTime(2000, 1, 1)
    # The steps on which a schedule fires over a run, its dates read through the ecosystem's epoch
    # and calendar.
    function datesteps(schedule, eco, timestep, nsteps)
        elapsed = 0.0s
        steps = Int[]
        for step in 1:nsteps
            elapsed += uconvert(s, float(timestep))
            EcoSISTEM._fires(schedule, eco, elapsed, timestep) &&
                push!(steps, step)
        end
        return steps
    end
    exact = _eco(epoch = epoch)
    counted = _eco(epoch = epoch, calendar = MeanMonths())

    # 1 March 2000 is day 60 of a leap year, so a daily run reaches it on step 60; dates may be
    # `Date`s or `DateTime`s, in any order.
    @test datesteps(AtDates([Dates.Date(2000, 3, 1)]), exact, 1.0day, 100) ==
          [60]
    @test datesteps(AtDates([Dates.Date(2000, 3, 1),
                                Dates.DateTime(2000, 1, 11)]), exact, 1.0day,
                    100) == [10, 60]

    # Every 1 July for three years at monthly steps: the real day under exact dates (day 182 of
    # 2000 is 5.98 mean months in), six mean months into each year when months are counted.
    july = EveryYear(month = 7)
    @test datesteps(july, exact, 1.0month_mean_duration, 36) == [6, 18, 30]
    @test datesteps(july, counted, 1.0month_mean_duration, 36) == [6, 18, 30]

    # 1 January 2001 is 366 days in, a little over twelve mean months, so the two calendars part by a
    # step there; the anniversary on the epoch itself fires at the start under both.
    @test datesteps(EveryYear(), exact, 1.0month_mean_duration, 13) == [13]
    @test datesteps(EveryYear(), counted, 1.0month_mean_duration, 13) == [12]
    @test EcoSISTEM._fires(EveryYear(), exact, 0.0s, 1.0month_mean_duration)
    @test EcoSISTEM._fires(EveryYear(), counted, 0.0s, 1.0month_mean_duration)

    # An anniversary before the epoch does not fire: a March start sees its first 1 January ten
    # months later, not on step one.
    later = _eco(epoch = Dates.DateTime(2000, 3, 1))
    # At the start every earlier instant is in the window, so this is where one before the epoch
    # would fire if it were not refused.
    @test !EcoSISTEM._fires(EveryYear(month = 1), later, 0.0s,
                            1.0month_mean_duration)
    @test datesteps(EveryYear(month = 1), later, 1.0month_mean_duration, 12) ==
          [11]

    # A date that is not in every year, or not in any, is refused where it is written.
    @test_throws ArgumentError EveryYear(month = 2, day = 29)
    @test_throws ArgumentError EveryYear(month = 4, day = 31)
    @test_throws ArgumentError EveryYear(month = 13)
    @test sprint(show, EveryYear(month = 7)) == "EveryYear(month = 7, day = 1)"
    @test occursin("2000-03-01",
                   sprint(show, AtDates([Dates.Date(2000, 3, 1)])))

    # A run with no epoch has no dates, so a date schedule is refused before its first step...
    clearance = Intervention(AtDates([Dates.Date(2000, 3, 1)]), AllCells(),
                             Deactivate())
    @test_throws "epoch" simulate!(_eco(), 2.0month_mean_duration,
                                   1.0month_mean_duration,
                                   intervention = clearance)
    # A date inside a step longer than a day would be acted on late, by a different amount each year,
    # so it is refused too, naming the date it would have been acted on. 1 January 2001 is 366 days
    # in, which no mean-month or mean-year step ends on under exact dates.
    refusal(run) =
        try
            run()
            ""
        catch err
            sprint(showerror, err)
        end
    yearly = Intervention(EveryYear(), AllCells(), RemoveAbundance(1, 1))
    monthly = refusal(() -> simulate!(_eco(epoch = epoch),
                                      13.0month_mean_duration,
                                      1.0month_mean_duration,
                                      intervention = yearly))
    @test occursin("2001-01-01", monthly) && occursin("2001-01-30", monthly)
    @test occursin("2001-01-01",
                   refusal(() -> simulate!(_eco(epoch = epoch), 3.0year,
                                           1.0year, intervention = yearly)))
    # Counting months, every twelfth mean month and every mean year ends on 1 January.
    @test isnothing(simulate!(_eco(epoch = epoch, calendar = MeanMonths()),
                              13.0month_mean_duration, 1.0month_mean_duration,
                              intervention = yearly))
    @test isnothing(simulate!(_eco(epoch = epoch, calendar = MeanMonths()),
                              3.0year, 1.0year, intervention = yearly))
    # A step of a day shows every date on its own day, under either calendar.
    for calendar in (ExactDates(), MeanMonths())
        @test isnothing(simulate!(_eco(epoch = epoch, calendar = calendar),
                                  400.0day, 1.0day, intervention = yearly))
    end
    # No monthly step ends on the fifteenth under either calendar; counting months, one ends on the
    # first.
    midmonth = Intervention(AtDates([Dates.Date(2000, 3, 15)]), AllCells(),
                            RemoveAbundance(1, 1))
    for calendar in (ExactDates(), MeanMonths())
        @test occursin("2000-03-15",
                       refusal(() -> simulate!(_eco(epoch = epoch,
                                                    calendar = calendar),
                                               6.0month_mean_duration,
                                               1.0month_mean_duration,
                                               intervention = midmonth)))
    end
    firstofmonth = Intervention(AtDates([Dates.Date(2000, 3, 1)]), AllCells(),
                                RemoveAbundance(1, 1))
    @test isnothing(simulate!(_eco(epoch = epoch, calendar = MeanMonths()),
                              6.0month_mean_duration, 1.0month_mean_duration,
                              intervention = firstofmonth))

    # ...and on a run whose steps end on its date it acts, once, on that step.
    eco = _eco(epoch = epoch, calendar = MeanMonths())
    simulate!(eco, 7.0month_mean_duration, 1.0month_mean_duration,
              intervention = Intervention(EveryYear(month = 7), RandomCells(20),
                                          Deactivate()))
    @test _active(eco) == 80
end

@testset "Regions resolve to the cells they name" begin
    eco = _eco()
    # The timestep is what turns a *rate*-valued count into a probability, so region resolution
    # needs it even for the regions that ignore it.
    function cells(r)
        return EcoSISTEM._regioncells(r, eco, Random.Xoshiro(1),
                                      1.0month_mean_duration)
    end
    @test length(cells(AllCells())) == 100
    @test length(cells(ActiveCells())) == 100

    # Deactivating half makes the two differ - which is the whole reason both exist.
    parent(eco.habitat.active)[1:50] .= false
    @test length(cells(AllCells())) == 100
    @test length(cells(ActiveCells())) == 50

    mask = falses(10, 10)
    mask[1:7] .= true
    @test length(cells(CellMask(mask))) == 7
    @test_throws ErrorException cells(CellMask(falses(3, 3)))

    # Random regions draw only from **active** cells, and never take more than exist.
    @test length(cells(RandomCells(20))) == 20
    @test all(in(cells(ActiveCells())), cells(RandomCells(20)))
    @test length(cells(RandomCells(500))) == 50      # capped at what is available
    @test length(cells(SpreadingCells(20))) == 20
    @test length(cells(SpreadingCells(500))) == 50

    # A **rate** draws binomially over the step instead of taking a fixed count - the process the
    # v0.4.0 scenarios had (`jbinom(1, npos, rate)`), which a fixed count cannot express.
    @test length(cells(RandomCells(0.0 / year))) == 0        # certain not to happen
    @test length(cells(RandomCells(1.0e6 / year))) == 50     # clamped at certainty, and at supply
    drawn = cells(RandomCells(0.5 / year))
    @test 0 <= length(drawn) <= 50
    # ...and the draw is reproducible, being keyed on the step rather than the global RNG.
    @test cells(RandomCells(0.5 / year)) == drawn
end

@testset "CellsAt names cells by place, moving a place out of an inactive cell within snap" begin
    eco = _placedeco()
    cells(r) = EcoSISTEM._regioncells(r, eco, nothing, nothing)
    grid = LinearIndices(parent(eco.habitat.active))
    here = LatLong(55.7°, -4.3°)
    cell = EcoSISTEM.getcellat(eco, here)
    @test all(2 .<= Tuple(cell) .<= (7, 8))
    # A place, a grid key, and the place's own cell by key: each cell once, in grid order, whatever
    # order the places come in.
    @test cells(CellsAt(here, CartesianIndex(1, 1), cell)) ==
          [grid[1, 1], grid[cell]]
    @test cells(CellsAt(cell, here, CartesianIndex(1, 1))) ==
          cells(CellsAt(here, CartesianIndex(1, 1)))
    # A place in an inactive cell is that cell, as a mask's are, unless `snap` reaches an active one:
    # its four neighbours are equally near, and the lower `(y, x)` is taken. No neighbour is within
    # less than a cell.
    parent(eco.habitat.active)[cell] = false
    @test cells(CellsAt(here)) == [grid[cell]]
    @test cells(CellsAt(here, snap = 3.0km)) ==
          [grid[cell - CartesianIndex(1, 0)]]
    @test cells(CellsAt(here, snap = 2.0km)) == [grid[cell]]
    # So a site can be reactivated by naming it.
    simulate!(eco, 1month_mean_duration, 1month_mean_duration,
              intervention = Intervention(AtTime(0.0year), CellsAt(here),
                                          Reactivate()))
    @test parent(eco.habitat.active)[cell]

    # A place off the grid is refused before the run's first step, by position or by key.
    for (region, message) in ((CellsAt(LatLong(57.0°, -4.3°)),
                               "which is off the grid"),
                              (CellsAt(CartesianIndex(9, 1)),
                               "the cell (9, 1), which is off the 8 by 9 grid"))
        fresh = _placedeco()
        @test_throws message simulate!(fresh, 2month_mean_duration,
                                       1month_mean_duration,
                                       intervention = Intervention(AtTime(1.0year),
                                                                   region,
                                                                   Deactivate()))
        @test iszero(EcoSISTEM.simulationtime(fresh))
    end
    @test_throws "at least one place" CellsAt()
    @test_throws "a finite, non-negative length; got -1 km" CellsAt(here,
                                                                    snap = -1km)
    @test_throws TypeError CellsAt(here, snap = 3)
    @test_throws MethodError CellsAt((1, 2))
    @test sprint(show, CellsAt(CartesianIndex(1, 2), snap = 5km)) ==
          "CellsAt(CartesianIndex(1, 2), snap = 5 km)"
    @test sprint(show, CellsAt(here)) == "CellsAt(" * sprint(show, here) * ")"
end

@testset "A shape is a region, as within takes one: the cells it covers, found once an ecosystem" begin
    path = joinpath(mktempdir(), "square.geojson")
    write(path,
          """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
             "geometry":{"type":"Polygon","coordinates":[[[-4.39,55.675],[-4.27,55.675],
             [-4.27,55.745],[-4.39,55.745],[-4.39,55.675]]]}}]}""")
    square = ShapeSpec(path)
    eco = _placedeco()
    cells(r) = EcoSISTEM._regioncells(EcoSISTEM._asregion(r), eco, nothing,
                                      nothing)
    # A bare shape is the cells it covers at least half of, and a rule picks others: every cell it
    # touches, or only those wholly inside.
    covered = materialise(ShapeCoverage(square, axis = EcoSISTEM.NicheAxis),
                          StudyArea(eco.habitat))
    half = cells(square)
    @test half == findall(vec(parent(covered.matrix) .>= 0.5))
    @test cells(ShapeMaskSpec(square)) == half
    touched = cells(ShapeMaskSpec(square, AnyOverlap()))
    inside = cells(ShapeMaskSpec(square, FullyWithin()))
    @test issubset(inside, half) && issubset(half, touched)
    @test 0 < length(inside) < length(half) < length(touched) < 72
    @test_throws "`WholeBoundingBox()` says where a study area's grid is cut" Intervention(AtTime(1.0year),
                                                                                           ShapeMaskSpec(square,
                                                                                                         WholeBoundingBox()),
                                                                                           Deactivate())

    # In a run, placed before the first step and acted on when it fires, the records of the shape's
    # file then joining the ecosystem's.
    iv = Intervention(AtTime(1month_mean_duration), square, Deactivate())
    @test sprint(show, iv) ==
          "Intervention(" * sprint(show, AtTime(1month_mean_duration)) * ", " *
          sprint(show, square) * ", " * sprint(show, Deactivate()) * ")"
    before = length(provenance(eco).inputs)
    simulate!(eco, 2month_mean_duration, 1month_mean_duration,
              intervention = iv)
    @test findall(!, vec(parent(eco.habitat.active))) == half
    @test length(provenance(eco).inputs) == before + 1
    @test any(r -> r.path == "square.geojson", provenance(eco).inputs)
    # Found once for the ecosystem, and again for another.
    found = iv.region.found[]
    @test found.owner.value === eco
    simulate!(eco, 1month_mean_duration, 1month_mean_duration,
              intervention = iv)
    @test iv.region.found[] === found
    other = _placedeco()
    simulate!(other, 2month_mean_duration, 1month_mean_duration,
              intervention = iv)
    @test iv.region.found[].owner.value === other
    @test findall(!, vec(parent(other.habitat.active))) == half
    # An intervention that never fires records nothing.
    never = _placedeco()
    simulate!(never, 1month_mean_duration, 1month_mean_duration,
              intervention = Intervention(AtTime(10.0year), square,
                                          Deactivate()))
    @test length(provenance(never).inputs) == before
    # A grid with no coordinate reference system cannot place a shape, refused before the first step.
    synthetic = _eco()
    @test_throws "needs a grid with a coordinate reference system" simulate!(synthetic,
                                                                             1month_mean_duration,
                                                                             1month_mean_duration,
                                                                             intervention = Intervention(AtTime(1.0year),
                                                                                                         square,
                                                                                                         Deactivate()))
    @test iszero(EcoSISTEM.simulationtime(synthetic))
end

@testset "The six operations do what they say" begin
    # Deactivate / Reactivate.
    eco = _eco()
    simulate!(eco, 2.0month_mean_duration, 1.0month_mean_duration,
              intervention = Intervention(AtTime(1.0month_mean_duration),
                                          RandomCells(20), Deactivate()))
    @test _active(eco) == 80
    simulate!(eco, 1.0month_mean_duration, 1.0month_mean_duration,
              intervention = Intervention(EveryStep(), AllCells(),
                                          Reactivate()))
    @test _active(eco) == 100

    # AddAbundance / RemoveAbundance, on a frozen population so the demography cannot confuse the
    # arithmetic - `NoGrowth` zeroes birth and death.
    eco2 = _eco()
    before = sum(eco2.abundances.matrix[1, :])
    EcoSISTEM.applyinterventions!(eco2,
                                  Intervention(EveryStep(), AllCells(),
                                               AddAbundance(1, 5)),
                                  1.0month_mean_duration,
                                  1.0month_mean_duration, 1)
    @test sum(eco2.abundances.matrix[1, :]) == before + 5 * 100

    # Removal is clamped per cell: taking more than are present cannot go negative.
    EcoSISTEM.applyinterventions!(eco2,
                                  Intervention(EveryStep(), AllCells(),
                                               RemoveAbundance(1, 10^9)),
                                  1.0month_mean_duration,
                                  1.0month_mean_duration, 1)
    @test all(iszero, eco2.abundances.matrix[1, :])
    @test all(>=(0), eco2.abundances.matrix)

    # Species may be named as well as indexed.
    eco3 = _eco()
    name = first(eco3.spplist.names)
    EcoSISTEM.applyinterventions!(eco3,
                                  Intervention(EveryStep(), AllCells(),
                                               AddAbundance(name, 3)),
                                  1.0month_mean_duration,
                                  1.0month_mean_duration, 1)
    @test sum(eco3.abundances.matrix[1, :]) ==
          sum(_eco().abundances.matrix[1, :]) + 300
    @test_throws ErrorException EcoSISTEM.applyinterventions!(eco3,
                                                              Intervention(EveryStep(),
                                                                           AllCells(),
                                                                           AddAbundance("nosuch",
                                                                                        1)),
                                                              1.0month_mean_duration,
                                                              1.0month_mean_duration,
                                                              1)

    # `SetLandCover` is the ONLY direct matrix write, and only on a categorical layer - a
    # continuous layer's values belong to its change rule.
    eco4 = _eco()
    @test_throws ErrorException EcoSISTEM.applyinterventions!(eco4,
                                                              Intervention(EveryStep(),
                                                                           AllCells(),
                                                                           SetLandCover(:open_water)),
                                                              1.0month_mean_duration,
                                                              1.0month_mean_duration,
                                                              1)
end

@testset "SetChange installs a rule, and bites the same step" begin
    # The ordering fix: interventions run after the dynamics and *before* the layer update, so a
    # change installed on a step is applied on that step rather than one later.
    eco = _eco()
    @test eco.habitat.regime.change isa NoLayerChange
    simulate!(eco, 1.0month_mean_duration, 1.0month_mean_duration,
              intervention = Intervention(AtTime(0.0s), AllCells(),
                                          SetChange(nothing,
                                                    IncrementBy(1.0K / year))))
    @test eco.habitat.regime.change isa SteadyLayerChange
    # Warmed already - not still sitting at its initial 285 K waiting for the next step.
    @test all(>(285.0K), eco.habitat.regime.matrix)

    # A named target addresses a sub-layer of a collection; an unknown name says so.
    @test_throws ErrorException EcoSISTEM._targetlayer(eco, :nosuchlayer)
end

@testset "An intervention due at the start acts before the first step's dynamics" begin
    # A species arriving at the start has been through the first step's births and deaths by the end
    # of it; one arriving on that step's own instant comes after them, and is still exactly as many
    # as arrived.
    arrival(time) = Intervention(AtTime(time), AllCells(),
                                 AddSpecies(abundance = 500))
    atstart = _eco()
    simulate!(atstart, 1.0month_mean_duration, 1.0month_mean_duration,
              intervention = arrival(0.0s))
    @test size(atstart.abundances.matrix, 1) == 4
    @test sum(atstart.abundances.matrix[4, :]) != 500
    onstep = _eco()
    simulate!(onstep, 1.0month_mean_duration, 1.0month_mean_duration,
              intervention = arrival(1.0month_mean_duration))
    @test size(onstep.abundances.matrix, 1) == 4
    @test sum(onstep.abundances.matrix[4, :]) == 500
end

@testset "Interventions are reproducible and order-stable" begin
    # Counter-based per step (seeded from `(seed, :intervention, k, step)`), so a run replays exactly and
    # every MPI rank and thread computes the same selection without communicating. v0.4.0's
    # `RandHabitatLoss!` drew from the **global** RNG and could not be reproduced at all.
    run() = begin
        e = _eco()
        simulate!(e, 5.0month_mean_duration, 1.0month_mean_duration,
                  intervention = Intervention(AtTime(2.0month_mean_duration),
                                              RandomCells(20), Deactivate()))
        (findall(vec(parent(e.habitat.active))), copy(e.abundances.matrix))
    end
    a, b = run(), run()
    @test a[1] == b[1]
    @test a[2] == b[2]

    # A set applies its members in the order written.
    eco = _eco()
    set = InterventionSet(Intervention(EveryStep(), AllCells(), Deactivate()),
                          Intervention(EveryStep(), AllCells(), Reactivate()))
    simulate!(eco, 1.0month_mean_duration, 1.0month_mean_duration,
              intervention = set)
    @test _active(eco) == 100                 # deactivated, then reactivated
    reversed = InterventionSet(Intervention(EveryStep(), AllCells(),
                                            Reactivate()),
                               Intervention(EveryStep(), AllCells(),
                                            Deactivate()))
    eco2 = _eco()
    simulate!(eco2, 1.0month_mean_duration, 1.0month_mean_duration,
              intervention = reversed)
    @test _active(eco2) == 0                  # reactivated, then deactivated
end

@testset " a StudyArea is reusable: environments must not share `active`" begin
    # Two environments built from one `StudyArea` must not share the **same** `active` array
    # object, which is what the synthetic path does if it passes `area.report.active` straight
    # through without copying. Deactivating cells in one ecosystem then silently deactivates them in
    # every other - and in the area itself, so the next `GridHabitat` inherits the damage. Harmless
    # until something
    # mutated `active`, which is why it survived until interventions arrived.
    area = StudyArea(extent = (10.0km, 10.0km), cellsize = 1.0km,
                     verbosity = :silent)
    env() = GridHabitat(regime = UniformSpec(285.0K,
                                             axis = Temperature),
                        supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                             axis = SolarRadiation),
                        area = area)
    a, b = env(), env()
    @test !(parent(a.active) === parent(b.active))
    parent(a.active)[1:20] .= false
    @test count(parent(b.active)) == 100
    @test count(parent(env().active)) == 100   # ...and the area itself is undamaged
end

@testset "an operation's species may be a Symbol" begin
    eco = _eco()
    before = sum(eco.abundances.matrix[1, :])
    EcoSISTEM.applyinterventions!(eco,
                                  Intervention(EveryStep(), AllCells(),
                                               AddAbundance(Symbol(first(eco.spplist.names)),
                                                            2)),
                                  1.0month_mean_duration,
                                  1.0month_mean_duration, 1)
    @test sum(eco.abundances.matrix[1, :]) == before + 200
end

@testset "AddAbundance adds a placement's individuals over the region at each firing" begin
    month = 1.0month_mean_duration
    fire!(eco, iv; step = 1) = EcoSISTEM.applyinterventions!(eco, iv, month,
                                                             month, step)
    mask = falses(10, 10)
    mask[2:3, 4:6] .= true
    linear = LinearIndices((10, 10))
    added = Intervention(AtTime(month), CellMask(mask),
                         AddAbundance(BySupply(total = 600)))
    eco = _eco(abundance = 0)
    eco.habitat.active[2, 4] = false
    fire!(eco, added)
    @test sum(eco.abundances.matrix) == 600
    @test sum(eco.abundances.matrix[:, linear[.!mask]]) == 0
    @test sum(eco.abundances.matrix[:, linear[2, 4]]) == 0
    # The ecosystem keeps its own placement; the firing's streams are its own, step by step.
    @test sprint(show, eco.placement) == "BySupply(total = 0)"
    again = _eco(abundance = 0)
    again.habitat.active[2, 4] = false
    fire!(again, added)
    @test again.abundances.matrix == eco.abundances.matrix
    later = _eco(abundance = 0)
    later.habitat.active[2, 4] = false
    fire!(later, added, step = 2)
    @test later.abundances.matrix != eco.abundances.matrix
    # A placement draws nothing from the stream the region and the other operations draw from.
    plain = Intervention(AtTime(month), RandomCells(5), AddAbundance("1", 3))
    sown = Intervention(AtTime(month), RandomCells(5),
                        AddAbundance(ToCapacity("2")), AddAbundance("1", 3))
    one, two = _eco(abundance = 0), _eco(abundance = 0)
    fire!(one, plain)
    fire!(two, sown)
    @test one.abundances.matrix[1, :] == two.abundances.matrix[1, :]
    @test findall(>(0), two.abundances.matrix[2, :]) ==
          findall(>(0), one.abundances.matrix[1, :])
    # What it refuses: a dated table, and what is not a placement.
    timed = (_, _, _) -> (species = ["1"], cell = [CartesianIndex(1, 1)],
                          time = [0.0s])
    @test_throws "its table is dated" fire!(_eco(abundance = 0),
                                            Intervention(AtTime(month),
                                                         AllCells(),
                                                         AddAbundance(timed)))
    @test_throws "a placement says how many individuals start where" AddAbundance(3)
    @test sprint(show, AddAbundance(ToCapacity("wheat"))) ==
          "AddAbundance(ToCapacity(\"wheat\"))"
end

@testset "AddAbundance serves a dated placement's rows as the clock reaches them" begin
    month = 1.0month_mean_duration
    fire!(eco, iv, elapsed, step) = EcoSISTEM.applyinterventions!(eco, iv,
                                                                  elapsed,
                                                                  month,
                                                                  step)
    linear = LinearIndices((10, 10))
    table = (species = ["1", "2", "1", "3"],
             cell = [CartesianIndex(1, 1), CartesianIndex(2, 2),
                 CartesianIndex(3, 3), CartesianIndex(4, 4)],
             count = [5, 7, 3, 4],
             time = [-1.0, 0.0, 1.0, 2.5] .* month)
    served = Intervention(EveryStep(), AllCells(),
                          AddAbundance(FromTable(table)))
    eco = _eco(abundance = 0)
    # The start takes every row at or before it; each step, those in its window.
    fire!(eco, served, 0.0s, 0)
    @test eco.abundances.matrix[1, linear[1, 1]] == 5
    @test eco.abundances.matrix[2, linear[2, 2]] == 7
    @test sum(eco.abundances.matrix) == 12
    fire!(eco, served, month, 1)
    @test eco.abundances.matrix[1, linear[3, 3]] == 3
    fire!(eco, served, 2month, 2)
    @test sum(eco.abundances.matrix) == 15
    rows = only(served.operations).served[].rows
    fire!(eco, served, 3month, 3)
    @test eco.abundances.matrix[3, linear[4, 4]] == 4
    @test only(served.operations).served[].rows === rows
    # The region filters each step's rows.
    mask = trues(10, 10)
    mask[2, 2] = false
    inregion = Intervention(EveryStep(), CellMask(mask),
                            AddAbundance(FromTable(table)))
    eco = _eco(abundance = 0)
    fire!(eco, inregion, 0.0s, 0)
    @test sum(eco.abundances.matrix) == 5
    # A run's end notes what it did not add: rows after it, and rows due in cells deactivated since.
    eco = _eco(abundance = 0)
    eco.habitat.active[3, 3] = false
    @test_logs((:info,
                r"^3 individuals due in cells deactivated since the ecosystem was built were not added\.$"),
               (:info,
                r"^4 individuals dated after the run's end were not added\.$"),
               match_mode = :any,
               simulate!(eco, 2month, month, intervention = served))
    # A run continuing another adds only what falls due as it runs.
    eco = _eco(abundance = 0)
    simulate!(eco, 1month, month)
    @test_logs((:info,
                r"^15 individuals dated at or before the run's start were not added"),
               match_mode = :any,
               simulate!(eco, 2month, month, intervention = served))
    # What is refused: a dated placement scheduled otherwise or sharing its intervention, and a sum
    # mixing dated and undated placements.
    @test_throws "carries its own times" Intervention(AtTime(month), AllCells(),
                                                      AddAbundance(FromTable(table)))
    @test_throws "carries its own times" Intervention(EveryStep(), AllCells(),
                                                      AddAbundance(FromTable(table)),
                                                      Deactivate())
    @test_throws "mixes dated and undated ones" AddAbundance(FromTable(table) +
                                                             BySupply())
end

@testset "A dated placement follows the run's calendar, any steps, and any ecosystem" begin
    month = 1.0month_mean_duration
    fire!(eco, iv, elapsed, step;
          timestep = month) = EcoSISTEM.applyinterventions!(eco,
                                                            iv,
                                                            elapsed,
                                                            timestep,
                                                            step)
    # Stepped on the run's own clock, each first of the month lands on its own step.
    epoch = Dates.Date(2000, 1, 1)
    dates = [epoch + Dates.Month(k) for k in 0:23]
    monthly = Intervention(EveryStep(), AllCells(),
                           AddAbundance(FromTable((species = fill("1", 24),
                                                   cell = fill(CartesianIndex(10,
                                                                              1),
                                                               24),
                                                   count = fill(7, 24),
                                                   time = dates))))
    eco = _eco(abundance = 0, epoch = epoch, calendar = MeanMonths())
    added = map(0:24) do k
        before = eco.abundances.matrix[1, 10]
        fire!(eco, monthly, EcoSISTEM.simulationtime(eco), k)
        EcoSISTEM._advanceclock!(eco, month)
        return (eco.abundances.matrix[1, 10] - before) / 7
    end
    @test added == [ones(24); 0]
    @test_throws "no epoch" fire!(_eco(abundance = 0), monthly, 0.0s, 0)
    # One two-month step adds what two one-month steps do, and any ecosystem gets the same rows.
    table = (species = ["1", "2", "3"],
             cell = [CartesianIndex(1, 1), CartesianIndex(2, 2),
                 CartesianIndex(3, 3)],
             count = [3, 4, 5], time = [0.0, 1.0, 2.0] .* month)
    timed = Intervention(EveryStep(), AllCells(),
                         AddAbundance(FromTable(table)))
    short, long = _eco(abundance = 0), _eco(abundance = 0)
    for (k, t) in enumerate((0.0, 1.0, 2.0) .* month)
        fire!(short, timed, t, k - 1)
    end
    fire!(long, timed, 0.0s, 0, timestep = 2month)
    fire!(long, timed, 2month, 1, timestep = 2month)
    @test long.abundances.matrix == short.abundances.matrix
    @test sum(short.abundances.matrix) == 12
    # The starting state, as a run's callback first sees it, holds the rows at the start.
    eco = _eco()
    before = sum(eco.abundances.matrix[1, :])
    atstart = Ref(0)
    simulate!(eco, month, month,
              intervention = Intervention(EveryStep(), AllCells(),
                                          AddAbundance(FromTable((species = ["1"],
                                                                  cell = [CartesianIndex(1,
                                                                                         1)],
                                                                  count = [40],
                                                                  time = [0.0s]))))
              ) do occurrence
        occurrence.count == 1 && (atstart[] = sum(eco.abundances.matrix[1, :]))
        return nothing
    end
    @test atstart[] == before + 40
end

@testset "ToCapacity at a firing tops each cell up, counting what is there" begin
    month = 1.0month_mean_duration
    fire!(eco, iv) = EcoSISTEM.applyinterventions!(eco, iv, month, month, 1)
    # A cell holds 100 of any species: 1e11 kJ a day supplied, 1e9 demanded by each.
    eco = _eco(abundance = 3000)
    before = copy(eco.abundances.matrix)
    fire!(eco,
          Intervention(AtTime(month), AllCells(),
                       AddAbundance(ToCapacity("2"))))
    @test eco.abundances.matrix[[1, 3], :] == before[[1, 3], :]
    @test all(==(100), sum(eco.abundances.matrix, dims = 1))
    # To half capacity: a cell holding more already gets none.
    half = _eco(abundance = 5000)
    present = vec(sum(half.abundances.matrix, dims = 1))
    @test any(>(50), present) && any(<(50), present)
    fire!(half,
          Intervention(AtTime(month), AllCells(),
                       AddAbundance(ToCapacity("2", fraction = 0.5))))
    @test vec(sum(half.abundances.matrix, dims = 1)) == max.(present, 50)
    # A cell already full gets none.
    full = _eco(abundance = 0)
    repopulate!(full, ToCapacity("1"), verbosity = :silent)
    @test all(==(100), full.abundances.matrix[1, :])
    fire!(full,
          Intervention(AtTime(month), AllCells(),
                       AddAbundance(ToCapacity("2"))))
    @test sum(full.abundances.matrix[2, :]) == 0
end

end
