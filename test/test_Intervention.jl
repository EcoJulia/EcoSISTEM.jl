# SPDX-License-Identifier: LGPL-3.0-or-later

module TestIntervention

using EcoSISTEM
# `[C7-VIS]` C: these are `public` rather than exported - a spec is what a user writes,
# and these are what it materialises into.
using EcoSISTEM: SteadyLayerChange
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Test
using Random
import Dates

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
    spp = build_species(3, tolerance = (285.0K, 5.0K),
                        toleranceaxis = Temperature, demand = 1.0e9kJ / day,
                        demandaxis = SolarRadiation,
                        abundance = abundance, seed = seed)
    return build_ecosystem(spp, env; seed = seed, kw...)
end
_active(eco) = count(parent(eco.habitat.active))

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

end
