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

# A small synthetic ecosystem on a 10 × 10 grid. A **fresh** `StudyArea` each time: an area is
# reusable by design, and sharing one between ecosystems is exactly the bug this file first exposed
# (see the regression test at the end).
function _eco(; abundance = 3000, seed = 1)
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
    return build_ecosystem(spp, env, seed = seed)
end
_active(eco) = count(parent(eco.habitat.active))

@testset "Schedules fire when they say" begin
    fires(sch, elapsed) = EcoSISTEM._fires(sch, elapsed, 1.0month_mean_duration)
    @test fires(EveryStep(), 0.0s)
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

    # ...and `AtTime(0)` fires on the FIRST step, not never. The first window is `(0, timestep]`,
    # which excludes zero, so without a closed lower end on step one a schedule at the start would
    # silently never fire at all.
    start = AtTime(0.0s)
    @test EcoSISTEM._fires(start, 1.0month_mean_duration,
                           1.0month_mean_duration)
    @test !EcoSISTEM._fires(start, 2.0month_mean_duration,
                            1.0month_mean_duration)

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
                                  0.0s, 1.0month_mean_duration, 0)
    @test sum(eco2.abundances.matrix[1, :]) == before + 5 * 100

    # Removal is clamped per cell: taking more than are present cannot go negative.
    EcoSISTEM.applyinterventions!(eco2,
                                  Intervention(EveryStep(), AllCells(),
                                               RemoveAbundance(1, 10^9)),
                                  0.0s, 1.0month_mean_duration, 0)
    @test all(iszero, eco2.abundances.matrix[1, :])
    @test all(>=(0), eco2.abundances.matrix)

    # Species may be named as well as indexed.
    eco3 = _eco()
    name = first(eco3.spplist.names)
    EcoSISTEM.applyinterventions!(eco3,
                                  Intervention(EveryStep(), AllCells(),
                                               AddAbundance(name, 3)),
                                  0.0s, 1.0month_mean_duration, 0)
    @test sum(eco3.abundances.matrix[1, :]) ==
          sum(_eco().abundances.matrix[1, :]) + 300
    @test_throws ErrorException EcoSISTEM.applyinterventions!(eco3,
                                                              Intervention(EveryStep(),
                                                                           AllCells(),
                                                                           AddAbundance("nosuch",
                                                                                        1)),
                                                              0.0s,
                                                              1.0month_mean_duration,
                                                              0)

    # `SetLandCover` is the ONLY direct matrix write, and only on a categorical layer - a
    # continuous layer's values belong to its change rule.
    eco4 = _eco()
    @test_throws ErrorException EcoSISTEM.applyinterventions!(eco4,
                                                              Intervention(EveryStep(),
                                                                           AllCells(),
                                                                           SetLandCover(:open_water)),
                                                              0.0s,
                                                              1.0month_mean_duration,
                                                              0)
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

@testset "Interventions are reproducible and order-stable" begin
    # Counter-based per step (`hash((seed, :intervention, k, step))`), so a run replays exactly and
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
