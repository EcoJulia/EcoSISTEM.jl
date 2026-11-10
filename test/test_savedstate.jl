# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/savedstate.jl`: `write_state` saving the state a run has reached, and
# `build_ecosystem(species, habitat; from)` starting an ecosystem from it - exactly as the saved run
# would have gone on, on its own seed and with the change rules its interventions installed in force
# again; with fresh streams on another; with its clock at zero; with more species than it saved; and
# refused on another grid, species, epoch or calendar.
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_savedstate.jl"])'
module TestSavedState

using EcoSISTEM
using EcoSISTEM: simulationtime
using EcoSISTEM.Units
using DataFrames: DataFrame
using Dates: DateTime
using DimensionalData: DimArray, Ti, X, Y
using DimensionalData.Lookups: NoLookup
using Diversity
using Distributions
using JLD2: jldsave
using Random: Random
using Test
using TOML
using Unitful, Unitful.DefaultSymbols

include("TestCases.jl")

const STEP = 1.0month_mean_duration

# A 5 by 7 grid warming steadily, whose temperature therefore depends on the run's history.
function _habitat(; cells = (5.0km, 7.0km))
    area = StudyArea(extent = cells, cellsize = 1.0km, verbosity = :silent)
    return GridHabitat(regime = Varying(UniformSpec(10.0°C, axis = Temperature),
                                        IncrementBy(0.5K / year)),
                       supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                            axis = SolarRadiation),
                       area = area)
end

# Four species, the fourth left out of the run, from a pool of `names`.
function _species(names = ["1", "2", "3", "4"])
    pool = SpeciesPool(names)
    "4" in names && deactivate_species!(pool, "4", verbosity = :silent)
    return build_species(pool,
                         tolerance = ToleranceSpec(Temperature, μ = 10.0°C,
                                                   σ = 3.0K),
                         demand = DemandSpec(SolarRadiation, 1.0e9kJ / day),
                         seed = 1, verbosity = :silent)
end

# A cell destroyed and the fourth species brought in, during the first year.
const EVENTS = InterventionSet(Intervention(AtTime(6 * STEP),
                                            CellsAt(CartesianIndex(2, 2)),
                                            Deactivate()),
                               Intervention(AtTime(9 * STEP), AllCells(),
                                            AddSpecies("4",
                                                       BySupply(total = 50))))

# An ecosystem run for `months` from the start, with the events unless given other interventions.
function _run(months; seed = 7, intervention = EVENTS, kw...)
    eco = build_ecosystem(BySupply(total = 2000), _species(), _habitat(),
                          seed = seed, verbosity = :silent; kw...)
    simulate!(eco, months * STEP, STEP, intervention = intervention)
    return eco
end

# A rule of steady warming installed on the temperature layer.
_rule(rate) = SetChange(nothing, IncrementBy(rate))

# A faster warming from the first month.
const ONERULE = Intervention(AtTime(STEP), AllCells(), _rule(5.0K / year))

# Two rules due in the third month, at instants named the other way round from the set's order.
const TWORULES = InterventionSet(Intervention(AtTime(2.5 * STEP), AllCells(),
                                              _rule(5.0K / year)),
                                 Intervention(AtTime(2.2 * STEP), AllCells(),
                                              _rule(-1.0K / year)))

# The state of a year's run, saved in a directory of its own.
function _saved(; kw...)
    path = joinpath(mktempdir(), "state.jld2")
    return write_state(path, _run(12; kw...))
end

@testset "a run started from a state goes on as the saved run would have" begin
    straight = _run(24)
    path = _saved()
    @test isfile(path)
    @test isfile(joinpath(dirname(path), "state.provenance.toml"))
    resumed = build_ecosystem(_species(), _habitat(), from = path,
                              verbosity = :silent)
    # What the state carries: the clock, the species brought in, the destroyed cell, and the
    # temperature the steady warming had reached, which no clock could recompute.
    @test simulationtime(resumed) ≈ 12 * STEP
    @test resumed.spplist.traits.active == [true, true, true, true]
    @test !parent(resumed.habitat.active)[2, 2]
    @test count(parent(resumed.habitat.active)) == 34
    @test parent(resumed.habitat.regime.matrix) ≈
          parent(_run(12).habitat.regime.matrix)
    @test any(r -> r.role === :state && r.path == "state.jld2",
              resumed.inputs)
    simulate!(resumed, 12 * STEP, STEP, intervention = EVENTS)
    @test resumed.abundances.matrix == straight.abundances.matrix
    @test parent(resumed.habitat.regime.matrix) ≈
          parent(straight.habitat.regime.matrix)
    # Its record says where it is, and starts a run as the path does.
    record = EcoSISTEM.ProvenanceRecord(role = :state, source = "saved state",
                                        path = "state.jld2", location = path)
    again = build_ecosystem(_species(), _habitat(), from = record,
                            verbosity = :silent)
    simulate!(again, 12 * STEP, STEP, intervention = EVENTS)
    @test again.abundances.matrix == straight.abundances.matrix
end

@testset "another seed gives fresh streams, and the clock may start again at zero" begin
    path = _saved()
    straight = _run(24)
    fresh = build_ecosystem(_species(), _habitat(), from = path, seed = 99,
                            verbosity = :silent)
    simulate!(fresh, 12 * STEP, STEP, intervention = EVENTS)
    @test fresh.abundances.matrix != straight.abundances.matrix
    same = build_ecosystem(_species(), _habitat(), from = path, seed = 7,
                           verbosity = :silent)
    simulate!(same, 12 * STEP, STEP, intervention = EVENTS)
    @test same.abundances.matrix == straight.abundances.matrix
    # At zero: the population and the ground as saved, the clock afresh.
    restarted = build_ecosystem(_species(), _habitat(), from = path,
                                elapsed = 0.0s, verbosity = :silent)
    @test iszero(simulationtime(restarted))
    @test restarted.abundances.matrix == _run(12).abundances.matrix
    @test parent(restarted.habitat.regime.matrix) ≈
          parent(_run(12).habitat.regime.matrix)
end

@testset "a run may hold species the state does not, and must hold those it does" begin
    path = _saved()
    more = build_ecosystem(_species(["1", "2", "3", "4", "5"]), _habitat(),
                           from = path, verbosity = :silent)
    @test size(more.abundances.matrix, 1) == 5
    @test more.abundances.matrix[1:4, :] == _run(12).abundances.matrix
    @test iszero(sum(more.abundances.matrix[5, :]))
    @test more.spplist.traits.active == [true, true, true, true, true]
    @test_throws "holds species this run does not: `4` (key 4)" build_ecosystem(_species([
                                                                                             "1",
                                                                                             "2",
                                                                                             "3"
                                                                                         ]),
                                                                                _habitat(),
                                                                                from = path,
                                                                                verbosity = :silent)
end

@testset "a run from a state is refused on another grid, epoch or calendar, or with a placement" begin
    path = _saved()
    @test_throws "was saved on another grid" build_ecosystem(_species(),
                                                             _habitat(cells = (6.0km,
                                                                               7.0km)),
                                                             from = path,
                                                             verbosity = :silent)
    @test_throws "a run starting from a state keeps its dates" build_ecosystem(_species(),
                                                                               _habitat(),
                                                                               from = path,
                                                                               epoch = DateTime(2000),
                                                                               verbosity = :silent)
    @test_throws "a run starting from a state keeps its dates" build_ecosystem(_species(),
                                                                               _habitat(),
                                                                               from = path,
                                                                               calendar = MeanMonths(),
                                                                               verbosity = :silent)
    @test_throws "so it takes no placement" build_ecosystem(NoAbundance(),
                                                            _species(),
                                                            _habitat(),
                                                            from = path)
    @test_throws "so it goes with `from`" build_ecosystem(_species(),
                                                          _habitat(),
                                                          elapsed = 0.0s,
                                                          verbosity = :silent)
    @test_throws "ends `.jld2`" write_state(joinpath(mktempdir(), "state"),
                                            _run(1))
    @test_throws "there is no saved state" build_ecosystem(_species(),
                                                           _habitat(),
                                                           from = joinpath(mktempdir(),
                                                                           "none.jld2"))
    # A file of abundances alone - what `SaveAbundance` wrote before - says how to carry it across.
    old = joinpath(mktempdir(), "old.jld2")
    jldsave(old; abun = _run(1).abundances.matrix)
    @test_throws "copy the saved matrix into `eco.abundances.matrix`" build_ecosystem(_species(),
                                                                                      _habitat(),
                                                                                      from = old)
    @test_throws "names no file" build_ecosystem(_species(), _habitat(),
                                                 from = EcoSISTEM.ProvenanceRecord(role = :state,
                                                                                   source = "saved state"))
end

@testset "the change rules a run's interventions installed are installed again when it goes on" begin
    # Saved after the rule was installed, the run is given the same interventions and goes on as one
    # never stopped, the rule its own again from its first step, and says nothing.
    for (interventions, saved, months) in ((ONERULE, 2, 6), (TWORULES, 4, 8))
        straight = _run(months, intervention = interventions)
        path = write_state(joinpath(mktempdir(), "state.jld2"),
                           _run(saved, intervention = interventions))
        again = @test_logs build_ecosystem(_species(), _habitat(),
                                           from = path, verbosity = :silent)
        @test !isnothing(again.savedchanges)
        @test_logs simulate!(again, (months - saved) * STEP, STEP,
                             intervention = interventions)
        @test again.abundances.matrix == straight.abundances.matrix
        @test parent(again.habitat.regime.matrix) ==
              parent(straight.habitat.regime.matrix)
        @test isnothing(again.savedchanges)
    end
    # Two rules due in one step are installed in the set's order, as the run installed them, and not
    # in the order of the times they name - so the later in the set is the one in force.
    path = write_state(joinpath(mktempdir(), "state.jld2"),
                       _run(4, intervention = TWORULES))
    again = build_ecosystem(_species(), _habitat(), from = path,
                            verbosity = :silent)
    simulate!(again, 2 * STEP, 2 * STEP, intervention = TWORULES) do _
        return nothing
    end
    @test again.habitat.regime.change.value ≈ -1.0K / year
    # An intervention whose rule is in force again is among what the run was built from.
    scenario = EcoSISTEM.ProvenanceRecord(role = :intervention,
                                          source = "warming scenario")
    recorded = Intervention(AtTime(STEP), AllCells(), _rule(5.0K / year),
                            provenance = scenario)
    path = write_state(joinpath(mktempdir(), "state.jld2"),
                       _run(2, intervention = recorded))
    again = build_ecosystem(_species(), _habitat(), from = path,
                            verbosity = :silent)
    named(eco) = any(r -> r.source == "warming scenario", eco.inputs)
    @test !named(again)
    simulate!(again, STEP, STEP, intervention = recorded)
    @test named(again)
end

@testset "a run whose rules differ from those saved is told so at its first step" begin
    path = write_state(joinpath(mktempdir(), "state.jld2"),
                       _run(2, intervention = ONERULE))
    warning = (:warn,
               r"changed by a steady 5\.0 K yr.*, and this run starts with it changing by a steady 0\.5 K yr.*installed at its first step\. Where this is meant, name the layer in `build_ecosystem`'s `changing`\.$")
    bare = build_ecosystem(_species(), _habitat(), from = path,
                           verbosity = :silent)
    @test_logs warning match_mode=:any simulate!(bare, STEP, STEP)
    @test bare.habitat.regime.change.value ≈ 0.5K / year
    # Restarted at zero, a run installs its rules as it reaches them: told at the start, and given the
    # rule a month later.
    afresh = build_ecosystem(_species(), _habitat(), from = path,
                             elapsed = 0.0s, verbosity = :silent)
    @test_logs warning match_mode=:any simulate!(afresh, STEP, STEP,
                                                 intervention = ONERULE)
    @test afresh.habitat.regime.change.value ≈ 5.0K / year
end

@testset "a run told which layers change rule is quiet of them, and records every change" begin
    path = write_state(joinpath(mktempdir(), "state.jld2"),
                       _run(2, intervention = ONERULE))
    # Named, the changed rule is not warned of, and the run's provenance says what it was and is.
    told = build_ecosystem(_species(), _habitat(), from = path,
                           changing = (:Temperature,), verbosity = :silent)
    @test_logs simulate!(told, STEP, STEP)
    change = only(provenance(told).context.run.changes)
    @test change.role == "regime" && change.layer == "Temperature" &&
          change.declared
    @test occursin("5.0 K yr", change.saved) && occursin("0.5 K yr", change.now)
    @test isnothing(told.savedchanges) && told.rulechanges == [change]
    shown = sprint(show, MIME("text/plain"), provenance(told))
    @test occursin(r"changed +regime `Temperature` from a steady 5\.0 K yr.* to a steady 0\.5 K yr.*, as declared",
                   shown)
    written = TOML.parsefile(write_provenance(joinpath(mktempdir(),
                                                       "run.toml"), told))
    @test only(written["run"]["changes"])["layer"] == "Temperature" &&
          only(written["run"]["changes"])["declared"] === true
    # Not named, it is warned of and recorded all the same.
    untold = build_ecosystem(_species(), _habitat(), from = path,
                             verbosity = :silent)
    @test_logs (:warn, r"`changing`\.$") simulate!(untold, STEP, STEP)
    @test !only(provenance(untold).context.run.changes).declared
    # Named and given the saved rule again, it is warned of as unchanged, and nothing is recorded.
    same = build_ecosystem(_species(), _habitat(), from = path,
                           changing = :Temperature, verbosity = :silent)
    @test_logs (:warn,
                r"^`changing` names `Temperature`, but this run starts with its regime layer changing by a steady 5\.0 K yr.*drop `Temperature` from `changing`\.$") simulate!(same,
                                                                                                                                                                               STEP,
                                                                                                                                                                               STEP,
                                                                                                                                                                               intervention = ONERULE)
    @test isempty(provenance(same).context.run.changes)
    # A run from no state, or one whose rules are its saved ones, records no change, nor writes one.
    @test isempty(provenance(_run(1)).context.run.changes)
    plain = TOML.parsefile(write_provenance(joinpath(mktempdir(), "run.toml"),
                                            _run(1)))
    @test !haskey(plain["run"], "changes")
    # Refused: a name that is no layer of the habitat, a name that is not a symbol, and `changing`
    # without `from`.
    @test_throws "`Rainfall`, which is no layer of the habitat; its layers are `Temperature` (regime), `SolarRadiation` (supply)" build_ecosystem(_species(),
                                                                                                                                                  _habitat(),
                                                                                                                                                  from = path,
                                                                                                                                                  changing = (:Rainfall,),
                                                                                                                                                  verbosity = :silent)
    @test_throws "as symbols such as `:Temperature`" build_ecosystem(_species(),
                                                                     _habitat(),
                                                                     from = path,
                                                                     changing = ("Temperature",),
                                                                     verbosity = :silent)
    @test_throws "so it goes with `from`" build_ecosystem(BySupply(total = 10),
                                                          _species(),
                                                          _habitat(),
                                                          changing = (:Temperature,),
                                                          verbosity = :silent)
end

# A name covers every layer of that name, regime and supply alike: each is quiet where it changed,
# and warned of where it did not.
@testset "a name in `changing` covers a regime and a supply of that name" begin
    # A regime and a supply both named `Precipitation`, neither changing now; the supply saved
    # changing by a rule since replaced, or by none.
    none = EcoSISTEM.NoLayerChange()
    function saved(supplyrule)
        return EcoSISTEM._SavedChanges("state.jld2",
                                       Dict("Precipitation" =>
                                                (change = "a replaced rule",
                                                 phrase = "a SeriesLayerChange")),
                                       Dict("Precipitation" =>
                                                (change = supplyrule,
                                                 phrase = "a SeriesLayerChange")),
                                       [:Precipitation])
    end
    layers = (Precipitation = (change = none,),)
    both = saved("a replaced rule")
    compared = EcoSISTEM._comparechanges(both, layers, layers)
    @test [c.role for c in compared] == ["regime", "supply"] &&
          all(c -> c.changed && c.declared, compared)
    @test_logs EcoSISTEM._warnchanges(both, compared)
    one = saved(EcoSISTEM._changeidentity(none))
    mixed = EcoSISTEM._comparechanges(one, layers, layers)
    @test [c.changed for c in mixed] == [true, false]
    @test_logs (:warn,
                r"^`changing` names `Precipitation`, but this run starts with its supply layer") EcoSISTEM._warnchanges(one,
                                                                                                                        mixed)
end

@testset "a series a layer follows takes its values at the run's start" begin
    # A twelve-month climatology whose slice k holds k °C, so the value names the month.
    slices = cat((fill(float(i) * °C, 5, 7) for i in 1:12)..., dims = 3)
    stack = DimArray(slices,
                     (Y(NoLookup()), X(NoLookup()), Ti((1:12) .* STEP)))
    seasonal() = GridHabitat(regime = Varying(UniformSpec(1.0°C,
                                                          axis = Temperature),
                                              ReplaceWith(SeriesChange(stack,
                                                                       calendar = MonthOfYearSeries(),
                                                                       atend = RepeatAtEnd()))),
                             supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                                  axis = SolarRadiation),
                             area = StudyArea(extent = (5.0km, 7.0km),
                                              cellsize = 1.0km,
                                              verbosity = :silent))
    eco = build_ecosystem(BySupply(total = 2000), _species(), seasonal(),
                          seed = 7, verbosity = :silent)
    simulate!(eco, 5 * STEP, STEP)
    path = write_state(joinpath(mktempdir(), "state.jld2"), eco)
    continued = build_ecosystem(_species(), seasonal(), from = path,
                                verbosity = :silent)
    @test parent(continued.habitat.regime.matrix) ≈
          parent(eco.habitat.regime.matrix)
    restarted = build_ecosystem(_species(), seasonal(), from = path,
                                elapsed = 0.0s, verbosity = :silent)
    @test all(≈(uconvert(K, 1.0°C)), restarted.habitat.regime.matrix)
    @test !all(≈(uconvert(K, 1.0°C)), eco.habitat.regime.matrix)
end

@testset "a categorical layer written during a run comes back as written" begin
    # The fixture's land types are laid out at random each time its habitat is built, so a run on a
    # habitat built afresh has them only from the state.
    eco = Test1Ecosystem(seed = 3)
    classes = parent(eco.habitat.regime.matrix)
    classes[1, 1] = classes[1, 1] == 1 ? 2 : 1
    path = write_state(joinpath(mktempdir(), "state.jld2"), eco)
    fresh = Test1Ecosystem(seed = 3).habitat
    resumed = build_ecosystem(eco.spplist, fresh, from = path,
                              verbosity = :silent)
    @test parent(resumed.habitat.regime.matrix) == classes
    @test resumed.abundances.matrix == eco.abundances.matrix
end

@testset "a dated placement's rows before a state's time are its, and land on it at zero" begin
    rows = FromTable((species = ["1", "1"],
                      cell = [CartesianIndex(1, 1), CartesianIndex(1, 1)],
                      count = [5, 7], time = [3.0, 18.0] .* STEP))
    dated = Intervention(EveryStep(), ActiveCells(), AddAbundance(rows))
    eco = build_ecosystem(NoAbundance(), _species(), _habitat(), seed = 7,
                          verbosity = :silent)
    simulate!(eco, 12 * STEP, STEP, intervention = dated)
    path = write_state(joinpath(mktempdir(), "state.jld2"), eco)
    # Continuing: the row the saved run added is not added again, and nothing is said of it; the one
    # after is added when it falls due.
    continued = build_ecosystem(_species(), _habitat(), from = path,
                                verbosity = :silent)
    @test_logs simulate!(continued, 12 * STEP, STEP, intervention = dated)
    # From zero: the rows at or before it land on the state, and the run's end says so.
    restarted = build_ecosystem(_species(), _habitat(), from = path,
                                elapsed = 0.0s, verbosity = :silent)
    @test_logs((:info,
                r"^4 individuals dated at or before the run's start were added to the saved state it started from\.$"),
               match_mode = :any,
               simulate!(restarted, STEP, STEP,
                         intervention = Intervention(EveryStep(), ActiveCells(),
                                                     AddAbundance(FromTable((species = ["1"],
                                                                             cell = [CartesianIndex(1,
                                                                                                    1)],
                                                                             count = [4],
                                                                             time = [0.0s]))))))
end

@testset "SaveAbundance writes states a run can start from" begin
    folder = mktempdir()
    eco = build_ecosystem(BySupply(total = 2000), _species(), _habitat(),
                          seed = 7, verbosity = :silent)
    simulate!(SaveAbundance(folder, "run"), eco, 12 * STEP, STEP,
              every = EveryInterval(6 * STEP), intervention = EVENTS)
    @test isfile(joinpath(folder, "run02.jld2")) &&
          isfile(joinpath(folder, "run02.provenance.toml"))
    resumed = build_ecosystem(_species(), _habitat(),
                              from = joinpath(folder, "run02.jld2"),
                              verbosity = :silent)
    @test resumed.abundances.matrix == eco.abundances.matrix
    @test simulationtime(resumed) ≈ 12 * STEP
end

end
