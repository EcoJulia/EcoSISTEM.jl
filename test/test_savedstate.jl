# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/savedstate.jl`: `write_state` saving the state a run has reached, and
# `build_ecosystem(species, habitat; from)` starting an ecosystem from it - exactly as the saved run
# would have gone on, on its own seed; with fresh streams on another; with its clock at zero; with
# more species than it saved; and refused on another grid, species, epoch or calendar.
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

# An ecosystem run for `months` from the start, with the events.
function _run(months; seed = 7, kw...)
    eco = build_ecosystem(BySupply(total = 2000), _species(), _habitat(),
                          seed = seed, verbosity = :silent; kw...)
    simulate!(eco, months * STEP, STEP, intervention = EVENTS)
    return eco
end

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

@testset "a change rule an intervention installed is not saved, and the run says so" begin
    eco = build_ecosystem(BySupply(total = 2000), _species(), _habitat(),
                          seed = 7, verbosity = :silent)
    simulate!(eco, 2 * STEP, STEP,
              intervention = Intervention(AtTime(STEP), AllCells(),
                                          SetChange(nothing,
                                                    IncrementBy(5.0K / year))))
    path = write_state(joinpath(mktempdir(), "state.jld2"), eco)
    @test_logs((:warn,
                r"changed by a steady 5\.0 K yr.*, and here it changes by a steady 0\.5 K yr.*give that intervention to this run too\.$"),
               match_mode = :any,
               build_ecosystem(_species(), _habitat(), from = path,
                               verbosity = :silent))
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
