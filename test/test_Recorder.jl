# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests for `src/Recorder.jl` - the values a run records into as it goes.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_Recorder.jl"])'

module TestRecorder

using EcoSISTEM
using EcoSISTEM: AbstractRecorder
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Diversity
using JLD2
using TOML
using Test

include("TestCases.jl")

const STEP = 1.0month_mean_duration

@testset "RecordAbundance keeps the starting state and each multiple in its slice" begin
    eco = Test1Ecosystem()
    start = copy(eco.abundances.matrix)
    recording = RecordAbundance(eco, 3)
    @test recording isa AbstractRecorder
    @test isnothing(provenance(recording))
    # Four steps; every two months is the starting state, two months and four.
    @test isnothing(simulate!(recording, eco, 4STEP, STEP,
                              every = EveryInterval(2STEP)))
    @test recording.storage[:, :, 1] == start
    @test recording.storage[:, :, 3] == eco.abundances.matrix
    @test provenance(recording).context.run.elapsed ≈
          EcoSISTEM.simulationtime(eco)
    @test startswith(sprint(show, recording), "RecordAbundance(")

    # A run reaching more occurrences than the storage holds is refused at the one it cannot keep.
    short = Test1Ecosystem()
    @test_throws "room for 2" simulate!(RecordAbundance(short, 2), short, 4STEP,
                                        STEP, every = EveryInterval(2STEP))
    # Every species of the pool has its row, an arrival's included.
    @test size(RecordAbundance(eco, 3).storage, 1) == length(eco.spplist.names)
end

@testset "a recorder's provenance is the run's as it stood at the last write" begin
    eco = Test1Ecosystem()
    recording = RecordAbundance(eco, 3)
    cull = Intervention(AtTime(STEP), AllCells(), RemoveAbundance(1, 1),
                        provenance = EcoSISTEM.ProvenanceRecord(role = :intervention,
                                                                source = "cull plan"))
    simulate!(recording, eco, STEP, STEP, intervention = cull)
    record = provenance(recording)
    @test any(r -> r.source == "cull plan", record.inputs)
    @test record.context.run.elapsed ≈ uconvert(s, STEP)
    path = write_provenance(joinpath(mktempdir(), "run.toml"), recording)
    @test haskey(TOML.parsefile(path), "run")
    @test_throws "not recorded anything" write_provenance(joinpath(mktempdir(),
                                                                   "none.toml"),
                                                          RecordAbundance(eco,
                                                                          1))
end

@testset "SaveAbundance writes each occurrence and its provenance" begin
    dir = joinpath(mktempdir(), "out")
    eco = Test1Ecosystem()
    saving = SaveAbundance(dir, "run")
    simulate!(saving, eco, STEP, STEP)
    @test isfile(joinpath(dir, "run00.jld2")) &&
          isfile(joinpath(dir, "run01.jld2"))
    @test !isfile(joinpath(dir, "run02.jld2"))
    # Each a state: its abundances, as a later run starting from it takes them.
    @test JLD2.load(joinpath(dir, "run01.jld2"), "abundances") ==
          eco.abundances.matrix
    @test JLD2.load(joinpath(dir, "run01.jld2"), "format") == "EcoSISTEM state"
    written = TOML.parsefile(joinpath(dir, "run01.provenance.toml"))
    @test written["software"]["package"] == "EcoSISTEM"
    # Named apart from an input's own record, so it is never read back as one.
    @test !isfile(joinpath(dir, "run01.jld2.provenance.toml"))
    @test sprint(show, saving) == "SaveAbundance($(repr(dir)), \"run\")"
end

@testset "RecordDiversity keeps a measure at each order" begin
    eco = Test1Ecosystem()
    qs = [0.0, 1.0, 2.0]
    ncells = size(eco.abundances.matrix, 2)
    recording = RecordDiversity(generate_storage(eco, length(qs), 2, 1),
                                norm_sub_alpha, qs)
    simulate!(recording, eco, STEP, STEP)
    @test recording.storage[:, :, 2] ≈
          reshape(norm_sub_alpha(eco, qs)[!, :diversity], ncells, length(qs))

    # Several things kept from one run: each recorder is called from a callback of your own.
    both = Test1Ecosystem()
    abundance = RecordAbundance(both, 2)
    diversity = RecordDiversity(generate_storage(both, 1, 2, 1), norm_sub_alpha,
                                [1.0])
    simulate!(both, STEP, STEP) do occurrence
        abundance(both, occurrence)
        return diversity(both, occurrence)
    end
    @test abundance.storage[:, :, 2] == both.abundances.matrix
    @test provenance(diversity).context.run == provenance(abundance).context.run

    # A metacommunity measure has no value per cell to keep, and is refused whatever the storage's
    # shape - including a single row, which its values would otherwise fit.
    occurrence = (count = 1, elapsed = 0.0s, date = nothing)
    @test_throws "subcommunity diversity measure" simulate!(RecordDiversity(generate_storage(eco,
                                                                                             length(qs),
                                                                                             3,
                                                                                             1),
                                                                            meta_gamma,
                                                                            qs),
                                                            Test1Ecosystem(),
                                                            STEP, STEP)
    @test_throws "gives metacommunity diversity" RecordDiversity(zeros(1,
                                                                       length(qs),
                                                                       1),
                                                                 norm_meta_alpha,
                                                                 qs)(Test1Ecosystem(),
                                                                     occurrence)
end

@testset "RecordTotals keeps each species' and each cell's individuals" begin
    eco = Test1Ecosystem()
    start = copy(eco.abundances.matrix)
    totals = RecordTotals(eco, 3)
    @test totals isa AbstractRecorder
    @test size(totals.species) == (length(eco.spplist.names), 3)
    @test size(totals.cells) == (size(start, 2), 3)
    simulate!(totals, eco, 4STEP, STEP, every = EveryInterval(2STEP))
    @test totals.species[:, 1] == vec(sum(start, dims = 2))
    @test totals.cells[:, 1] == vec(sum(start, dims = 1))
    @test totals.species[:, 3] == vec(sum(eco.abundances.matrix, dims = 2))
    @test totals.cells[:, 3] == vec(sum(eco.abundances.matrix, dims = 1))
    @test sum(totals.species[:, 3]) == sum(totals.cells[:, 3])
    @test provenance(totals).context.run.elapsed ≈
          EcoSISTEM.simulationtime(eco)
    @test sprint(show, totals) ==
          "RecordTotals(15 species, 35 cells, 3 occurrences)"
    short = Test1Ecosystem()
    @test_throws "room for 2" simulate!(RecordTotals(short, 2), short, 4STEP,
                                        STEP, every = EveryInterval(2STEP))
end

# Diversity's measure types and the functions computing each at each level, which a recording must
# match.
const SUBMEASURES = (RawAlpha => raw_sub_alpha,
                     NormalisedAlpha => norm_sub_alpha,
                     RawBeta => raw_sub_beta, NormalisedBeta => norm_sub_beta,
                     RawRho => raw_sub_rho, NormalisedRho => norm_sub_rho,
                     Diversity.Gamma => sub_gamma)
const METAMEASURES = (RawAlpha => raw_meta_alpha,
                      NormalisedAlpha => norm_meta_alpha,
                      RawBeta => raw_meta_beta,
                      NormalisedBeta => norm_meta_beta,
                      RawRho => raw_meta_rho, NormalisedRho => norm_meta_rho,
                      Diversity.Gamma => meta_gamma)

# Species with no similarity between them, on a non-square grid with a supply that keeps them all,
# placed at random so that abundances vary by species and by cell.
function uniqueecosystem(numspecies, grid)
    species = build_species(SpeciesPool(numspecies, verbosity = :silent),
                            tolerance = ToleranceSpec(Temperature, μ = 298.0K,
                                                      σ = 2.0K),
                            demand = DemandSpec(SolarRadiation, 2.0kJ / day),
                            verbosity = :silent)
    area = StudyArea(extent = (grid[1] * 1.0km, grid[2] * 1.0km),
                     cellsize = 1.0km, verbosity = :silent)
    habitat = build_habitat(area = area,
                            supply = UniformSpec(5.0e3kJ / km^2 / day,
                                                 axis = SolarRadiation),
                            verbosity = :silent)
    return build_ecosystem(BySupply(total = 1000 * numspecies), species,
                           habitat, seed = 1, verbosity = :silent)
end

@testset "RecordDiversity keeps Diversity's measures at both levels" begin
    qs = [0.0, 1.0, 2.0, Inf]
    # Without and with a phylogeny, every value is the one Diversity's own function gives.
    for eco in (uniqueecosystem(30, (6, 8)), Test1Ecosystem())
        ncells = size(eco.abundances.matrix, 2)
        recording = RecordDiversity(eco, 3, qs = qs,
                                    subcommunity = first.(SUBMEASURES),
                                    metacommunity = first.(METAMEASURES))
        @test size(recording.subcommunity) ==
              (length(SUBMEASURES), ncells, length(qs), 3)
        @test size(recording.metacommunity) == (length(METAMEASURES),
               length(qs), 3)
        @test all(isnan, recording.metacommunity)
        simulate!(recording, eco, 2STEP, STEP)
        for (i, (_, f)) in enumerate(SUBMEASURES)
            @test recording.subcommunity[i, :, :, 3] ==
                  reshape(f(eco, qs)[!, :diversity], ncells, length(qs))
        end
        for (i, (_, f)) in enumerate(METAMEASURES)
            @test recording.metacommunity[i, :, 3] == f(eco, qs)[!, :diversity]
        end
        @test !any(isnan, recording.metacommunity)
    end

    # A measure listed at one level only is kept at that level only, in the order listed, and a
    # level listing none keeps nothing.
    eco = uniqueecosystem(30, (6, 8))
    listed = RecordDiversity(eco, 1, qs = [1.0],
                             metacommunity = [Diversity.Gamma, RawBeta])
    @test isnothing(listed.subcommunity)
    listed(eco, (count = 1, elapsed = 0.0s, date = nothing))
    @test listed.metacommunity[:, 1, 1] ==
          [
        meta_gamma(eco, 1.0)[1, :diversity],
        raw_meta_beta(eco, 1.0)[1, :diversity]
    ]
    @test sprint(show, listed) ==
          "RecordDiversity(subcommunity = [], metacommunity = [Gamma, RawBeta], qs = [1.0])"
    cells = RecordDiversity(eco, 1, qs = [0.0],
                            subcommunity = [NormalisedAlpha])
    @test isnothing(cells.metacommunity) && isnothing(cells.storage)
    @test cells.measures.subcommunity == [NormalisedAlpha]

    # Refused: no measure, an empty order list, a function or anything else that is not a concrete
    # measure type, a measure listed twice, and an occurrence with no room.
    @test_throws "needs a measure to record" RecordDiversity(eco, 1, qs = qs)
    @test_throws "at least one order" RecordDiversity(eco, 1, qs = Float64[],
                                                      subcommunity = [RawAlpha])
    @test_throws "`RecordDiversity(storage, measure, qs)`" RecordDiversity(eco,
                                                                           1,
                                                                           qs = qs,
                                                                           subcommunity = [norm_sub_alpha])
    @test_throws "was given Int64" RecordDiversity(eco, 1, qs = qs,
                                                   metacommunity = [Int])
    @test_throws "DiversityMeasure" RecordDiversity(eco, 1, qs = qs,
                                                    metacommunity = [Diversity.DiversityMeasure])
    @test_throws "RawRho` twice" RecordDiversity(eco, 1, qs = qs,
                                                 metacommunity = [RawRho,
                                                     RawRho])
    @test_throws "room for 1" simulate!(RecordDiversity(eco, 1, qs = qs,
                                                        metacommunity = [RawRho]),
                                        eco, STEP, STEP)
end

# One snapshot per recording: what a recording allocates is what Diversity's measures need for
# their per-cell values, never a copy of the abundances. Building the snapshot from
# `abun / sum(abun)` makes this fail; a cost that grows with the cells, measures and orders but stays
# below one copy of the abundances is not caught.
@testset "a diversity recording allocates less than one copy of the abundances" begin
    eco = uniqueecosystem(200, (40, 50))
    recording = RecordDiversity(eco, 1, qs = [0.0, 1.0, 2.0, Inf],
                                subcommunity = [NormalisedAlpha],
                                metacommunity = first.(METAMEASURES))
    occurrence = (count = 1, elapsed = 0.0s, date = nothing)
    function allocated(recorder, eco, occurrence)
        recorder(eco, occurrence)
        return @allocated recorder(eco, occurrence)
    end
    @test allocated(recording, eco, occurrence) < sizeof(eco.abundances.matrix)
end

end
