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
    @test isnothing(simulate!(recording, eco, 3STEP, STEP,
                              every = EveryInterval(2STEP)))
    @test recording.storage[:, :, 1] == start
    @test recording.storage[:, :, 3] == eco.abundances.matrix
    @test provenance(recording).run.elapsed ≈ EcoSISTEM.simulationtime(eco)
    @test startswith(sprint(show, recording), "RecordAbundance(")

    # A run reaching more occurrences than the storage holds is refused at the one it cannot keep.
    short = Test1Ecosystem()
    @test_throws "room for 2" simulate!(RecordAbundance(short, 2), short, 3STEP,
                                        STEP, every = EveryInterval(2STEP))
    # ...and storage with no room for the species already there is refused before the run.
    @test_throws "maxspecies" RecordAbundance(eco, 3, maxspecies = 1)
end

@testset "a recorder's provenance is the run's as it stood at the last write" begin
    eco = Test1Ecosystem()
    recording = RecordAbundance(eco, 3)
    cull = Intervention(AtTime(STEP), AllCells(), RemoveAbundance(1, 1),
                        provenance = EcoSISTEM.InputRecord(role = :intervention,
                                                           dataset = "cull plan"))
    simulate!(recording, eco, STEP, STEP, intervention = cull)
    record = provenance(recording)
    @test any(r -> r.dataset == "cull plan", record.inputs)
    @test record.run.elapsed ≈ uconvert(s, 2STEP)
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
          isfile(joinpath(dir, "run02.jld2"))
    @test !isfile(joinpath(dir, "run03.jld2"))
    @test JLD2.load(joinpath(dir, "run02.jld2"), "abun") ==
          eco.abundances.matrix
    written = TOML.parsefile(joinpath(dir, "run02.provenance.toml"))
    @test written["software"]["package"] == "EcoSISTEM"
    # Named apart from an input's own record, so it is never read back as one.
    @test !isfile(joinpath(dir, "run02.jld2.provenance.toml"))
    @test sprint(show, saving) == "SaveAbundance($(repr(dir)), \"run\")"
end

@testset "RecordDiversity keeps a measure at each order" begin
    eco = Test1Ecosystem()
    qs = [0.0, 1.0, 2.0]
    ncells = size(eco.abundances.matrix, 2)
    recording = RecordDiversity(generate_storage(eco, length(qs), 3, 1),
                                norm_sub_alpha, qs)
    simulate!(recording, eco, STEP, STEP)
    @test recording.storage[:, :, 3] ≈
          reshape(norm_sub_alpha(eco, qs)[!, :diversity], ncells, length(qs))

    # Several things kept from one run: each recorder is called from a callback of your own.
    both = Test1Ecosystem()
    abundance = RecordAbundance(both, 3)
    diversity = RecordDiversity(generate_storage(both, 1, 3, 1), norm_sub_alpha,
                                [1.0])
    simulate!(both, STEP, STEP) do occurrence
        abundance(both, occurrence)
        return diversity(both, occurrence)
    end
    @test abundance.storage[:, :, 3] == both.abundances.matrix
    @test provenance(diversity).run == provenance(abundance).run

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

end
