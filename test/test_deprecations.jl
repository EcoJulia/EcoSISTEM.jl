# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDeprecations

using EcoSISTEM
# `[C7-VIS]` C: these are `public` rather than exported - a spec is what a user writes,
# and these are what it materialises into.
using EcoSISTEM: RateChange, SeriesLayerChange, AbsoluteChange
# `[C7-VIS]` B1/B2/B3: these are `public` rather than exported, so they must be named.
using EcoSISTEM: getdist
using Extents: Extent
using Test
using Distributions
using Unitful, Unitful.DefaultSymbols
using EcoSISTEM.Units
using RasterDataSources
using DimensionalData: DimensionalData, DimArray, Y, X, Ti, Dim
using Diversity: norm_sub_alpha, norm_sub_beta

include("TestCases.jl")

# Coverage for `src/deprecations.jl` and the dataset-typed shims in the RasterDataSources extension.
# Every shim is checked on *both* halves: it warns (`@test_deprecated`) **and** its result matches
# the current API it forwards to.

@testset "Deprecations" begin
    @testset "simulate_action! keeps its own timing" begin
        # One step after each multiple of the interval, over the steps `simulate!` takes, handed a
        # bare count; with `offset` the grid starts at the timestep and the run is a step shorter.
        step = 1.0month_mean_duration
        fired(offset) = begin
            eco = Test1Ecosystem()
            calls = Tuple{Int, typeof(1.0s)}[]
            @test_deprecated simulate_action!(eco, 3step, step, step,
                                              offset = offset) do counting
                return push!(calls,
                             (counting, EcoSISTEM.simulationtime(eco)))
            end
            calls
        end
        plain = fired(false)
        @test first.(plain) == 1:4
        @test last.(plain) ≈ [uconvert(s, k * step) for k in 1:4]
        shifted = fired(true)
        @test first.(shifted) == 1:3
        @test last.(shifted) ≈ [uconvert(s, k * step) for k in 1:3]
    end

    @testset "the recording functions keep their own timing" begin
        step = 1.0month_mean_duration
        # The starting state, then each multiple up to `times`, over `times / timestep` steps.
        eco = Test1Ecosystem()
        start = copy(eco.abundances.matrix)
        storage = generate_storage(eco, 4, 1)
        @test_deprecated simulate_record!(storage, eco, 3step, step, step)
        @test storage[:, :, 1] == start
        @test storage[:, :, 4] == eco.abundances.matrix
        @test EcoSISTEM.simulationtime(eco) ≈ uconvert(s, 3step)

        # The caching `simulate!` saves on the step after each multiple, numbered from `00`.
        dir = mktempdir()
        @test_deprecated simulate!(Test1Ecosystem(), 3step, step, step, dir,
                                   "testrun")
        @test isfile(joinpath(dir, "testrun00.jld2")) &&
              isfile(joinpath(dir, "testrun03.jld2"))

        # The three diversity forms, on `simulate_action!`'s timing.
        qs = collect(1.0:3)
        eco = Test1Ecosystem()
        ncells = size(eco.abundances.matrix, 2)
        @test_deprecated simulate_record_diversity!(generate_storage(eco,
                                                                     length(qs),
                                                                     4, 1),
                                                    eco, 3step, step, step,
                                                    norm_sub_alpha, qs)
        divfuns = Function[norm_sub_alpha, norm_sub_beta]
        @test_deprecated simulate_record_diversity!(generate_storage(Test1Ecosystem(),
                                                                     length(divfuns),
                                                                     4, 1),
                                                    Test1Ecosystem(), 3step,
                                                    step, step, divfuns, 1.0)
        sub = zeros(Float64, ncells, 3, 3, 4)
        meta = zeros(Float64, 3, 3, 4)
        result = @test_deprecated simulate_record_diversity!(sub, meta,
                                                             Test1Ecosystem(),
                                                             3step, step, step,
                                                             qs)
        @test result.subcommunity === sub && result.metacommunity === meta
    end

    # The reads need downloaded raster data, which is unavailable or slow on Windows CI.
    if !Sys.iswindows()
        # The dataset-typed `read` methods and `readfile` are spellings of a `RasterSpec` read now,
        # each kept for one release. `read(T, layers)` keeps its old contract - bare magnitudes -
        # so it is compared to the spec form with the unit stripped.
        @testset "dataset-typed reads -> read(::RasterSpec)" begin
            wind = getraster(WorldClim{Climate}, :wind, month = 1:12)
            bio1 = getraster(WorldClim{BioClim}, :bio1)
            @test_deprecated read(WorldClim{Climate}, :wind, month = 1:2)
            @test isequal(read(WorldClim{Climate}, :wind, month = 1:2).array,
                          ustrip.(read(SourceSpec(WorldClim{Climate}, :wind,
                                                  month = 1:2)).array))
            winddir = dirname(first(wind))
            @test_deprecated read(CRUTS, winddir, "tavg")
            @test isequal(read(CRUTS, winddir, "tavg").array,
                          read(SourceSpec(CRUTS, "tavg", directory = winddir)).array)
            # The CHELSA directory reader keeps its own body: a directory of another provider's
            # files is not checked against CHELSA's catalogue row, as a spec would be.
            @test_deprecated read(CHELSA{Climate}, winddir, "wind")
            @test unit(read(CHELSA{Climate}, winddir, "wind").array[1]) == m / s
            @test_deprecated readfile(bio1)
            @test isequal(readfile(bio1).array,
                          read(RasterFileSpec(bio1,
                                              axis = EcoSISTEM.NicheAxis)).array)
            @test_deprecated readfile(bio1, unit = K)
            @test unit(readfile(bio1, unit = K).array[1]) == K
        end
    end
end

@testset "demographics: boost is gone, the birth cap is 1" begin
    birth = fill(0.6 / year, 3)
    death = fill(0.6 / year, 3)
    U = typeof(unit(first(birth)))
    # Each five-argument form warns and gives the four-field value.
    @test_deprecated EqualPop(0.6 / year, 0.6 / year, 1.0, 0.2, 1.0)
    @test_deprecated PopGrowth{U}(birth, death, 1.0, 0.2, 100.0)
    @test_deprecated NoGrowth{U}(birth, death, 1.0, 0.2, 1.0)
    @test EqualPop(0.6 / year, 0.6 / year, 1.0, 0.2, 1.0) ==
          EqualPop(0.6 / year, 0.6 / year, 1.0, 0.2)
    @test fieldnames(EqualPop) == (:birth, :death, :longevity, :survival)
    @test fieldnames(PopGrowth) == fieldnames(EqualPop)
    @test fieldnames(NoGrowth) == fieldnames(EqualPop)
    # The keyword warns too, and the default does not.
    @test_deprecated build_species(DefaultEcosystem(), boost = 10.0,
                                   verbosity = :silent)
    @test_nowarn build_species(DefaultEcosystem(), verbosity = :silent)
end

end
