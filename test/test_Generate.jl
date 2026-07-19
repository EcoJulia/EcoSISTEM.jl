# SPDX-License-Identifier: LGPL-3.0-or-later

module TestGenerate

using EcoSISTEM
using Test
using Unitful
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using JLD2

include("TestCases.jl")

@testset "Update functions" begin
    eco = Test1Ecosystem()
    @test_nowarn update!(eco, 1month)
    @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                              1, 1, 1, eco, 10)
    @test typeof(EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                              1, 1, 1, eco, 10)) ==
          Vector{Int64}
    @test_nowarn populate!(EcoSISTEM.emptygridlandscape(eco.abenv, eco.spplist),
                           eco.spplist, eco.abenv, eco.relationship, eco.rngs)
    @test_nowarn repopulate!(eco)

    # Test Cylinder
    @test_nowarn EcoSISTEM.calc_lookup_moves!(Cylinder(), 1, 1, 1, eco, 10)
    # Test Torus
    @test_nowarn EcoSISTEM.calc_lookup_moves!(Torus(), 1, 1, 1, eco, 10)

    # Test ecosystem with multiple budgets
    eco = TestMultiEcosystem()
    @test_nowarn update!(eco, 1month)
    @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                              1, 1, 1, eco, 10)
    @test typeof(EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                              1, 1, 1, eco, 10)) ==
          Vector{Int64}
    @test_nowarn populate!(EcoSISTEM.emptygridlandscape(eco.abenv, eco.spplist),
                           eco.spplist, eco.abenv, eco.relationship, eco.rngs)
    @test_nowarn repopulate!(eco)
end

@testset "NoGrowth freezes with one or two budgets" begin
    # NoGrowth must zero the birth/death adjustment regardless of how many energy
    # budgets the environment has (the two-budget path previously skipped it).
    N = 8
    grid = (4, 4)
    area = 100.0km^2
    birth = 0.1 / month
    nogrowth = NoGrowth{typeof(unit(birth))}(fill(birth, N), fill(birth, N),
                                             1.0, 0.0, 1.0)
    traits = GaussTrait(fill(274.0K, N), fill(0.5K, N))
    movement = BirthOnlyMovement(fill(GaussianKernel(1.0km, 1.0e-3), N))
    native = fill(true, N)
    abun = fill(10, N)
    rel = Gauss{typeof(1.0K)}()

    # single budget
    sppl1 = SpeciesList(N, traits, abun, SolarRequirement(fill(10.0kJ, N)),
                        movement, nogrowth, native)
    abenv1 = simplehabitatAE(274.0K, grid, 10000.0kJ / km^2, area)
    eco1 = Ecosystem(sppl1, abenv1, rel)
    @test EcoSISTEM.energy_adjustment(eco1, eco1.abenv.budget, 1, 1) ==
          (0.0, 0.0)

    # two budgets (the previously-buggy path)
    energy2 = ReqCollection2(SolarRequirement(fill(10.0kJ, N)),
                             WaterRequirement(fill(2.0mm, N)))
    sppl2 = SpeciesList(N, traits, abun, energy2, movement, nogrowth, native)
    aenv1 = simplehabitatAE(274.0K, grid, 10000.0kJ / km^2, area)
    aenv2 = simplehabitatAE(274.0K, grid, 10.0mm / km^2, area)
    budget = BudgetCollection2(aenv1.budget, aenv2.budget)
    abenv2 = GridAbioticEnv{typeof(aenv1.habitat), typeof(budget)}(aenv1.habitat,
                                                                   aenv1.active,
                                                                   budget,
                                                                   aenv1.names)
    eco2 = Ecosystem(sppl2, abenv2, rel)
    @test EcoSISTEM.energy_adjustment(eco2, eco2.abenv.budget, 1, 1) ==
          (0.0, 0.0)
end

@testset "Multithreaded reproducibility" begin
    # A seeded run must give identical results regardless of the number of
    # threads. Run the same seeded simulation in child processes forced to use
    # 1, 2 and 4 threads (so the parallel `update!` path is genuinely exercised)
    # and check that all three final abundance matrices are identical. Running
    # the 2-thread case twice also confirms same-thread-count repeatability, so a
    # failure distinguishes "not reproducible at all" from "not reproducible
    # across thread counts".
    script = pkgdir(EcoSISTEM, "test", "threading_reproducibility.jl")
    dir = mktempdir()
    function run_with(nthreads, tag)
        out = joinpath(dir, "repro_$tag.jld2")
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -t $nthreads $script $out`
        @test success(pipeline(cmd; stdout = stdout, stderr = stderr))
        return load(out, "matrix")
    end
    m1 = run_with(1, "t1")
    m2 = run_with(2, "t2")
    m4 = run_with(4, "t4")
    m2b = run_with(2, "t2b")
    @test m1 == m2 == m4
    @test m2 == m2b
end

end
