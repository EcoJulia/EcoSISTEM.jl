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
    @test_nowarn populate!(EcoSISTEM.emptygridlandscape(eco.habitat,
                                                        eco.spplist),
                           eco.spplist, eco.habitat, eco.nichefit, eco.rngs)
    @test_nowarn repopulate!(eco)

    # Test Cylinder
    @test_nowarn EcoSISTEM.calc_lookup_moves!(Cylinder(), 1, 1, 1, eco, 10)
    # Test Torus
    @test_nowarn EcoSISTEM.calc_lookup_moves!(Torus(), 1, 1, 1, eco, 10)

    # Test ecosystem with multiple supplies
    eco = TestMultiEcosystem()
    @test_nowarn update!(eco, 1month)
    @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                              1, 1, 1, eco, 10)
    @test typeof(EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                              1, 1, 1, eco, 10)) ==
          Vector{Int64}
    @test_nowarn populate!(EcoSISTEM.emptygridlandscape(eco.habitat,
                                                        eco.spplist),
                           eco.spplist, eco.habitat, eco.nichefit, eco.rngs)
    @test_nowarn repopulate!(eco)
end

@testset "NoGrowth freezes with one or two supplies" begin
    # NoGrowth must zero the birth/death adjustment regardless of how many resource
    # supplies the environment has (the two-supply path previously skipped it).
    N = 8
    grid = (4, 4)
    area = 100.0km^2
    birth = 0.1 / month
    nogrowth = NoGrowth{typeof(unit(birth))}(fill(birth, N), fill(birth, N),
                                             1.0, 0.0, 1.0)
    tolerance = NicheTolerance(Temperature, Normal, fill(274.0K, N),
                               fill(0.5K, N))
    movement = BirthOnlyMovement(fill(GaussianKernel(1.0km, 1.0e-3), N))
    native = fill(true, N)
    abun = fill(10, N)
    nichefit = NicheSuitability{typeof(1.0K)}()

    # single supply
    sppl1 = SpeciesList(N, tolerance, abun, SolarDemand(fill(10.0kJ, N)),
                        movement, nogrowth, native)
    habitat1 = simplehabitat(274.0K, grid, 10000.0kJ / km^2, area)
    eco1 = Ecosystem(sppl1, habitat1, nichefit)
    @test EcoSISTEM.resource_adjustment(eco1, eco1.habitat.supply, 1, 1) ==
          (0.0, 0.0)

    # two supplies (the previously-buggy path)
    resource2 = DemandCollection2(SolarDemand(fill(10.0kJ, N)),
                                  WaterDemand(fill(2.0mm, N)))
    sppl2 = SpeciesList(N, tolerance, abun, resource2, movement, nogrowth,
                        native)
    habitat_solar = simplehabitat(274.0K, grid, 10000.0kJ / km^2, area)
    habitat_water = simplehabitat(274.0K, grid, 10.0mm / km^2, area)
    supply = SupplyCollection2(habitat_solar.supply, habitat_water.supply)
    habitat2 = GridHabitat{typeof(habitat_solar.regime), typeof(supply)}(habitat_solar.regime,
                                                                         habitat_solar.active,
                                                                         supply,
                                                                         habitat_solar.names)
    eco2 = Ecosystem(sppl2, habitat2, nichefit)
    @test EcoSISTEM.resource_adjustment(eco2, eco2.habitat.supply, 1, 1) ==
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
