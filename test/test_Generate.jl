# SPDX-License-Identifier: LGPL-3.0-or-later

module TestGenerate

using EcoSISTEM
using Test
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units

include("TestCases.jl")

@testset "Update functions" begin
    eco = Test1Ecosystem()
    @test_nowarn update!(eco, 1month)
    @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary, 1,
                                              1, 1, eco, 10)
    @test typeof(EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary, 1,
                                              1, 1, eco,
                                              10)) == Vector{Int64}
    @test_nowarn populate!(EcoSISTEM.emptygridlandscape(eco.abenv, eco.spplist),
                           eco.spplist, eco.abenv, eco.relationship)
    @test_nowarn repopulate!(eco)

    # Test Cylinder
    @test_nowarn EcoSISTEM.calc_lookup_moves!(Cylinder(), 1, 1, 1, eco, 10)
    # Test Torus
    @test_nowarn EcoSISTEM.calc_lookup_moves!(Torus(), 1, 1, 1, eco, 10)

    # Test ecosystem with multiple budgets
    eco = TestMultiEcosystem()
    @test_nowarn update!(eco, 1month)
    @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary, 1,
                                              1, 1, eco, 10)
    @test typeof(EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary, 1,
                                              1, 1, eco,
                                              10)) == Vector{Int64}
    @test_nowarn populate!(EcoSISTEM.emptygridlandscape(eco.abenv, eco.spplist),
                           eco.spplist, eco.abenv, eco.relationship)
    @test_nowarn repopulate!(eco)
end

@testset "Multithreaded reproducibility" begin
    # Run the reproducibility check in a child process forced to use several
    # threads, so the parallel `update!` path is exercised even when the parent
    # test suite is launched single-threaded.
    script = EcoSISTEM.path("threading_reproducibility.jl")
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -t 4 $script`
    @test success(pipeline(cmd; stdout = stdout, stderr = stderr))
end

end
