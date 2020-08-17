using Simulation
using Test
using Unitful.DefaultSymbols
using Distributions
using Simulation.Units

include("TestCases.jl")
@testset "Update functions" begin
    eco = TestEcosystem()
    @test_nowarn update!(eco, 1month)
    @test_nowarn Simulation.calc_lookup_moves!(eco.spplist.movement.boundary, 1, 1, 1, eco, 10)
    @test typeof(Simulation.calc_lookup_moves!(eco.spplist.movement.boundary, 1, 1, 1, eco,
        10)) == Vector{Int64}
    @test_nowarn populate!(Simulation.emptygridlandscape(eco.abenv, eco.spplist), eco.spplist, eco.abenv, eco.relationship)
    @test_nowarn repopulate!(eco)

    # Test Cylinder
    @test_nowarn Simulation.calc_lookup_moves!(Cylinder(), 1, 1, 1, eco, 10)
    # Test Torus
    @test_nowarn Simulation.calc_lookup_moves!(Torus(), 1, 1, 1, eco, 10)

    # Test ecosystem with multiple budgets
    eco = TestMultiEcosystem()
    @test_nowarn update!(eco, 1month)
    @test_nowarn Simulation.calc_lookup_moves!(eco.spplist.movement.boundary, 1, 1, 1, eco, 10)
    @test typeof(Simulation.calc_lookup_moves!(eco.spplist.movement.boundary, 1, 1, 1, eco,
        10)) == Vector{Int64}
    @test_nowarn populate!(Simulation.emptygridlandscape(eco.abenv, eco.spplist), eco.spplist, eco.abenv, eco.relationship)
    @test_nowarn repopulate!(eco)
end
