module TestEpiSystem

using EcoSISTEM
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units

include("TestCases.jl")

@testset "Episystem" begin
    epi = Test1EpiSystem()

    @test sum(epi.abundances.matrix, dims = 2)[:, 1] == epi.spplist.species.abun
    @test_nowarn gettraitrel(epi)
    @test gettraitrel(epi) == epi.relationship
    @test_nowarn gethabitat(epi)
    @test gethabitat(epi) == epi.abenv.habitat
    @test_nowarn getsize(epi)
    @test getsize(epi) ==
          size(epi.abundances.matrix, 2) .* epi.abenv.habitat.size^2
    @test_nowarn getgridsize(epi)
    @test getgridsize(epi) == epi.abenv.habitat.size
    @test_nowarn getdispersaldist(epi, 1)
    @test getdispersaldist(epi, 1) ==
          epi.spplist.species.movement.localmoves.kernels[1].dist
    @test_nowarn getdispersaldist(epi, "Susceptible")
    @test getdispersaldist(epi, "Susceptible") ==
          epi.spplist.species.movement.localmoves.kernels[1].dist
    @test_nowarn getdispersalvar(epi, 1)
    @test_nowarn getdispersalvar(epi, "Susceptible")

    @testset "Approximate equality" begin
        atol = 1
        epi_2 = deepcopy(epi)
        @test isapprox(epi_2, epi; atol = atol)
        # Make some small change within the approximate equality range
        epi_2.abundances.matrix[1, 1] += 1
        @test isapprox(epi_2, epi; atol = atol)
        # Make some larger change outside the approximate equality range
        epi_2.abundances.matrix[1, 1] += 1000
        @test !isapprox(epi_2, epi; atol = atol)
    end

    @testset "save and load (with JLSO)" begin
        # Test saving Function
        EcoSISTEM.save("testepi.jlso", epi)
        loaded_epi = EcoSISTEM.load("testepi.jlso", Ecosystem)
        # Ideally we'd compare epi and loaded_epi, but `EpiSystem` still does not support comparisons
        @test loaded_epi isa Ecosystem
        rm("testepi.jlso") # Delete temporary file
    end
end

end
