using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units

include("TestCases.jl")

epi = TestEpiSystem()

@test sum(epi.abundances.matrix, dims = 2)[:, 1] == epi.epilist.human.abun
@test_nowarn gettraitrel(epi)
@test gettraitrel(epi) == epi.relationship
@test_nowarn gethabitat(epi)
@test gethabitat(epi) == epi.epienv.habitat
@test_nowarn getsize(epi)
@test getsize(epi) == size(epi.abundances.matrix, 2) .* epi.epienv.habitat.size^2
@test_nowarn getgridsize(epi)
@test getgridsize(epi) == epi.epienv.habitat.size
@test_nowarn getdispersaldist(epi, 1)
@test getdispersaldist(epi, 1) == epi.epilist.human.movement.kernels[1].dist
@test_nowarn getdispersaldist(epi, "Susceptible")
@test getdispersaldist(epi, "Susceptible") == epi.epilist.human.movement.kernels[1].dist
@test_nowarn getdispersalvar(epi, 1)
@test_nowarn getdispersalvar(epi, "Susceptible")

@testset "EpiSystem from initial population" begin
    initial_pop = rand(10, 10) * 10
    epi = TestEpiSystemFromPopulation(initial_pop)
    @test epi.epienv.initial_population == round.(initial_pop)
    @test epi.abundances.grid[1, :, :] == round.(initial_pop)
end

@testset "Approximate equality" begin
    atol = 1
    epi_2 = deepcopy(epi)
    @test isapprox(epi_2, epi; atol=atol)
    # Make some small change within the approximate equality range
    epi_2.abundances.matrix[1, 1] += 1
    @test isapprox(epi_2, epi; atol=atol)
    # Make some larger change outside the approximate equality range
    epi_2.abundances.matrix[1, 1] += 1000
    @test !isapprox(epi_2, epi; atol=atol)
end

@testset "save and load (with JLSO)" begin
    # Test saving Function
    Simulation.save("testepi.jlso", epi)
    loaded_epi = Simulation.load("testepi.jlso", EpiSystem)
    # Ideally we'd compare epi and loaded_epi, but `EpiSystem` still does not support comparisons
    @test loaded_epi isa EpiSystem
    rm("testepi.jlso") # Delete temporary file
end
