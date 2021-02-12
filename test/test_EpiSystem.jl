using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units

include("TestCases.jl")

@testset "Episystem" begin
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
    @test getdispersaldist(epi, 1) == epi.epilist.human.movement.home.kernels[1].dist
    @test_nowarn getdispersaldist(epi, "Susceptible")
    @test getdispersaldist(epi, "Susceptible") == epi.epilist.human.movement.home.kernels[1].dist
    @test_nowarn getdispersalvar(epi, 1)
    @test_nowarn getdispersalvar(epi, "Susceptible")

    @testset "EpiSystem from initial population" begin
        @testset "All active" begin
            initial_pop = rand(10, 10) * 10
            epi = TestEpiSystemFromPopulation(initial_pop)
            @test epi.abundances.grid[1, :, :] == round.(initial_pop)
        end

        @testset "Some inactive" begin
            initial_pop = Matrix{Union{Float64, Missing}}(rand(10, 15) * 10)
            # Fill some NaNs and missing, these should indicate inactive areas
            missing_idx = CartesianIndex.([1, 4, 8], [5, 8, 2])
            initial_pop[missing_idx[1]] = NaN
            initial_pop[missing_idx[2]] = missing
            initial_pop[missing_idx[3]] = missing
            initial_pop_ref = copy(initial_pop)

            expected_active = fill(true, size(initial_pop))
            expected_active[missing_idx] .= false
            expected_size = size(initial_pop)
            expected_pop = round.(replace(initial_pop, NaN => 0, missing => 0))

            epi = TestEpiSystemFromPopulation(initial_pop)
            @test size(epi.epienv.habitat.matrix) == expected_size
            @test epi.epienv.active == expected_active
            @test epi.abundances.grid[1, :, :] == expected_pop
            # Should not have modified initial_pop
            @test isequal(initial_pop, initial_pop_ref)

            # If we specify inactive regions when constructing the EpiEnv, these should be
            # preserved
            epienv_active = fill(true, size(initial_pop))
            epienv_active[1, 1] = false
            epienv_active[2, 2] = false
            epienv_active[3, 3] = false
            expected_active = copy(epienv_active)
            expected_active[missing_idx] .= false
            @test expected_active != epienv_active

            epi = TestEpiSystemFromPopulation(initial_pop; epienv_active=epienv_active)
            @test size(epi.epienv.habitat.matrix) == expected_size
            @test epi.epienv.active == expected_active
            @test epi.abundances.grid[1, :, :] == expected_pop
            # Should not have modified initial_pop
            @test isequal(initial_pop, initial_pop_ref)
        end

        @testset "Initial population is provided as a AxisArray Matrix" begin
            grid = (10, 10)
            active = fill(true, grid)
            # Set inactive cells on the outer layers
            # These should not be trimmed off by the EpiSystem constructor
            for row in (1, 2, 3, 10)
                active[row, :] .= false
            end
            for col in (1, 2, 10)
                active[:, col] .= false
            end
            initial_pop = AxisArray(
                zeros(grid...);
                x=Symbol.("x_", 1:grid[1]),
                y=Symbol.("y_", 1:grid[2]),
            )
            initial_pop.data[.!active] .= NaN

            expected_active = active
            expected_initial_pop = Int.(zeros(grid))
            expected_initial_pop[.!active] .= 0

            epi = TestEpiSystemFromPopulation(initial_pop)

            @test epi.epienv.active == expected_active
            @test size(epi.epienv.habitat.matrix) == grid
            @test epi.abundances.grid[1, :, :] == expected_initial_pop
        end
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
end
