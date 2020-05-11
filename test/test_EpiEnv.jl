using Simulation
using Unitful.DefaultSymbols
using Test
using Simulation.Units
using Simulation.ClimatePref
using AxisArrays

grid = (5, 5)
area = 25.0km^2
active = fill(true, grid)
control = NoControl()

# TEST simplehabitatAE
fillval = 1.0
abenv = simplehabitatAE(fillval, grid, area, control)
@test_nowarn simplehabitatAE(fillval, grid, area, control)
@test_nowarn simplehabitatAE(fillval, grid, area, active, control)
@test all(abenv.habitat.matrix .== fillval)
@test size(abenv.habitat.matrix) == grid
@test abenv.active == active
@test all(abenv.active)
@test abenv.habitat.size^2 * grid[1] * grid[2] == area

@testset "Constructing from initial population" begin
    fillval = 1.0
    area = 25.0km^2
    control = NoControl()
    initial_pop = rand(10, 15) * 10
    # Fill some NaNs, these should indicate inactive areas
    missing_idx = CartesianIndex.([1, 4, 8], [5, 8, 2])
    initial_pop[missing_idx] .= NaN

    expected_active = Matrix{Bool}(.!isnan.(initial_pop))
    expected_size = size(initial_pop)
    expected_pop = round.(replace(initial_pop, NaN => 0))

    epienv = simplehabitatAE(fillval, area, control, initial_pop)
    @test all(epienv.habitat.matrix .== fillval)
    @test size(epienv.habitat.matrix) == expected_size
    @test epienv.active == expected_active
    @test epienv.initial_population == expected_pop
end

@testset "Shrink grid" begin
    grid = (10, 10)
    area = 100.0km^2
    fillval = 1.0

    active = fill(true, grid)
    # Set inactive cells on the outer two layers
    # These should be trimmed off
    for row in (1, 2, 9, 10)
        active[row, :] .= false
    end
    for col in (1, 2, 9, 10)
        active[:, col] .= false
    end
    # Set some more inactive cells, these shouldn't result in any further trimming
    active[3, 3] = false
    active[4, 3] = false
    active[5, 5] = false
    active[8, 8] = false
    # Should be trimmed from 10^2 to 8^2
    expected_grid= (8, 8)
    expected_active = active[3:8, 3:8]
    expected_area = area * (8/10)^2
    # Shouldn't change
    expected_gridlength = 1km
    control = NoControl()
    expected_matrix = fill(fillval, expected_grid)

    @testset "Construct directly" begin
        epienv = simplehabitatAE(fillval, grid, area, active, control)

        @test epienv.active == expected_active
        @test size(epienv.habitat.matrix) == expected_grid
        @test epienv.habitat.size ≈ expected_gridlength
        @test epienv.habitat.matrix ≈ expected_matrix
    end

    @testset "Construct from initial population" begin
        # Set an initial population of zeros
        # The zeros should not be masked out (but the NaNs should)
        initial_pop = zeros(grid...)
        initial_pop[.!active] .= NaN

        epienv = simplehabitatAE(fillval, area, control, initial_pop)
        expected_initial_pop = zeros(expected_grid)

        @test epienv.active == expected_active
        @test size(epienv.habitat.matrix) == expected_grid
        @test epienv.habitat.size ≈ expected_gridlength
        @test epienv.habitat.matrix ≈ expected_matrix
        @test epienv.initial_pop ≈ expected_initial_pop
    end
end
