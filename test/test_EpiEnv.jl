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
