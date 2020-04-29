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
