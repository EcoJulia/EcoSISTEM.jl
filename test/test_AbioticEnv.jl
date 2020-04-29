using Simulation
using Unitful.DefaultSymbols
using Compat.Test
using Simulation.Units
using Simulation.ClimatePref
using AxisArrays

grid = (5, 5)
area = 25.0km^2
totalK = 10000.0kJ/km^2
numNiches = 4
active = fill(true, grid)

# TEST simplehabitatAE
fillval = 0.0
abenv = simplehabitatAE(fillval, grid, totalK, area)
@test_nowarn simplehabitatAE(fillval, grid, totalK, area)
@test_nowarn simplehabitatAE(fillval, grid, totalK, area, active)
@test all(abenv.habitat.matrix .== fillval)
@test size(abenv.habitat.matrix) == grid
@test sum(abenv.budget.matrix) == totalK * area
@test abenv.active == active
@test all(abenv.active)

# TEST tempgradAE
abenv = tempgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K/month)
@test_nowarn tempgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K/month)
@test minimum(abenv.habitat.matrix) == -10.0K
@test maximum(abenv.habitat.matrix) == 10.0K
@test size(abenv.habitat.matrix) == grid
@test sum(abenv.budget.matrix) == totalK * area
@test abenv.active == active
@test all(abenv.active)

# TEST simplenicheAE
abenv = simplenicheAE(numNiches, grid, totalK, area)
@test_nowarn simplenicheAE(numNiches, grid, totalK, area)
@test_nowarn simplenicheAE(numNiches, grid, totalK, area, active)
@test maximum(abenv.habitat.matrix) <= numNiches
@test size(abenv.habitat.matrix) == grid
@test sum(abenv.budget.matrix) == totalK * area
@test abenv.active == active
@test all(abenv.active)
