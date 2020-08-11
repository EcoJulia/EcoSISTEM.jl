using Simulation
using Unitful.DefaultSymbols
using Test
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

@test Simulation._getsubcommunitynames(abenv) == abenv.names
@test Simulation.getavailableenergy(abenv) == sum(abenv.budget.matrix)


# TEST tempgradAE
abenv = tempgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K/month)
@test_nowarn tempgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K/month)
@test minimum(abenv.habitat.matrix) == -10.0K
@test maximum(abenv.habitat.matrix) == 10.0K
@test size(abenv.habitat.matrix) == grid
@test sum(abenv.budget.matrix) == totalK * area
@test abenv.active == active
@test all(abenv.active)

# TEST peakedgradAE
abenv = peakedgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K/month)
@test_nowarn peakedgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K/month)
@test minimum(abenv.habitat.matrix) == -10.0K
@test maximum(abenv.habitat.matrix) == 10.0K
@test size(abenv.habitat.matrix) == grid
@test sum(abenv.budget.matrix) == totalK * area
@test abenv.active == active
@test all(abenv.active)

# TEST raingradAE
abenv = raingradAE(0.0mm, 100.0mm, grid, totalK, area, 0.01mm/month)
@test_nowarn raingradAE(0.0mm, 100.0mm, grid, totalK, area, 0.01mm/month)
@test minimum(abenv.habitat.matrix) == 0.0mm
@test maximum(abenv.habitat.matrix) == 100.0mm
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

# TEST eraAE
temp = AxisArray(fill(1.0K, 10, 10, 3), Axis{:latitude}(1:10), Axis{:longitude}(1:10), Axis{:time}(collect(1:3) .* s))
eratemp = ERA(temp)
active = fill(true, 10, 10)
totalK = 1000.0kJ/m^2; area = 100.0km^2
ea = eraAE(eratemp, totalK, area)
@test_nowarn eraAE(eratemp, totalK, area)
@test_nowarn eraAE(eratemp, totalK, area, active)
@test size(ea.habitat.matrix) == size(temp)
@test sum(ea.budget.matrix) == totalK * area
@test ea.active == active
@test all(ea.active)

solar = SolarTimeBudget(fill(10.0kJ, 10, 10, 3), 1)
ea = eraAE(eratemp, solar, active)

# TEST worldclimAE
temp = AxisArray(fill(1.0K, 10, 10, 12), Axis{:latitude}(collect(1:10) .* m), Axis{:longitude}(collect(1:10) .* m), Axis{:time}(collect(1:12) .* month))
wctemp = Worldclim(temp)
active = fill(true, 10, 10)
totalK = 1000.0kJ/m^2; area = 100.0km^2
wc = worldclimAE(wctemp, totalK, area)
@test_nowarn worldclimAE(wctemp, totalK, area)
@test_nowarn worldclimAE(wctemp, totalK, area, active)
@test size(wc.habitat.matrix) == size(temp)
@test sum(wc.budget.matrix) == totalK * area
@test wc.active == active
@test all(wc.active)

solar = SolarTimeBudget(fill(10.0kJ, 10, 10, 3), 1)
wc = worldclimAE(wctemp, solar, active)
