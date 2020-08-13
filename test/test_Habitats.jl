using Simulation
using Unitful.DefaultSymbols
using Test
using Simulation.Units
using Simulation.ClimatePref
using AxisArrays

grid = (10, 10)
area = 25.0km^2
totalK = 10000.0kJ/km^2
numNiches = 4
active = fill(true, grid)

# TEST simplehabitatAE
fillval = 0.0
abenv = simplehabitatAE(fillval, grid, totalK, area)
@test Simulation.iscontinuous(abenv.habitat) == true
@test eltype(abenv.habitat) == typeof(abenv.habitat.matrix[1])
@test size(abenv.habitat, 1) == grid[1]
@test Simulation.xmin(abenv.habitat) == 0
@test Simulation.ymin(abenv.habitat) == 0
@test Simulation.xcellsize(abenv.habitat) == sqrt(area/prod(grid))/km
@test Simulation.ycellsize(abenv.habitat) == sqrt(area/prod(grid))/km
@test Simulation.xcells(abenv.habitat) == size(abenv.habitat, 1)
@test Simulation.ycells(abenv.habitat) == size(abenv.habitat, 2)
@test Simulation.indices(abenv.habitat) == Simulation.coordinates(abenv.habitat)
@test Simulation.indices(abenv.habitat, 1) == repeat(collect(1:grid[1]), grid[2])

# TEST tempgradAE
abenv = tempgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K/month)
@test Simulation.iscontinuous(abenv.habitat) == true
@test eltype(abenv.habitat) == typeof(abenv.habitat.matrix[1])
@test size(abenv.habitat, 1) == grid[1]

# TEST simplenicheAE
abenv = simplenicheAE(numNiches, grid, totalK, area)
@test Simulation.iscontinuous(abenv.habitat) == false
@test eltype(abenv.habitat) == typeof(abenv.habitat.matrix[1])
@test size(abenv.habitat, 1) == grid[1]

# TEST eraAE
temp = AxisArray(fill(1.0K, 10, 10, 3), Axis{:latitude}(1:10), Axis{:longitude}(1:10), Axis{:time}(collect(1:3) .* s))
eratemp = ERA(temp)
active = fill(true, 10, 10)
solar = SolarTimeBudget(fill(10.0kJ, 10, 10, 3), 1)
ea = eraAE(eratemp, solar, active)
@test Simulation.iscontinuous(ea.habitat) == true
@test eltype(ea.habitat) == typeof(ea.habitat.matrix[1])
ea.habitat.time = 2
Simulation._resettime!(ea.habitat)
@test ea.habitat.time == 1
@test size(ea.habitat, 1) == grid[1]

# TEST worldclimAE
temp = AxisArray(fill(1.0K, 10, 10, 12), Axis{:latitude}(collect(1:10) .* m), Axis{:longitude}(collect(1:10) .* m), Axis{:time}(collect(1:12) .* month))
wctemp = Worldclim(temp)
active = fill(true, 10, 10)
solar = SolarTimeBudget(fill(10.0kJ, 10, 10, 3), 1)
wc = worldclimAE(wctemp, solar, active)
@test Simulation.iscontinuous(wc.habitat) == true
@test eltype(wc.habitat) == typeof(wc.habitat.matrix[1])
wc.habitat.time = 2
Simulation._resettime!(wc.habitat)
@test wc.habitat.time == 1
@test size(wc.habitat, 1) == grid[1]

# Test multi habitats
hab = HabitatCollection2(wc.habitat, ea.habitat)
@test Simulation.iscontinuous(hab) == [true, true]
@test eltype(hab) == [typeof(hab.h1.matrix[1]), typeof(hab.h2.matrix[1])]
@test size(hab, 1) == grid[1]
hab.h1.time = 2
Simulation._resettime!(hab)
@test hab.h1.time == 1
@test Simulation._getgridsize(hab) == hab.h1.size
@test isapprox(Simulation._getsize(hab), hab.h1.size^2 * prod(grid), rtol = 1e-5)
@test Simulation._getdimension(hab) == grid
@test Simulation._countsubcommunities(hab) == prod(grid)

hab = HabitatCollection3(abenv.habitat, wc.habitat, ea.habitat)
@test Simulation.iscontinuous(hab) == [false, true, true]
@test eltype(hab) == [typeof(hab.h1.matrix[1]), typeof(hab.h2.matrix[1]), typeof(hab.h3.matrix[1])]
@test size(hab, 1) == grid[1]
@test Simulation._getgridsize(hab) == hab.h1.size
@test isapprox(Simulation._getsize(hab), hab.h1.size^2 * prod(grid), rtol = 1e-5)
@test Simulation._getdimension(hab) == grid
@test Simulation._countsubcommunities(hab) == prod(grid)

hab = HabitatCollection3(wc.habitat, wc.habitat, ea.habitat)
hab.h2.time = 2
Simulation._resettime!(hab)
@test hab.h2.time == 1
