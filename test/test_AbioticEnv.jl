using Simulation
using Unitful.DefaultSymbols
using Base.Test

grid = (5, 5)
area = 25.0km^2
totalK = 10000.0
numNiches = 4

@test_nowarn abenv = simplehabitatAE(0.0, grid, totalK, area)
@test_nowarn abenv = tempgradAE(-10.0°C, 10.0°C, grid, totalK, area,
    0.01°C/month)
@test_nowarn abenv = simplenicheAE(numNiches, grid, totalK, area)
