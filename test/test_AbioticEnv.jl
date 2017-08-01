using Simulation
using Base.Test

grid = (5,5)
gridSize = 1.0
totalK = 1000.0
numNiches = 4

@test_nowarn abenv = simplehabitatAE(0.0, grid, totalK, gridSize)
@test_nowarn abenv = tempgradAE(-10.0, 10.0, grid, totalK, gridSize, 0.01)
@test_nowarn abenv = simplenicheAE(numNiches, grid, totalK, gridSize)
