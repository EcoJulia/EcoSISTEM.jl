using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units
import Simulation: emptyepilandscape

include("TestCases.jl")

epi = TestEpiSystem()

@test_nowarn emptyepilandscape(epi.epienv, epi.epilist, Int64(1))
el = emptyepilandscape(epi.epienv, epi.epilist, Int64(1))
@test size(el.grid, 2) * size(el.grid, 3) == size(el.matrix, 2)
@test size(el.grid, 1) == size(el.matrix, 1)
@test sum(el.grid) == sum(el.matrix)
el.matrix[1, 1] += 1
@test sum(el.grid) == sum(el.matrix)
el.grid[1, 1, 1] += 1
@test sum(el.grid) == sum(el.matrix)
