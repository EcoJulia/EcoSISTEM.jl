using Simulation
using Unitful.DefaultSymbols
using Test
using Simulation.Units
import Simulation.EpiGrowth

birth = [0.0001/day; fill(1e-10/day, 3)]
death = [0.07/day; fill(1e-10/day, 3)]
beta = 0.05/day
sigma = 0.05/day
@test_nowarn EpiGrowth{typeof(unit(beta))}(birth, death, beta, sigma)
params = EpiGrowth{typeof(unit(beta))}(birth, death, beta, sigma)
@test length(params.birth) == length(params.birth)
@test params.birth == birth
@test params.death == death
@test params.beta == beta
@test params.sigma == sigma
