using Simulation
using Test
using Unitful.DefaultSymbols
using Distributions
using Simulation.Units


birth = [0.0001/day; fill(1e-10/day, 3)]
death = [0.07/day; fill(1e-10/day, 3)]
beta = 0.05/day
sigma = 0.05/day
param = SIRGrowth{typeof(unit(beta))}(birth, death, beta, sigma)

abun = [10, 1000, 1, 0]
dispersal_dists = [2.0km; fill(0.01km, 3)]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

traits = GaussTrait(fill(298.0K, 4), fill(0.1K, 4))
@test_nowarn SIR(traits, abun, movement, param)
epilist = SIR(traits, abun, movement, param)
@test epilist.names[1] == "Virus"
@test epilist.names[2] == "Susceptible"
@test epilist.names[3] == "Infected"
@test epilist.names[4] == "Recovered"

@test length(epilist.names) == length(epilist.movement.kernels)
@test length(epilist.names) == length(epilist.abun)
@test length(epilist.names) == length(epilist.traits.mean)
@test length(epilist.names) == length(epilist.traits.var)
