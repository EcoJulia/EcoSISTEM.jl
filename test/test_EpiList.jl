using Simulation
using Test
using Unitful.DefaultSymbols
using Distributions
using Simulation.Units


birth = [0.0/day; fill(1e-5/day, 3); 0.0/day]
death = [0.0/day; fill(1e-5/day, 3); 0.0/day]
beta = 0.0005/day
sigma = 0.05/day
virus_growth = 0.0001/day
virus_decay = 0.07/day
param = SIRGrowth{typeof(unit(beta))}(birth, death, virus_growth, virus_decay, beta, sigma)
param = transition(param)

abun = [10, 1000, 1, 0, 0]
dispersal_dists = [1e-2km; fill(2.0km, 3); 1e-2km]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

traits = GaussTrait(fill(298.0K, 5), fill(0.1K, 5))
@test_nowarn SIR(traits, abun, movement, param)
epilist = SIR(traits, abun, movement, param)
@test epilist.names[1] == "Virus"
@test epilist.names[2] == "Susceptible"
@test epilist.names[3] == "Infected"
@test epilist.names[4] == "Recovered"
@test epilist.names[5] == "Dead"

@test length(epilist.names) == length(epilist.movement.kernels)
@test length(epilist.names) == length(epilist.abun)
@test length(epilist.names) == length(epilist.traits.mean)
@test length(epilist.names) == length(epilist.traits.var)
