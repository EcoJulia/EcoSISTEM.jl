using Simulation
using Test
using Unitful.DefaultSymbols
using Distributions
using Simulation.Units


birth = [0.0/day; fill(1e-5/day, 3); 0.0/day]
death = [0.0/day; fill(1e-5/day, 3); 0.0/day]
beta_force = 5.0/day
beta_env = 0.5/day
sigma = 0.05/day
virus_growth = 0.0001/day
virus_decay = 0.07/day
param = SIRGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
param = transition(param)

abun_h = [1000, 1, 0, 0]
abun_v = [10]
dispersal_dists = [fill(2.0km, 3); 1e-2km]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

traits = GaussTrait([298.0K], [0.1K])
@test_nowarn SIR(traits, abun_v, abun_h, movement, param)
epilist = SIR(traits, abun_v, abun_h, movement, param)
@test epilist.virus.names[1] == "Virus"
@test epilist.human.names[1] == "Susceptible"
@test epilist.human.names[2] == "Infected"
@test epilist.human.names[3] == "Recovered"
@test epilist.human.names[4] == "Dead"

@test length(epilist.human.names) == length(epilist.human.movement.kernels)
@test length(epilist.human.names) == length(epilist.human.abun)
@test length(epilist.virus.names) == length(epilist.virus.abun)
@test length(epilist.virus.names) == length(epilist.virus.traits.mean)
@test length(epilist.virus.names) == length(epilist.virus.traits.var)
