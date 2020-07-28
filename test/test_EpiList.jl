using Simulation
using Test
using Unitful.DefaultSymbols
using Distributions
using Simulation.Units

numvirus = 2
numclasses = 4
birth = [0.0/day; fill(1e-5/day, 3); 0.0/day]
death = [0.0/day; fill(1e-5/day, 3); 0.0/day]
beta_force = 5.0/day
beta_env = 0.5/day
sigma = 0.05/day
virus_growth = 0.0001/day
virus_decay = 0.07/day
param = SIRGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
param = transition(param)

sus = ["Susceptible"]
inf = ["Infected"]
abun_h = (
    Susceptible = 1000,
    Infected = 1,
    Recovered = 0,
    Dead = 0
)
disease_classes = (
    susceptible = ["Susceptible"],
    infectious = ["Infected"]
)
abun_v = (Environment = 10, Force = 10)

dispersal_dists = fill(2.0km, numclasses)
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = EpiMovement(kernel)

traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
@test_nowarn EpiList(traits, abun_v, abun_h, disease_classes, movement, param)
epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, param)
@test epilist.virus.names[1] == "Environment"
@test epilist.virus.names[2] == "Force"
@test epilist.human.names[1] == "Susceptible"
@test epilist.human.names[2] == "Infected"
@test epilist.human.names[3] == "Recovered"
@test epilist.human.names[4] == "Dead"

@test length(epilist.human.names) == length(epilist.human.abun)
@test length(epilist.virus.names) == length(epilist.virus.abun)
@test length(epilist.virus.names) == length(epilist.virus.traits.mean)
@test length(epilist.virus.names) == length(epilist.virus.traits.var)
