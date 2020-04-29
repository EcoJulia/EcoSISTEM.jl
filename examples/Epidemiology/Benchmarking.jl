using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

birth = [0.0/day; fill(1e-5/day, 3); 0.0/day]
death = [0.0/day; fill(1e-5/day, 3); 0.0/day]
beta = 0.05/day
sigma = 0.05/day
virus_growth = 0.0001/day
virus_decay = 0.07/day
param = SIRGrowth{typeof(unit(beta))}(birth, death, virus_growth, virus_decay, beta, sigma)
param = transition(param)

grid = (10, 10)
area = 100.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

abun = [10, 1000, 1, 0, 0]

dispersal_dists = [1e-2km; fill(2.0km, 3); 1e-2km]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

traits = GaussTrait(fill(298.0K, 5), fill(0.1K, 5))
epilist = SIR(traits, abun, movement, param)

rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)

times = 1year; interval = 1day; timestep = 1day
@time simulate!(epi, times, timestep)
