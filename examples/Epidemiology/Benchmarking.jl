using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

beta = 0.05/day
sigma = 0.05/day
birth = [0.001/day; fill(1e-10/day, 3)]
death = [0.001/day; fill(1e-10/day, 3)]
param = SIRGrowth{typeof(unit(beta))}(birth, death, beta, sigma)

grid = (10, 10)
area = 100.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

abun = [10, 1000, 1, 0]

dispersal_dists = [2.0km; fill(0.01km, 3)]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

traits = GaussTrait(fill(298.0K, 4), fill(0.1K, 4))
epilist = SIR(traits, abun, movement, param)

rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)

times = 1year; interval = 1day; timestep = 1day
@time simulate!(epi, times, timestep)
