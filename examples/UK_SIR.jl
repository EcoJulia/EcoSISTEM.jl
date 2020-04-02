using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

# Larger grid
birth = [0.001/day; fill(1e-10/day, 3)]
death = [0.001/day; fill(1e-10/day, 3)]
beta = 0.14/day
sigma = 1.0/14days
param = EpiGrowth{typeof(unit(beta))}(birth, death, beta, sigma)

grid = (500, 500)
area = 250_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

abun = [0, 60_000_000, 0, 0]

dispersal_dists = [2.0km; fill(0.01km, 3)]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

traits = GaussTrait(fill(298.0K, 4), fill(0.1K, 4))
epilist = SIR(traits, abun, movement, param)

rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)
epi.abundances.matrix[1, 1] = 10_000
epi.abundances.matrix[3, 1] = 100

abuns = zeros(Int64, 4, 250_000, 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

display(heatmap(reshape(abuns[3, :, 1], 500, 500), layout = (@layout [a b; c d]), subplot = 1, title = "Day 1", clim = (0, 100)))
display(heatmap!(reshape(abuns[3, :, 7], 500, 500), subplot = 2, title = "3 months", clim = (0, 100)))
display(heatmap!(reshape(abuns[3, :, 14], 500, 500), subplot = 3, title = "6 months", clim = (0, 100)))
display(heatmap!(reshape(abuns[3, :, 30], 500, 500), subplot = 4, title = "1 year", clim = (0, 100)))
