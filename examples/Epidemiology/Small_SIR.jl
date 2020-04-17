using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Plots
plotlyjs()

# Set simulation parameters
birth = [0.0001/day; fill(1e-10/day, 3)]
death = [0.07/day; fill(1e-10/day, 3)]
beta = 0.05/day
sigma = 0.05/day
param = EpiGrowth{typeof(unit(beta))}(birth, death, beta, sigma)

# Set up simple gridded environment
grid = (2, 2)
area = 10.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
abun = [10, 1000, 1, 0]

# Dispersal kernels for virus and disease classes
dispersal_dists = [2.0km; fill(0.01km, 3)]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, 4), fill(0.1K, 4))
epilist = SIR(traits, abun, movement, param)

# Create epi system with all information
rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)

# Run simulation
abuns = zeros(Int64, 4, 4, 731)
times = 2years; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

# Plot SIR dynamics
display(plot(mapslices(sum, abuns[2, :, :], dims = 1)[1, :], color = :Blue, label = "Susceptibles"))
display(plot!(mapslices(sum, abuns[3, :, :], dims = 1)[1, :], color = :Red, label = "Infecteds"))
display(plot!(mapslices(sum, abuns[4, :, :], dims = 1)[1, :], color = :Black, label = "Recovereds"))


# Again, but with larger grid
birth = [0.0/day; fill(1e-10/day, 3)]
death = [0.001/day; fill(1e-10/day, 3)]
beta = 0.1/day
sigma = 0.005/day
viralload = 1.0/day
param = EpiGrowth{typeof(unit(beta))}(birth, death, beta, sigma, viralload)

grid = (10, 10)
area = 100.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

abun = [0, 10_000, 0, 0]

dispersal_dists = [0.5km; fill(0.01km, 3)]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

traits = GaussTrait(fill(298.0K, 4), fill(0.1K, 4))
epilist = SIR(traits, abun, movement, param)

rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)
epi.abundances.matrix[1, 1] = 100
epi.abundances.matrix[3, 1] = 1

abuns = zeros(Int64, 4, 100, 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

# Plot as heatmaps of infected individuals over space
display(heatmap(reshape(abuns[3, :, 1], 10, 10), layout = (@layout [a b; c d]), subplot = 1, title = "Day 0", clim = (0, 120)))
display(heatmap!(reshape(abuns[3, :, 7], 10, 10), subplot = 2, title = "Day 7", clim = (0, 120)))
display(heatmap!(reshape(abuns[3, :, 14], 10, 10), subplot = 3, title = "Day 14", clim = (0, 120)))
display(heatmap!(reshape(abuns[3, :, 30], 10, 10), subplot = 4, title = "Day 30", clim = (0, 120)))
