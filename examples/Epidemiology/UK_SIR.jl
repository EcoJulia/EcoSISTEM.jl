using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

# Set simulation parameters
birth = [0.1/day; fill(1e-5/day, 3)]
death = [1/day; fill(1e-5/day, 3)]
beta = 1e-2/day
sigma = 1/14days
param = EpiGrowth{typeof(unit(beta))}(birth, death, beta, sigma)

# Set up simple gridded environment
grid = (500, 500)
area = 250_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
abun = [0, 60_000_000, 0, 0]

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

# Add in initial infections randomly across grid
samp = rand(1:250_000, 100)
epi.abundances.matrix[1, samp] .= 100
epi.abundances.matrix[3, samp] .= 10

# Run simulation
abuns = zeros(Int64, 4, 250_000, 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

# Plot heatmap of spatial disease progression
using Plots
plotlyjs()
display(heatmap(reshape(abuns[3, :, 1], 500, 500), layout = (@layout [a b; c d]), subplot = 1, title = "Day 1", clim = (0, 200)))
display(heatmap!(reshape(abuns[3, :, 7], 500, 500), subplot = 2, title = "Day 7", clim = (0, 200)))
display(heatmap!(reshape(abuns[3, :, 14], 500, 500), subplot = 3, title = "Day 14", clim = (0, 200)))
display(heatmap!(reshape(abuns[3, :, 30], 500, 500), subplot = 4, title = "Day 30", clim = (0, 200)))

# View summed SIR dynamics for whole area
display(plot(mapslices(sum, abuns[2, :, :], dims = 1)[1, :], color = :Blue, label = "Susceptible"))
display(plot!(mapslices(sum, abuns[3, :, :], dims = 1)[1, :], color = :Red, label = "Infected"))
display(plot!(mapslices(sum, abuns[4, :, :], dims = 1)[1, :], color = :Black, label = "Recovered"))
