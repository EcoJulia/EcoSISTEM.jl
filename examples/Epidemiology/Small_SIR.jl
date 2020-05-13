using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Plots
plotlyjs()

do_plot = false

# Set simulation parameters
birth = [0.0/day; fill(1e-5/day, 3); 0.0/day]
death = [0.0/day; fill(1e-5/day, 3); 0.0/day]
beta_force = 0.0005/day
beta_env = 0.0005/day
sigma = 0.05/day
virus_growth = 0.0001/day
virus_decay = 0.07/day
param = SIRGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
param = transition(param)

# Set up simple gridded environment
grid = (2, 2)
area = 10.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
abun = [10, 1000, 1, 0, 0]

# Dispersal kernels for virus and disease classes
dispersal_dists = [1e-2km; fill(2.0km, 3); 1e-2km]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, 5), fill(0.1K, 5))
epilist = SIR(traits, abun, movement, param)

# Create epi system with all information
rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)

# Run simulation
abuns = zeros(Int64, 5, 4, 731)
times = 2years; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

if do_plot
    # Plot SIR dynamics
    display(plot(mapslices(sum, abuns[2, :, :], dims = 1)[1, :], color = :Blue, label = "Susceptibles"))
    display(plot!(mapslices(sum, abuns[3, :, :], dims = 1)[1, :], color = :Red, label = "Infecteds"))
    display(plot!(mapslices(sum, abuns[4, :, :], dims = 1)[1, :], color = :Black, label = "Recovereds"))
end

# Again, but with larger grid
birth = [0.0/day; fill(1e-5/day, 3); 0.0/day]
death = [0.0/day; fill(1e-5/day, 3); 0.0/day]
beta_force = 0.05/day
beta_env = 0.05/day
sigma = 0.05/day
virus_growth = 0.0001/day
virus_decay = 0.07/day
param = SIRGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
param = transition(param)

grid = (10, 10)
area = 100.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

abun = [0, 10_000, 0, 0, 0]

dispersal_dists = [1e-2km; fill(0.5km, 3); 1e-2km]
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

traits = GaussTrait(fill(298.0K, 5), fill(0.1K, 5))
epilist = SIR(traits, abun, movement, param)

rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)
epi.abundances.matrix[1, 1] = 100
epi.abundances.matrix[3, 1] = 1

abuns = zeros(Int64, 5, 100, 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

if do_plot
    # Plot as heatmaps of infected individuals over space
    display(heatmap(reshape(abuns[3, :, 1], 10, 10), layout = (@layout [a b; c d]), subplot = 1, title = "Day 0", clim = (0, 120)))
    display(heatmap!(reshape(abuns[3, :, 7], 10, 10), subplot = 2, title = "Day 7", clim = (0, 120)))
    display(heatmap!(reshape(abuns[3, :, 14], 10, 10), subplot = 3, title = "Day 14", clim = (0, 120)))
    display(heatmap!(reshape(abuns[3, :, 30], 10, 10), subplot = 4, title = "Day 30", clim = (0, 120)))
end
