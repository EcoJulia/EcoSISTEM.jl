using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Test

# Set simulation parameters
birth = fill(0.0/day, 4)
death = fill(0.0/day, 4)
beta = 0.01/day
sigma = 0.02/day
virus_growth = 1e-3/day
virus_decay = 1.0/day
param = SISGrowth{typeof(unit(beta))}(birth, death, virus_growth, virus_decay, beta, sigma)
param = transition(param)

# Set up simple gridded environment
grid = (4, 4)
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
virus = 10
susceptible = 5_000_000
infected = 1
dead = 0
abun = [virus, susceptible, infected, dead]

# Dispersal kernels for virus and disease classes
dispersal_dists = fill(1.0km, 8)
dispersal_dists[3] = 50.0km
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, 4), fill(0.1K, 4))
epilist = SIS(traits, abun, movement, param)

# Create epi system with all information
rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)

# Run simulation
abuns = zeros(Int64, size(epi.abundances.matrix, 1), size(epi.abundances.matrix, 2), 731)
times = 2years; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

@test sum(abuns[end, :, :]) == 0
@test all(sum(abuns[2:3, :, :], dims = (1, 2)) .== (susceptible + infected))
