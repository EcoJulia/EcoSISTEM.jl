using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

# Set simulation parameters
birth = [0.0/day; fill(1e-5/day, 2); 0.0/day]
death = [0.0/day; fill(1e-5/day, 2); 0.0/day]
beta = 0.0005/day
sigma = 0.05/day
virus_growth = 0.0001/day
virus_decay = 0.07/day
param = SISGrowth{typeof(unit(beta))}(birth, death, virus_growth, virus_decay, beta, sigma)
param = transition(param)

# Set up simple gridded environment
grid = (2, 2)
area = 10.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
abun = [10, 1000, 1, 0]

# Dispersal kernels for virus and disease classes
dispersal_dists = [1e-2km; fill(2.0km, 2); 1e-2km]
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
