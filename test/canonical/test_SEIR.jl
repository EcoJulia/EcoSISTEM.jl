using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Test

# Set simulation parameters
birth = fill(0.0/day, numclasses)
death = fill(0.0/day, numclasses)
beta = 5.0days
mu = 1/7days
sigma = 1/7days
virus_growth = 1e-3/day
virus_decay = Inf/day
param = SEIRGrowth{typeof(unit(beta))}(birth, death, virus_growth, virus_decay, beta, mu, sigma)
param = transition(param)

# Set up simple gridded environment
grid = (4, 4)
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
virus = 10_000
susceptible = 5_000_000
exposed = 0
infected = 1_000
recovered = 0
dead = 0
abun = [virus, susceptible, exposed, infected, recovered, dead]

# Dispersal kernels for virus and disease classes
dispersal_dists = fill(100.0km, numclasses)
dispersal_dists[3] = 700.0km
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numclasses), fill(0.1K, numclasses))
epilist = SEIR(traits, abun, movement, param)

# Create epi system with all information
rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)

# Run simulation
abuns = zeros(Int64, numclasses, 16, 731)
times = 2years; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

# Test no-one dies (death rate = 0)
@test sum(thisabun[end, :, :]) == 0
# Test overall population size stays constant (birth rate = death rate = 0)
@test all(sum(thisabun[2:3, :, :], dims = (1, 2)) .== (susceptible + infected))
