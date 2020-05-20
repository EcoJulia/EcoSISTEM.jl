using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Test

# sort out settings to potentially save inputs/outputs of `simulate`
do_save = (@isdefined do_save) ? do_save : false
save_path = (@isdefined save_path) ? save_path : pwd()

numclasses = 5
# Set simulation parameters
birth = fill(0.0/day, numclasses)
death = fill(0.0/day, numclasses)
beta_force = 1.0/day
beta_env = 1.0/day
mu = 1/7days
sigma = 0.1/days
epsilon = 0.01/day
virus_growth = 1.0/day
virus_decay = 1.0/2days
param = SEIRSGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, mu, sigma, epsilon)
param = transition(param)

# Set up simple gridded environment
grid = (8, 8)
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
virus = 0
susceptible = 500_000 * prod(grid)
exposed = 0
infected = 100 * prod(grid)
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
epilist = SEIRS(traits, abun, movement, param)

# Create epi system with all information
rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)

# Run simulation
times = 2years; interval = 1day; timestep = 1day
abuns = zeros(Int64, size(epi.abundances.matrix, 1), grid[1]*grid[2], convert(Int64, floor(times / interval)) + 1)
@time simulate_record!(abuns, epi, times, interval, timestep; save=do_save, save_path=save_path)

# Test no-one dies (death rate = 0)
@test sum(abuns[end, :, :]) == 0
# Test overall population size stays constant (birth rate = death rate = 0)
@test all(sum(abuns[2:5, :, :], dims = (1, 2)) .== (susceptible + infected))
