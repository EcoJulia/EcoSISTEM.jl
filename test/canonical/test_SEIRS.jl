using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Test

# sort out settings to potentially save inputs/outputs of `simulate`
do_save = (@isdefined do_save) ? do_save : false
save_path = (@isdefined save_path) ? save_path : pwd()

numclasses = 5
numvirus = 2
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
abun_h = (
  Susceptible = 500_000 * prod(grid),
  Exposed = 0,
  Infected = 100 * prod(grid),
  Recovered = 0,
  Dead = 0
)
disease_classes = (
  susceptible = ["Susceptible"],
  infectious = ["Infected"]
)
abun_v = (Environment = 0, Force = 0)
# Dispersal kernels for virus and disease classes
dispersal_dists = fill(700.0km, prod(grid))
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = EpiMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, param)

# Create epi system with all information
rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)

# Run simulation
times = 2years; interval = 1day; timestep = 1day
abuns = zeros(Int64, size(epi.abundances.matrix, 1), grid[1]*grid[2], convert(Int64, floor(times / interval)) + 1)
@time simulate_record!(abuns, epi, times, interval, timestep; save=do_save, save_path=save_path)

# Test no-one dies (death rate = 0)
@test sum(abuns[end, :, :]) == abun_h.Dead
# Test overall population size stays constant (birth rate = death rate = 0)
@test all(sum(abuns[1:4, :, :], dims = (1, 2)) .== (abun_h.Susceptible + abun_h.Infected))

### TEST OUTPUTS
# TODO: When shifting out virus from Epilist, these indexes will need updating.

idx_sus = 1
idx_rec = 4
idx_dead = 5

# Test susceptible population decreasing or constant only [Source]
# https://github.com/ScottishCovidResponse/Simulation.jl/pull/37
@test sum(abuns[idx_sus, :, 1]) == abun_h.Susceptible

@test sum(abuns[idx_rec, :, 1]) == abun_h.Recovered

# Test dead population increasing or constant only [Sink]
@test all(diff(vec(sum(abuns[idx_dead, :, :], dims = 1))) .>= 0)
@test sum(abuns[idx_dead, :, 1]) == abun_h.Dead
