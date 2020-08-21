using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Test
using DataFrames

@testset "SEIRS" begin

# sort out settings to potentially save inputs/outputs of `simulate`
do_save = (@isdefined do_save) ? do_save : false
save_path = (@isdefined save_path) ? save_path : pwd()

# Set up simple gridded environment
grid = (8, 8)
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all pathogen categories
abun_v = DataFrame([
    (name="Environment", initial=0),
    (name="Force", initial=0),
])
numvirus = nrow(abun_v)

# Set initial population sizes for all human categories
susceptible = 500_000 * prod(grid)
exposed = 0
infected = 100 * prod(grid)
abun_h = DataFrame([
    (name="Susceptible", type=Susceptible, initial=susceptible),
    (name="Exposed", type=OtherDiseaseState, initial=exposed),
    (name="Infected", type=Infectious, initial=infected),
    (name="Recovered", type=Removed, initial=0),
    (name="Dead", type=Removed, initial=0),
])
numclasses = nrow(abun_h)

# Set non-pathogen mediated transitions
mu = 1/7days
sigma = 0.1/days
epsilon = 0.01/day
transitions = DataFrame([
  (from="Exposed", to="Infected", prob=mu),
  (from="Infected", to="Recovered", prob=sigma),
  (from="Recovered", to="Susceptible", prob=epsilon),
])


# Set simulation parameters
birth = fill(0.0/day, numclasses)
death = fill(0.0/day, numclasses)
beta_force = 1.0/day
beta_env = 1.0/day
virus_growth = 1.0/day
virus_decay = 1.0/2days
param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force)

# Dispersal kernels for virus and disease classes
dispersal_dists = fill(700.0km, prod(grid))
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = EpiMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param)

# Create epi system with all information
rel = Gauss{eltype(epienv.habitat)}()
epi = EpiSystem(epilist, epienv, rel)

# Run simulation
times = 2years; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), div(times, interval) + 1)
@time simulate_record!(abuns, epi, times, interval, timestep; save=do_save, save_path=save_path)

# Find correct indices in arrays
idx_sus = findfirst(==("Susceptible"), abun_h.name)
idx_rec = findfirst(==("Recovered"), abun_h.name)
idx_dead = findfirst(==("Dead"), abun_h.name)

# Test no-one dies (death rate = 0)
@test all(sum(abuns[idx_dead, :, :], dims = (1)) .== abun_h.initial[idx_dead])
# Test overall population size stays constant (birth rate = death rate = 0)
@test all(sum(abuns[1:4, :, :], dims = (1, 2)) .== sum(abun_h.initial))

### TEST OUTPUTS

# Test susceptible population decreasing or constant only [Source]
# https://github.com/ScottishCovidResponse/Simulation.jl/pull/37
@test sum(abuns[idx_sus, :, 1]) == abun_h.initial[idx_sus]

@test sum(abuns[idx_rec, :, 1]) == abun_h.initial[idx_rec]

# Test dead population increasing or constant only [Sink]
@test all(diff(vec(sum(abuns[idx_dead, :, :], dims = 1))) .>= 0)
@test sum(abuns[idx_dead, :, 1]) == abun_h.initial[idx_dead]
end
