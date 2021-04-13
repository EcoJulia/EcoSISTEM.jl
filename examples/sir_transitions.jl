using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using DataFrames
using Plots

# Set up simple gridded environment
grid = (10, 10)
area = 1_000.0km^2
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
sigma = 1/7days

# Set simulation parameters
birth = fill(0.0/day, numclasses)
death = fill(0.0/day, numclasses)
beta_force = 10.0/day
beta_env = 10.0/day
virus_growth = 0.1/day
virus_decay = 1.0/2days
param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force)

transitions = DataFrame([
  (from="Susceptible", from_id=[1], to="Exposed", to_id=[2], prob = (env = beta_env, force = beta_force)),
  (from="Exposed", from_id=[2], to="Infected", to_id=[3], prob=[mu]),
  (from="Infected", from_id=[3], to="Recovered", to_id=[4], prob=[sigma]),
])

# Dispersal kernels for virus and disease classes
dispersal_dists = fill(10.0km, prod(grid))
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = EpiMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param)

# Create epi system with all information
rel = Gauss{eltype(epienv.habitat)}()
transitions = create_transition_list_SEIR(epilist, epienv)
epi = Ecosystem(epilist, epienv, rel, transitions = transitions)

# Run simulation
times = 1month; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), floor(Int, times/timestep) + 1)
@time simulate_record!(abuns, epi, times, interval, timestep);
plot_epidynamics(epi, abuns)

epi = Ecosystem(epilist, epienv, rel, transitions = transitions)

# Run simulation
times = 1month; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), floor(Int, times/timestep) + 1)
@time epi_simulate_record!(abuns, epi, times, interval, timestep);
plot_epidynamics(epi, abuns)

# Simulation Parameters
burnin = 1month; times = 5years; timestep = 1day; record_interval = 1day; repeats = 1
lensim = length(0years:record_interval:times)
abuns = zeros(Int64, numSpecies, prod(grid), lensim)
# Burnin
@time simulate!(epi, burnin, timestep);
@time simulate_record!(abuns, epi, times, record_interval, timestep);

# Benchmark
using BenchmarkTools
epi = Episystem(sppl, abenv, rel)
@benchmark simulate!(epi, burnin, timestep)
epi = Episystem(sppl, abenv, rel)
@benchmark wpi_simulate!(epi, burnin, timestep)

using Plots
@gif for i in 1:lensim
    heatmap(abuns[:, :, i], clim = (50, 150))
end

using ProfileView
epi = Episystem(sppl, abenv, rel)
@profview simulate!(epi, burnin, timestep)
