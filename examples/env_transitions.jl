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
ages = 5
abun_v = DataFrame([
    (name="Environment", initial=0),
    (name="Force", initial=fill(0, ages)),
])
numvirus = sum(length.(abun_v.initial))

# Set initial population sizes for all human categories
susceptible = fill(round(Int64, 500_000 * prod(grid)/ages), ages)
exposed = fill(0, ages)
infected = fill(100 * prod(grid), ages)
abun_h = DataFrame([
    (name="Susceptible", type=Susceptible, initial=susceptible),
    (name="Exposed", type=OtherDiseaseState, initial=exposed),
    (name="Infected", type=Infectious, initial=infected),
    (name="Recovered", type=Removed, initial=fill(0, ages)),
    (name="Dead", type=Removed, initial=fill(0, ages)),
])
numclasses = nrow(abun_h) * ages

# Set non-pathogen mediated transitions
mu = fill(1/7days, ages)
sigma = fill(1/7days, ages)

# Set simulation parameters
birth = fill(0.0/day, numclasses, ages)
death = fill(0.0/day, numclasses, ages)
beta_force = fill(10.0/day, ages)
beta_env = fill(10.0/day, ages)
virus_growth = fill(0.1/day, ages)
virus_decay = 1.0/2days
age_mixing = fill(1.0, ages, ages)
param = (birth = birth, death = death, virus_growth = virus_growth,
        virus_decay = virus_decay, beta_env = beta_env,
        beta_force = beta_force, age_mixing = age_mixing)

cat_idx = reshape(1:(nrow(abun_h) * ages), ages, nrow(abun_h))
transitions = DataFrame([
  (from="Susceptible", from_id=cat_idx[:, 1], to="Exposed", to_id=cat_idx[:, 2], prob = (env = beta_env, force = beta_force)),
  (from="Exposed", from_id=cat_idx[:, 2], to="Infected", to_id=cat_idx[:, 3], prob=mu),
  (from="Infected", from_id=cat_idx[:, 3], to="Recovered", to_id=cat_idx[:,4], prob=sigma),
])

# Dispersal kernels for virus and disease classes
dispersal_dists = fill(10.0km, prod(grid))
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = EpiMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param, ages)

category_map = (
    "Susceptible" => cat_idx[:, 1],
    "Exposed" => cat_idx[:, 2],
    "Infected" => cat_idx[:, 3],
    "Recovered" => cat_idx[:, 4],
    "Dead" => cat_idx[:, 5],
)

extra_params = (env_exposure = 2.5,)
# Create epi system with all information
rel = Gauss{eltype(epienv.habitat)}()
transitions = Simulation.create_transition_list_env_SEIR(epilist, epienv, extra_params)
epi = EpiSystem(epilist, epienv, rel, transitions = transitions)

# Run simulation
times = 1month; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), floor(Int, times/timestep) + 1)
@time new_simulate_record!(abuns, epi, times, interval, timestep);
plot_epidynamics(epi, abuns, category_map = category_map)

transitions = Simulation.create_transition_list_SEIR(epilist, epienv)
epi = EpiSystem(epilist, epienv, rel, transitions = transitions)

# Run simulation
times = 1month; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), floor(Int, times/timestep) + 1)
@time new_simulate_record!(abuns, epi, times, interval, timestep);
plot_epidynamics(epi, abuns, category_map = category_map)
