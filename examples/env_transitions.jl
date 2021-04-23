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
        beta_force = beta_force, age_mixing = age_mixing,
        env_exposure = 2.0/mean(epienv.habitat.matrix))

cat_idx = reshape(1:(nrow(abun_h) * ages), ages, nrow(abun_h))
transitiondat = DataFrame([
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
epilist = SpeciesList(traits, abun_v, abun_h, movement, transitiondat, param, ages)

category_map = (
    "Susceptible" => cat_idx[:, 1],
    "Exposed" => cat_idx[:, 2],
    "Infected" => cat_idx[:, 3],
    "Recovered" => cat_idx[:, 4],
    "Dead" => cat_idx[:, 5],
)

rel = Gauss{eltype(epienv.habitat)}()


transitions = create_transition_list()
addtransition!(transitions, UpdateEpiEnvironment(update_epi_environment!))
for loc in eachindex(epienv.habitat.matrix)
    addtransition!(transitions, ViralLoad(loc, param.virus_decay))
    for age in 1:ages
        addtransition!(transitions, ForceProduce(cat_idx[age, 3], loc, param.virus_growth[age]))
        addtransition!(transitions, ForceDisperse(cat_idx[age, 3], loc))
        addtransition!(transitions, EnvTransition(Exposure(transitiondat[1, :from_id][age], loc,
            transitiondat[1, :to_id][age], transitiondat[1, :prob].force[age], transitiondat[1, :prob].env[age]),
            param.env_exposure))
        addtransition!(transitions, Infection(transitiondat[2, :from_id][age], loc,
            transitiondat[2, :to_id][age], transitiondat[2, :prob][age]))
        addtransition!(transitions, Recovery(transitiondat[3, :from_id][age], loc,
            transitiondat[3, :to_id][age], transitiondat[3, :prob][age]))
    end
end

# Create epi system with all information
epi = Ecosystem(epilist, epienv, rel, transitions = transitions)

# Run simulation
times = 1month; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), floor(Int, times/timestep) + 1)
@time simulate_record!(abuns, epi, times, interval, timestep);
plot_epidynamics(epi, abuns, category_map = category_map)

# Compare to non environmental exposure
transitions = create_transition_list()#
addtransition!(transitions, UpdateEpiEnvironment(update_epi_environment!))
for loc in eachindex(epienv.habitat.matrix)
    addtransition!(transitions, ViralLoad(loc, param.virus_decay))
    for age in 1:ages
        addtransition!(transitions, ForceProduce(cat_idx[age, 3], loc, param.virus_growth[age]))
        addtransition!(transitions, ForceDisperse(cat_idx[age, 3], loc))
        addtransition!(transitions, Exposure(transitiondat[1, :from_id][age], loc,
            transitiondat[1, :to_id][age], transitiondat[1, :prob].force[age], transitiondat[1, :prob].env[age]))
        addtransition!(transitions, Infection(transitiondat[2, :from_id][age], loc,
            transitiondat[2, :to_id][age], transitiondat[2, :prob][age]))
        addtransition!(transitions, Recovery(transitiondat[3, :from_id][age], loc,
            transitiondat[3, :to_id][age], transitiondat[3, :prob][age]))
    end
end
epi = Ecosystem(epilist, epienv, rel, transitions = transitions)

# Run simulation
times = 1month; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), floor(Int, times/timestep) + 1)
@time simulate_record!(abuns, epi, times, interval, timestep);
plot_epidynamics(epi, abuns, category_map = category_map)
