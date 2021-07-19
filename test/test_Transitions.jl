using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using DataFrames
using Test

## Test generic biodiversity
# Set up initial parameters for ecosystem
numSpecies = 10; grid = (10, 10); req= 10.0kJ; individuals=10_000; area = 1000.0*km^2; totalK = 1000.0kJ/km^2

# Set up how much energy each species consumes
energy_vec = SolarRequirement(fill(req, numSpecies))


# Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.0
boost = 1000.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
movement = AlwaysMovement(kernel, Torus())


# Create species list, including their temperature preferences, seed abundance and native status
opts = fill(274.0K, numSpecies)
vars = fill(0.5K, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
sppl.params.birth

# Create abiotic environment - even grid of one temperature
abenv = simplehabitatAE(274.0K, grid, totalK, area)


# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

# Create transition list
transitions = create_transition_list()
addtransition!(transitions, UpdateEnergy(EcoSISTEM.update_energy_usage!))
addtransition!(transitions, UpdateEnvironment(update_environment!))
for spp in eachindex(sppl.species.names)
    for loc in eachindex(abenv.habitat.matrix)
        addtransition!(transitions, BirthProcess(spp, loc, sppl.params.birth[spp]))
        addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
        addtransition!(transitions, AllDisperse(spp, loc))
    end
end

#Create ecosystem
eco = Ecosystem(sppl, abenv, rel, transitions = transitions)


# Simulation Parameters
burnin = 5years; times = 50years; timestep = 1month; record_interval = 3months; repeats = 1
lensim = length(0years:record_interval:times)
abuns = zeros(Int64, numSpecies, prod(grid), lensim)
# Burnin
@test_nowarn simulate!(eco, burnin, timestep);
@test_nowarn simulate_record!(abuns, eco, times, record_interval, timestep);

## Test plant biodiversity
# Change movement to seed dispersal
movement = BirthOnlyMovement(kernel, Torus())
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)

# Create transition list
transitions = create_transition_list()
addtransition!(transitions, UpdateEnergy(EcoSISTEM.update_energy_usage!))
addtransition!(transitions, UpdateEnvironment(update_environment!))
for spp in eachindex(sppl.species.names)
    for loc in eachindex(abenv.habitat.matrix)
        addtransition!(transitions, GenerateSeed(spp, loc, sppl.params.birth[spp]))
        addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
        addtransition!(transitions, SeedDisperse(spp, loc))
    end
end

#Create ecosystem
eco = Ecosystem(sppl, abenv, rel, transitions = transitions)

# Simulation Parameters
burnin = 5years; times = 50years; timestep = 1month; record_interval = 3months; repeats = 1
lensim = length(0years:record_interval:times)
abuns = zeros(Int64, numSpecies, prod(grid), lensim)
# Burnin
@test_nowarn simulate!(eco, burnin, timestep);
@test_nowarn simulate_record!(abuns, eco, times, record_interval, timestep);



## Test SEIR
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

transitiondat = DataFrame([
  (from="Susceptible", from_id=1, to="Exposed", to_id=2, prob = (env = beta_env, force = beta_force)),
  (from="Exposed", from_id=2, to="Infected", to_id=3, prob=mu),
  (from="Infected", from_id=3, to="Recovered", to_id=4, prob=sigma),
])

# Dispersal kernels for virus and disease classes
dispersal_dists = fill(10.0km, prod(grid))
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = EpiMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = SpeciesList(traits, abun_v, abun_h, movement, transitiondat, param)
rel = Gauss{eltype(epienv.habitat)}()

# Create list of transitions for the simulation
transitions = create_transition_list()
addtransition!(transitions, SeedInfection(seedinfected!))
addtransition!(transitions, UpdateEpiEnvironment(update_epi_environment!))
for loc in eachindex(epienv.habitat.matrix)
    addtransition!(transitions, ForceProduce(3, loc, param.virus_growth))
    addtransition!(transitions, ViralLoad(loc, param.virus_decay))
    addtransition!(transitions, Exposure(transitiondat[1, :from_id], loc,
        transitiondat[1, :to_id], transitiondat[1, :prob].force, transitiondat[1, :prob].env))
    addtransition!(transitions, Infection(transitiondat[2, :from_id], loc,
        transitiondat[2, :to_id], transitiondat[2, :prob]))
    addtransition!(transitions, Recovery(transitiondat[3, :from_id], loc,
        transitiondat[3, :to_id], transitiondat[3, :prob]))
    for spp in eachindex(epilist.species.names)
        addtransition!(transitions, ForceDisperse(spp, loc))
    end
end

# Create epi system with all information
epi = Ecosystem(epilist, epienv, rel, transitions = transitions)

# Run simulation
times = 1month; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), floor(Int, times/timestep) + 1)
@test_nowarn simulate!(epi, times, timestep);
@test_nowarn simulate_record!(abuns, epi, times, interval, timestep);


## Test SEIR with environmental exposure
# Create list of transitions for the simulation
transitions = create_transition_list()
addtransition!(transitions, SeedInfection(seedinfected!))
addtransition!(transitions, UpdateEpiEnvironment(update_epi_environment!))
for loc in eachindex(epienv.habitat.matrix)
    addtransition!(transitions, ForceProduce(3, loc, param.virus_growth))
    addtransition!(transitions, ViralLoad(loc, param.virus_decay))
    addtransition!(transitions, Exposure(transitiondat[1, :from_id], loc,
        transitiondat[1, :to_id], transitiondat[1, :prob].force, transitiondat[1, :prob].env))
    addtransition!(transitions, Infection(transitiondat[2, :from_id], loc,
        transitiondat[2, :to_id], transitiondat[2, :prob]))
    addtransition!(transitions, Recovery(transitiondat[3, :from_id], loc,
        transitiondat[3, :to_id], transitiondat[3, :prob]))
    for spp in eachindex(epilist.species.names)
        addtransition!(transitions, ForceDisperse(spp, loc))
    end
end

import EcoSISTEM.getprob
getprob(eco::Ecosystem, rule::Exposure) = (rule.force_prob * 2.0/K * get_env(eco.abenv.habitat, rule.location), rule.virus_prob)

# Create epi system with all information
epi = Ecosystem(epilist, epienv, rel, transitions = transitions)

# Run simulation
times = 1month; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), floor(Int, times/timestep) + 1)
@test_nowarn simulate!(epi, times, timestep);
@test_nowarn simulate_record!(abuns, epi, times, interval, timestep);
