using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Profile
using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using Diversity

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

# Create abiotic environment - even grid of one temperature
abenv = simplehabitatAE(274.0K, grid, totalK, area)


# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

# Create transition list
transitions = TransitionList(true)
addtransition!(transitions, UpdateEnergy(EcoSISTEM.update_energy_usage!))
addtransition!(transitions, UpdateEnvironment(update_environment!))
for spp in eachindex(sppl.species.names)
    for loc in eachindex(abenv.habitat.matrix)
        addtransition!(transitions, BirthProcess(spp, loc, sppl.params.birth[spp]))
        addtransition!(transitions, DeathProcess(spp, loc, sppl.params.death[spp]))
        addtransition!(transitions, AllDisperse(spp, loc))
    end
end
transitions = specialise_transition_list(transitions)

#Create ecosystem
eco = Ecosystem(sppl, abenv, rel, transitions = transitions)

# Simulation Parameters
burnin = 5years; times = 50years; timestep = 1month; record_interval = 3months; repeats = 1
lensim = length(0years:record_interval:times)
abuns = zeros(Int64, numSpecies, prod(grid), lensim)
# Burnin
@time simulate!(eco, burnin, timestep, specialise = true);

# Create transition list
transitions = TransitionList()
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
@time simulate!(eco, burnin, timestep);


transitions = TransitionList(true)
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
transitions = specialise_transition_list(transitions)

# Simulation Parameters
burnin = 5years; times = 50years; timestep = 1month; record_interval = 3months; repeats = 1
lensim = length(0years:record_interval:times)
abuns = zeros(Int64, numSpecies, prod(grid), lensim)
# Burnin
@time simulate!(eco, burnin, timestep);

transitions = TransitionList()
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
@time simulate!(eco, burnin, timestep, specialise = true);