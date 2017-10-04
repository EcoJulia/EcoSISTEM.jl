using Simulation
using Base.Test
using Distributions
using Unitful.DefaultSymbols
using RCall
using Diversity

numSpecies=4
numTraits = 2
numNiches = 2
# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2.0], numSpecies))

# Set probabilities
birth = 6.0/year
death = 6.0/year
l = 1.0
s = 0.5
boost = 1.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = Simulation.equalpop(EqualPop(birth, death, l, s, boost), numSpecies)

grid = (10,10)
area = 100.0km^2
totalK = 100000.0
individuals=10000

# Create ecosystem
kernel = GaussianKernel(0.1km, numSpecies, 10e-4)
movement = AlwaysMovement(kernel)
abun = Multinomial(individuals, numSpecies)

niches = DiscreteTrait(rand([1, 2], 4))
opts = repmat([5.0°C], numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies) * °C
temp = ContinuousTrait(opts, vars)

trtcol = TraitCollection2(niches, temp)
sim = UniqueTypes(numSpecies)
sppl = SpeciesList{typeof(trtcol), typeof(energy_vec),
            typeof(movement), typeof(sim), typeof(param)}(trtcol, rand(abun),
             energy_vec, sim, movement, param)

abenv1 = simplenicheAE(numNiches, grid, totalK, area)
abenv2 = tempgradAE(-10.0°C, 10.0°C, grid, totalK, area,
0.0°C/month)
hab = Simulation.HabitatCollection2(abenv1.habitat, abenv2.habitat)
abenv = GridAbioticEnv{typeof(hab), typeof(abenv1.budget)}(hab, abenv1.active,
abenv1.budget, abenv1.names)
rel1 = TraitRelationship{eltype(abenv.habitat)[1]}(SimpleNiche)
rel2 = TraitRelationship{eltype(abenv.habitat)[2]}(GaussTemp)
rel = Simulation.multiplicativeTR2(rel1, rel2)
eco = Ecosystem(sppl, abenv, rel)

#Simulation.update!(eco, 1month)
times = 10.0year; burnin = 1.0year; interval = 1.0month
# Run simulation grid
lensim = length(0.0month:interval:times)

# Run simulations 10 times
abun = generate_storage(eco, lensim, 1)
simulate!(eco, burnin, interval, timestep)
simulate_record!(abun, eco, times, interval, timestep)


plot_abun(abun, numSpecies, grid)
plot_mean(abun, numSpecies, grid)
