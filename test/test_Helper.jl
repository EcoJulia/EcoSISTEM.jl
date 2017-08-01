using Simulation
using Base.Test
using Distributions
using RCall

numSpecies=4
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2], numSpecies))

# Set probabilities
birth = 0.6
death = 0.6
l = 1.0
s = 0.0
boost = 1000.0
timestep = 1.0

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, l, s, boost)

grid = (5,5)
gridSize = 1.0
totalK = 1000.0
individuals=100

# Create ecosystem
kernel = GaussianKernel(0.2, numSpecies, 10e-4)
movement = AlwaysMovement(kernel)

opts = repmat([5.0], numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies)
traits = TempTrait(opts, vars)
abun = Multinomial(individuals, numSpecies)
names = map(x -> "$x", 1:numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param)
abenv = simplehabitatAE(0.0, grid, totalK, gridSize)
rel = TraitRelationship(GaussTemp)
eco = Ecosystem(sppl, abenv, rel)

times = 1000; burnin = 500; interval = 10
# Run simulation grid
lensim = length(0:interval:times)

# Run simulations 10 times
@test_nowarn generate_storage(eco, lensim, 1)
abun = generate_storage(eco, lensim, 1)
@test_nowarn simulate!(eco, burnin, interval, timestep)
@test_nowarn simulate_record!(abun, eco, times, interval, timestep)
