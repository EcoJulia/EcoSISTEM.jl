using Simulation
using Base.Test
using Distributions
using Unitful.DefaultSymbols

@test_nowarn rel = TraitRelationship{Unitful.Temperature}(GaussTemp)
@test_nowarn rel = TraitRelationship{Int64}(SimpleNiche)


numSpecies=4
numTraits = 2
numNiches = 2
# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2], numSpecies))

# Set probabilities
birth = 6.0/year
death = 6.0/year
l = 1.0
s = 0.0
boost = 1000.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, l, s, boost)

grid = (5,5)
area = 10.0km^2
totalK = 1000.0
individuals=100

# Create ecosystem
kernel = GaussianKernel(2.0km, numSpecies, 10e-4)
movement = AlwaysMovement(kernel)
abun = Multinomial(individuals, numSpecies)

# Test out discrete trait update
sppl = SpeciesList(numSpecies, numTraits, Multinomial(individuals, numSpecies),
                   energy_vec, movement, param)
abenv = simplenicheAE(numNiches, grid, totalK, area)
rel = TraitRelationship{eltype(abenv.habitat)}(SimpleNiche)
eco = Ecosystem(sppl, abenv, rel)
@test_nowarn Simulation.update!(eco, 1.0year)


numSpecies=4
# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2], numSpecies))

# Set probabilities
birth = 6.0/year
death = 6.0/year
l = 1.0
s = 0.0
boost = 1000.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, l, s, boost)

grid = (5,5)
area = 10.0km^2
totalK = 1000.0
individuals=100

# Create ecosystem
kernel = GaussianKernel(2.0km, numSpecies, 10e-4)
movement = AlwaysMovement(kernel)
abun = Multinomial(individuals, numSpecies)
opts = repmat([5.0°C], numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies) * °C
traits = ContinuousTrait(opts, vars)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param)
abenv = tempgradAE(-10.0°C, 10.0°C, grid, totalK, area, 0.01°C/month)
rel = TraitRelationship{eltype(abenv.habitat)}(GaussTemp)
eco = Ecosystem(sppl, abenv, rel)
@test_nowarn Simulation.update!(eco, 1.0month)
