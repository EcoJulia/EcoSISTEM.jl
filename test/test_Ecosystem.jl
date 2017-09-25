using Simulation
using Base.Test
using Unitful.DefaultSymbols
using Distributions

## Run simulation over a grid and plot
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

grid = (5, 5)
area = 25.0km^2
totalK = 1000.0
individuals=100

# Create ecosystem
kernel = GaussianKernel(2.0km, numSpecies, 10e-4)
movement = AlwaysMovement(kernel)

opts = repmat([5.0°C], numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies) * °C
traits = TempTrait(opts, vars)
abun = Multinomial(individuals, numSpecies)
names = map(x -> "$x", 1:numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param)
abenv = tempgradAE(-10.0°C, 10.0°C, grid, totalK, area,
0.01°C/month)
rel = TraitRelationship(GaussTemp)
@test_nowarn eco = Ecosystem(sppl, abenv, rel)
eco = Ecosystem(sppl, abenv, rel)
@test_nowarn gettraitrel(eco)
@test_nowarn gethabitat(eco)
@test_nowarn gethabitat(eco, 1)
@test_nowarn getenvtype(eco)
@test_nowarn getsize(eco)
@test_nowarn getgridsize(eco)
@test_nowarn getdispersaldist(eco)
@test_nowarn getdispersaldist(eco, 1)
@test_nowarn getdispersaldist(eco, "1")
@test_nowarn getdispersalvar(eco)
@test_nowarn getdispersalvar(eco, 1)
@test_nowarn getdispersalvar(eco, "1")
@test_nowarn resetrate!(eco, 0.5°C/month)
@test_nowarn resetrate!(eco, 0.5/month)
