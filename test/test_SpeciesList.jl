using Simulation
using Base.Test
using Unitful.DefaultSymbols
using Distributions

## Run simulation over a grid and plot
numSpecies=4
numTraits = 2

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

individuals=100

# Create ecosystem
kernel = GaussianKernel(2.0km, numSpecies, 10e-4)
movement = AlwaysMovement(kernel)

opts = repmat([5.0°C], numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies) * °C
traits = ContinuousTrait(opts, vars)
abun = Multinomial(individuals, numSpecies)
@test_nowarn sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param)
@test_nowarn sppl = SpeciesList(numSpecies, numTraits, abun,
                   energy_vec, movement, param)
