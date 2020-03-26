using Simulation
using Compat.Test
using Unitful.DefaultSymbols
using Distributions
using Simulation.Units

## Run simulation over a grid and plot
numSpecies=4
numTraits = 2

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(fill(2, numSpecies))

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

opts = fill(5.0°C, numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies) * °C
traits = GaussTrait(opts, vars)
abun = rand(Multinomial(individuals, numSpecies))
native = fill(true, numSpecies)
@test_nowarn sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
@test_nowarn sppl = SpeciesList(numSpecies, numTraits, abun,
                   energy_vec, movement, param, native)
