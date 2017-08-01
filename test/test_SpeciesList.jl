using Simulation
using Base.Test
using Distributions

## Run simulation over a grid and plot
numSpecies=4
numTraits = 2

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

individuals=100

# Create ecosystem
kernel = GaussianKernel(0.2, numSpecies, 10e-4)
movement = AlwaysMovement(kernel)

opts = repmat([5.0], numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies)
traits = TempTrait(opts, vars)
abun = Multinomial(individuals, numSpecies)
@test_nowarn sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param)
@test_nowarn sppl = SpeciesList(numSpecies, numTraits, Multinomial(individuals, numSpecies),
                   energy_vec, movement, param)
