# Include simulation functions
using Simulation
using Base.Test
using Distributions
using RCall
using Unitful.DefaultSymbols


numSpecies=4

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2], numSpecies))

# Set probabilities
birth = 0.6/month
death = 0.6/month
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

opts = repmat([5.0°C], numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies)  * °C
traits = TempTrait(opts, vars)
abun = Multinomial(individuals, numSpecies)
names = map(x -> "$x", 1:numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param)
abenv = simplehabitatAE(0.0, grid, totalK, area)
rel = TraitRelationship(GaussTemp)
eco = Ecosystem(sppl, abenv, rel)

@test_nowarn plot_move(eco, 2, 2, 1, false)

times = 10.0year; burnin = 1.0year; interval = 1.0month
# Run simulation grid
lensim = length(0.0month:interval:times)

# Run simulations 10 times
abun = generate_storage(eco, lensim, 1)
simulate!(eco, burnin, interval, timestep)
simulate_record!(abun, eco, times, interval, timestep)


@test_nowarn plot_abun(abun, numSpecies, grid)
