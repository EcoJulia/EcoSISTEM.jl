## Code for benchmarking model

# Include simulation functions and other modules
using Diversity
using Simulation
using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols

## Run simulation over a grid and plot
numSpecies=100

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2], numSpecies))

# Set probabilities
birth = 0.6/month
death = 0.6/month
long = 1.0
surv = 0.0
boost = 1000.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, long, surv, boost)

grid = (10, 10)
area = 10000.0km^2
totalK = 1000000.0
individuals=10000

# Create ecosystem
kernel = GaussianKernel(10.0km, numSpecies, 10e-4)
movement = AlwaysMovement(kernel)

opts = repmat([5.0°C], numSpecies) #collect(linspace(minT, maxT, 8))
vars = rand(Uniform(0, 25/9), numSpecies) * °C
traits = TempTrait(opts, vars)
abun = Multinomial(individuals, numSpecies)
names = map(x -> "$x", 1:numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param)
abenv = simplehabitatAE(0.0, grid, totalK, area)
rel = TraitRelationship(GaussTemp)
eco = Ecosystem(sppl, abenv, rel)

#plot_move(eco, 4, 4, 1, true)

function runsim(eco::Ecosystem, times::Unitful.Time)
burnin = 1year; interval = 3month
lensim = length(0month:interval:times)
abun = generate_storage(eco, lensim, 1)
simulate!(eco, burnin, interval, timestep)
simulate_record!(abun, eco, times, interval, timestep)
end

times = 10year;
runsim(eco, 1year)
#
using BenchmarkTools
using ProfileView
Profile.clear()  # in case we have any previous profiling data
@profile runsim(eco, times)
#Profile.print(format := flat)
ProfileView.view()
