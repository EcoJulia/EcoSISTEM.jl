using Simulation
using Distributions
#using RCall
using Unitful
using Unitful.DefaultSymbols

## TEST TEMPERATURE GRADIENT
numSpecies=4

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

grid = (5, 5)
area = 25.0km^2
totalK = 10000.0
individuals=1000

# Create ecosystem
kernel = GaussianKernel(2.0km, numSpecies, 10e-04)
movement = AlwaysMovement(kernel)

opts = repmat([5.0°C], numSpecies) #collect(linspace(minT, maxT, 8))
vars = rand(Uniform(0, 25/9), numSpecies) * °C
traits = ContinuousTrait(opts, vars)
abun = Multinomial(individuals, numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param)
loss = 0.0/month
abenv = degradedhabitatAE(0.0, grid, totalK, area, loss)
rel = TraitRelationship(GaussTemp)
eco = Ecosystem(sppl, abenv, rel)

#plot_move(eco, 4, 4, 1, true)

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 3month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, 1)
    simulate!(eco, burnin, interval, timestep)
    resetrate!(eco, 0.01/month)
    simulate_record!(abun, eco, times, interval, timestep)
end
times = 10year;
abun = runsim(eco, times);
