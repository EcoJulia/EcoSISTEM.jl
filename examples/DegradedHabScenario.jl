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
abenv = tempgradAE(0.0°C, 10.0°C, grid, totalK, area, 0.0°C/month)
rel = TraitRelationship{eltype(abenv.habitat)}(GaussTemp)
eco = Ecosystem(sppl, abenv, rel)

# Set up scenario of total habitat loss at certain rate
loss = 5.0 /year
degradation = SimpleScenario(Simulation.ClustHabitatLoss!, loss)

#plot_move(eco, 4, 4, 1, true)

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 3month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, 1)
    simulate!(eco, burnin, interval, timestep)
    simulate_record!(abun, eco, times, interval, timestep, degradation)
end

times = 100year;
abun = runsim(eco, times);
plot_abun(abun, numSpecies, grid)