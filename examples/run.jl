# Include simulation functions
using Simulation
using Distributions
#using RCall
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
kernel = GaussianKernel(10.0km, numSpecies, 10e-04)
movement = AlwaysMovement(kernel)

opts = repmat([5.0°C], numSpecies) #collect(linspace(minT, maxT, 8))
vars = rand(Uniform(0, 25/9), numSpecies) * °C
traits = TempTrait(opts, vars)
abun = Multinomial(individuals, numSpecies)
names = map(x -> "$x", 1:numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param)
active = Array{Bool, 2}(grid)
fill!(active, true)
active[1, 1] = false
abenv = simplehabitatAE(0.0, grid, totalK, area, active)
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

#Profile.clear_malloc_data()
#runsim(eco, times)
@time runsim(eco, times)

abun = runsim(eco, times);
plot_abun(abun, numSpecies, 1, collect(1:16))
