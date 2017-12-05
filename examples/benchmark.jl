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


using Diversity
using Simulation
using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols

numSpecies = 10
numTraits = 2
numNiches = 2
# Set up how much energy each species consumes
energy_vec = SimpleRequirement(sample(2.0:10, numSpecies))

# Set probabilities
birth = 0.6/month
death = 0.6/month
long = 1.0
surv = 0.0
boost = 1000.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, long, surv, boost)

# Set up a large grid of 4km^2 grid squares and enough energy to sustain
# 1 million individuals
grid = (10, 10)
area = 1000.0km^2
totalK = 1000000.0 * numSpecies
individuals=20000 * numSpecies

# Create movement type - all individuals are allowed to move and have a wide range
kernel = GaussianKernel(4.0km, numSpecies, 10e-04)
movement = AlwaysMovement(kernel)

#
abun = Multinomial(individuals, numSpecies)
native = Vector{Bool}(numSpecies)
fill!(native, true)
sppl = SpeciesList(numSpecies, numTraits, abun,
                   energy_vec, movement,UniqueTypes(numSpecies),
                   param, native)
abenv = simplenicheAE(numNiches, grid, totalK, area)
rel = Match{eltype(abenv.habitat)}()
eco = Ecosystem(sppl, abenv, rel)

habloss = 1.0 /10year
declines = 1.0 /10year
scenario = [SimpleScenario(UniformDecline, declines),
    SimpleScenario(ProportionalDecline, declines),
    SimpleScenario(LargeDecline, declines),
    SimpleScenario(RareDecline, declines),
    SimpleScenario(CommonDecline, declines),
    SimpleScenario(Invasive, declines),
    SimpleScenario(HabitatReplacement, habloss),
    SimpleScenario(RandHabitatLoss!, habloss)]
divfuns = [norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta,
    norm_meta_rho, raw_meta_rho, meta_gamma, meta_speciesrichness, meta_shannon, meta_simpson]
q = 1.0

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 3month; reps = 1
    lensim = length(0month:interval:times)
    abun = zeros(1, length(divfuns), lensim, length(scenario))
    #abun = zeros(1, length(divfuns), lensim, length(scenario), reps)
    for i in 1:length(scenario)
            reenergise!(eco, totalK, grid)
            repopulate!(eco);
            thisabun = view(abun, :, :, :, i);
            simulate!(eco, burnin, timestep)
            simulate_record_diversity!(thisabun, eco, times, interval, timestep,
            scenario[i], divfuns, q)
    end
    abun
end

times = 10year;
runsim(eco, 1year)

using BenchmarkTools
using ProfileView
Profile.clear()  # in case we have any previous profiling data
@profile runsim(eco, times)
#Profile.print(format := flat)
ProfileView.view()
