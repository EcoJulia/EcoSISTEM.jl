## STATIC VERSION
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using RCall
using JLD
using StatsBase
using GLM
using Diversity
using Simulation
using DataFrames
using Distributions
using DataStructures

# Set up Ecosystem as discrete habitat with species having a trait preference
# for one of the two niche types
numSpecies = 150
numTraits = 2
numNiches = 2
numInvasive = 1

sus_mean = 1.0
sus_var = 0.5

size_mean = 5.0
size_var = 15.0

# Set probabilities
birth = fill(0.0/month, numSpecies+numInvasive)
death = fill(0.0/month, numSpecies+numInvasive)
long = 1.0
surv = 0.0
boost = 1000.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = NoGrowth{typeof(unit(0.0/month))}(birth, death, long, surv, boost)

# Set up a large grid of 4km^2 grid squares and enough energy to sustain
# 1 million individuals
grid = (50, 50)
area = 10000.0km^2
totalK = 1000000.0 * (numSpecies + numInvasive)
individuals=20000 * numSpecies

# Create movement type - all individuals are allowed to move and have a wide range
kernel = GaussianKernel(0.0km, numSpecies+numInvasive, 10e-04)
movement = NoMovement(kernel)

# Set up scenario of total habitat loss at certain rate
habloss = 1.0 /10year
declines = 1.0 /10year
scenario = [SimpleScenario(UniformDecline, declines),
    SimpleScenario(ProportionalDecline, declines),
    SimpleScenario(RareDecline, declines),
    SimpleScenario(CommonDecline, declines),
    SimpleScenario(LargeDecline, declines),
    SimpleScenario(Invasive, declines),
    SimpleScenario(Invasive, declines),
    SimpleScenario(SusceptibleDecline, declines),
    SimpleScenario(RandHabitatLoss!, habloss),
    SimpleScenario(ClustHabitatLoss!, habloss)]
divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson,
    mean_abun, geom_mean_abun, pd]
q = 1.0

native = fill(true, numSpecies + numInvasive)
native[numSpecies+numInvasive] = false
pop_mass = rand(Normal(-0.75, 0.1))
sppl = SpeciesList(numSpecies + numInvasive, 2, pop_mass, size_mean,
size_var, area, movement, param, native, [0.5, 0.5])
Simulation.resettraits!(sppl.types.tree)
sppl.susceptible = ContinuousEvolve(sus_mean, sus_var, sppl.types.tree).mean
abenv = simplenicheAE(numNiches, grid, totalK, area)
rel = Match{eltype(abenv.habitat)}()
eco = Ecosystem(trait_populate!, sppl, abenv, rel)

function runsim(times::Unitful.Time)
    interval = 1month; reps = 1
    lensim = length(0month:interval:times)
    abun = zeros(1, length(divfuns), lensim, length(scenario), reps)
    for j in 1:reps
        for i in 1:length(scenario)
            native = fill(true, numSpecies + numInvasive)
            native[numSpecies+numInvasive] = false
            pop_mass = rand(Normal(-0.75, 0.1))
            sppl = SpeciesList(numSpecies + numInvasive, 2, pop_mass, size_mean,
            size_var, area, movement, param, native, [0.5, 0.5])
            Simulation.resettraits!(sppl.types.tree)
            sppl.susceptible = ContinuousEvolve(sus_mean, sus_var, sppl.types.tree).mean
            if i == 7
                reroot!(sppl.types.tree, "151")
            end
            abenv = simplenicheAE(numNiches, grid, totalK, area)
            rel = Match{eltype(abenv.habitat)}()
            eco = Ecosystem(trait_populate!, sppl, abenv, rel)
            thisabun = view(abun, :, :, :, i, j);
            simulate_record_diversity!(thisabun, eco, times, interval, timestep,
            scenario[i], divfuns, q)
        end
    end
    abun
end


times = 1year
div = runsim(times)
