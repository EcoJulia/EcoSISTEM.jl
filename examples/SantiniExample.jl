## REPRODUCTION OF SANTINI PAPER ##


## STATIC VERSION
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using RCall
#using JLD
using StatsBase

# Start parallel processes
#addprocs(9)

# Load packages to all cores
using Diversity
using Simulation
using DataFrames
using Distributions
using DataStructures

# Set up Ecosystem as discrete habitat with species having a trait preference
# for one of the two niche types
numSpecies = 50
numTraits = 2
numNiches = 2
numInvasive = 1

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(sample(2.0:10, numSpecies+numInvasive))

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

#
abun = Multinomial(individuals, numSpecies)
native = fill(true, numSpecies + numInvasive)
native[numSpecies+numInvasive] = false
traits = sample([1,2,3], weights([0.4, 0.2, 0.4]), numSpecies)
push!(traits, 3)
sppl = SpeciesList(numSpecies + numInvasive, DiscreteTrait(traits), abun,
                   energy_vec, movement, param, native)
abenv = simplenicheAE(numNiches, grid, totalK, area)
rel = Match{eltype(abenv.habitat)}()
eco = Ecosystem(populate!, sppl, abenv, rel)


# Set up scenario of total habitat loss at certain rate
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
    burnin = 3year; interval = 3month; reps = 1000
    lensim = length(0month:interval:times)
    abun = SharedArray(zeros(1, length(divfuns), lensim, length(scenario), reps))
    #abun = zeros(1, length(divfuns), lensim, length(scenario), reps)
    @sync @parallel for i in 1:length(scenario)
        print(scenario)
    #for i in 1:length(scenario)
        for j in 1:reps
            reenergise!(eco, totalK, grid)
            repopulate!(eco);
            thisabun = view(abun, :, :, :, i, j);
            simulate!(eco, burnin, timestep)
            simulate_record_diversity!(thisabun, eco, times, interval, timestep,
            scenario[i], divfuns, q)
        end
    end
    abun
end

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 0year; interval = 3month; reps = 1
    lensim = length(0month:interval:times)
    #abun = SharedArray(zeros(1, length(divfuns), lensim, length(scenario), reps))
    abun = zeros(1, length(divfuns), lensim, length(scenario), reps)
    for i in 1:length(scenario)
        for j in 1:reps
            reenergise!(eco, totalK, grid);
            repopulate!(eco);
            thisabun = view(abun, :, :, :, i, j);
            #simulate!(eco, burnin, timestep)
            simulate_record_diversity!(thisabun, eco, times, interval, timestep,
            scenario[i], divfuns, q)
        end
    end
    abun
end

times = 1year;
#div = runsim(eco, times);
#save("SantiniRun50spp.jld", "div", div)
#div[isnan.(div)] = 0.0
div = load("SantiniRun50spp.jld", "div")

# Simple niche preference ecosystem with 150 species over 10,000 km2 area
#
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using RCall
using JLD

# Start parallel processes
addprocs(15)

# Load packages to all cores
@everywhere using Diversity
@everywhere using Simulation
@everywhere using DataFrames
@everywhere using Distributions
@everywhere using DataStructures

# Set up Ecosystem as discrete habitat with species having a trait preference
# for one of the two niche types
numSpecies = 50
numTraits = 2
numNiches = 2
numInvasive = 1

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(sample(2.0:10, numSpecies+numInvasive))

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
grid = (50, 50)
area = 10000.0km^2
totalK = 1000000.0 * (numSpecies + numInvasive)
individuals=20000 * numSpecies

# Create movement type - all individuals are allowed to move and have a wide range
kernel = GaussianKernel(8.0km, numSpecies+numInvasive, 10e-04)
movement = AlwaysMovement(kernel)

#
abun = Multinomial(individuals, numSpecies)
native = Vector{Bool}(numSpecies + numInvasive)
fill!(native, true)
native[numSpecies+numInvasive] = false
sppl = SpeciesList(numSpecies + numInvasive, numTraits, abun,
                   energy_vec, movement, UniqueTypes(numSpecies+numInvasive),
                   param, native)
abenv = simplenicheAE(numNiches, grid, totalK, area)
rel = Match{eltype(abenv.habitat)}()
eco = Ecosystem(sppl, abenv, rel)

# Set up scenario of total habitat loss at certain rate
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
    burnin = 3year; interval = 3month; reps = 1000
    lensim = length(0month:interval:times)
    abun = SharedArray(zeros(1, length(divfuns), lensim, length(scenario), reps))
    #abun = zeros(1, length(divfuns), lensim, length(scenario), reps)
    @sync @parallel for i in 1:length(scenario)
        print(scenario)
    #for i in 1:length(scenario)
        for j in 1:reps
            reenergise!(eco, totalK, grid)
            repopulate!(eco);
            thisabun = view(abun, :, :, :, i, j);
            simulate!(eco, burnin, timestep)
            simulate_record_diversity!(thisabun, eco, times, interval, timestep,
            scenario[i], divfuns, q)
        end
    end
    abun
end

times = 10year;
#div = runsim(eco, times);
#save("SantiniRun50spp.jld", "div", div)
#div[isnan.(div)] = 0.0
div = load("SantiniRun50spp.jld", "div")

function slopes(divarray::Array{Float64, 1})
    map(x -> (divarray[x] - divarray[x - 1])/divarray[x], 2:length(divarray))
end

function slopes(divarray::Array{Float64, 5})
    mapslices(slopes, divarray, 3)
end

slopemat = slopes(Array(div)) .* 100
meanslope = mapslices(x->mean(x[.!isnan.(x)]), slopemat, [1, 3, 5])[1,:,1,:, 1]
repslope = mapslices(x->mean(x[.!isnan.(x)]), slopemat, [1, 3])[1,:,1,:, :]
repslope = mapslices(x->all(x.>0) || all(x.<0), repslope, 3)[:,:,1]
repslope = reshape(repslope, 80, 1)
@rput meanslope
@rput repslope
R"
library(fields);
library(viridis);library(RColorBrewer)
png('Meanslope_static50.png', width = 1000, height = 800)
par(mfrow=c(1,1), mar=c(4,4,6,6));
image(meanslope, axes=FALSE, xlab='', ylab='', srt=45, col = magma(50),
    breaks =seq(-4, 4,length.out=51));
axis(1, at = seq(0,1, length.out=10),
labels = c('norm alpha q1', 'raw alpha q1', 'norm beta q1', 'raw beta q1', 'norm rho q1',
'raw rho q1', 'gamma q1', 'richness', 'shannon', 'simpson'));
axis(2, at = seq(0,1, length.out=8),
labels = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive','HabRep', 'HabLoss'));
image.plot(meanslope, col = magma(50), legend.only=TRUE,
    breaks =seq(-2, 2,length.out=51), legend.lab ='% change in diversity metric')
mat = expand.grid(seq(0,1, length.out=10), seq(0,1, length.out=8));
mat = mat[repslope, ]
points(mat[,1],mat[,2], pch=8, col ='white')
dev.off()
"



function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 3year; interval = 3month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, 1)
    simulate!(eco, burnin, timestep)
    simulate_record!(abun, eco, times, interval, timestep, scenario[1])
end

times = 1year;
abun = runsim(eco, times);
#cleanup!(abun)
plot_abun(abun, numSpecies, grid, 1)
plot_mean(abun, numSpecies, grid)
