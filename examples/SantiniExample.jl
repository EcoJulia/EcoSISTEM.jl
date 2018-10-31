## REPRODUCTION OF SANTINI PAPER ##


## STATIC VERSION
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using RCall
using JLD
using StatsBase
using GLM

# Start parallel processes
addprocs(9)

# Load packages to all cores
@everywhere using Diversity
@everywhere using Simulation
@everywhere using DataFrames
@everywhere using Distributions
@everywhere using DataStructures

# Set up Ecosystem as discrete habitat with species having a trait preference
# for one of the two niche types
numSpecies = 150
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
probs = rand(LogNormal(1.0, 0.5), numSpecies)
probs /= sum(probs)

# Create movement type - all individuals are allowed to move and have a wide range
kernel = GaussianKernel(0.0km, numSpecies+numInvasive, 10e-04)
movement = NoMovement(kernel)

# Set up scenario of total habitat loss at certain rate
habloss = 1.0 /10year
declines = 1.0 /10year
scenario = [SimpleScenario(UniformDecline, declines),
    SimpleScenario(ProportionalDecline, declines),
    SimpleScenario(LargeDecline, declines),
    SimpleScenario(RareDecline, declines),
    SimpleScenario(CommonDecline, declines),
    SimpleScenario(Invasive, declines),
    SimpleScenario(RandHabitatLoss!, habloss)]
divfuns = [norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta,
    norm_meta_rho, raw_meta_rho, meta_gamma]
q = 1.0


function runsim(times::Unitful.Time)
    burnin = 3year; interval = 1month; reps = 1000
    lensim = length(0month:interval:times)
    abun = SharedArray(zeros(1, length(divfuns), lensim, length(scenario), reps))
    @sync @parallel  for j in 1:reps
        for i in 1:length(scenario)
            abunvec = Multinomial(individuals, probs)
            native = fill(true, numSpecies + numInvasive)
            native[numSpecies+numInvasive] = false
            sppl = SpeciesList(numSpecies + numInvasive, 2, abunvec,
                               energy_vec, movement, param, native, [0.5, 0.5])
            abenv = simplenicheAE(numNiches, grid, totalK, area)
            rel = Match{eltype(abenv.habitat)}()
            eco = Ecosystem(trait_populate!, sppl, abenv, rel)
            thisabun = view(abun, :, :, :, i, j);
            simulate!(eco, burnin, timestep)
            simulate_record_diversity!(thisabun, eco, times, interval, timestep,
            scenario[i], divfuns, q)
        end
    end
    abun
end


times = 10year;
div = runsim(times);
save("SantiniRun150spp.jld", "div", div)
#div[isnan.(div)] = 0.0
div = load("SantiniRun50spp.jld", "div")
div = Array(div)
div[isnan.(div)] = 0
standardise(x) = (x .- mean(x))./std(x)
#stanmat = mapslices(standardise, Array(div), (4,5))[1, :, :, :, :]
#stanmat = mapslices(standardise, div, 2)
stanmat = copy(div)
for divfun in 1:7
    stanmat[:, divfun, :, :, :] = standardise(div[:, divfun, :, :, :])
end

function linmod(x)
    df = DataFrame(X = x, Y = 1:length(x))
    mod = GLM.lm(@formula(X ~ Y), df)
    return coef(mod)[2]
end
slopemat = mapslices(linmod, stanmat[1, :, :, :, :], 2)[:, 1, :, :]
meanslope = mapslices(mean, slopemat, 3)[:, :, 1] *10
repslope = mapslices(x-> (sum(x .> 0)/length(x)) > 0.95 || (sum(x .< 0)/length(x)) > 0.95, slopemat, 3)[:, :, 1]
repslope = reshape(repslope, 49, 1)

#slopemat = slopes(Array(div)) .* 100
#meanslope = mapslices(x->mean(x[.!isnan.(x)]), slopemat, [1, 3, 5])[1,:,1,:, 1]
#repslope = mapslices(x->mean(x[.!isnan.(x)]), slopemat, [1, 3])[1,:,1,:, :]
#repslope = mapslices(x->all(x.>0) || all(x.<0), repslope, 3)[:,:,1]
#repslope = reshape(repslope, 70, 1)
@rput meanslope
@rput repslope
R"
library(fields);
library(viridis);library(RColorBrewer)
png('Meanslope_static150_reeve.png', width = 1000, height = 800)
par(mfrow=c(1,1), mar=c(4,4,6,8));
image(meanslope, axes=FALSE, xlab='', ylab='', srt=45, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51),
    breaks =seq(-1, 1,length.out=52));
axis(1, at = seq(0,1, length.out=7),
labels = c('Norm alpha q1', 'Raw alpha q1', 'Norm beta q1', 'Raw beta q1', 'Norm rho q1',
'Raw rho q1', 'Gamma q1'));
axis(2, at = seq(0,1, length.out=7),
labels = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive', 'Habitat \n Loss'));
image.plot(meanslope, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51), legend.only=TRUE,
    breaks =seq(-1, 1,length.out=52), legend.lab ='% change in diversity metric')
mat = expand.grid(seq(0,1, length.out=7), seq(0,1, length.out=7));
mat = mat[repslope, ]
points(mat[,1],mat[,2], pch=8, col ='grey20')
dev.off()
"
## STATIC VERSION
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using RCall
using JLD
using StatsBase
using GLM

# Start parallel processes
addprocs(20)

# Load packages to all cores
@everywhere using Diversity
@everywhere using Simulation
@everywhere using DataFrames
@everywhere using Distributions
@everywhere using DataStructures

# Set up Ecosystem as discrete habitat with species having a trait preference
# for one of the two niche types
numSpecies = 150
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
probs = rand(LogNormal(1.0, 0.5), numSpecies)
probs /= sum(probs)

# Create movement type - all individuals are allowed to move and have a wide range
kernel = GaussianKernel(0.0km, numSpecies+numInvasive, 10e-04)
movement = NoMovement(kernel)

# Set up scenario of total habitat loss at certain rate
habloss = 1.0 /10year
declines = 1.0 /10year
scenario = [SimpleScenario(UniformDecline, declines),
    SimpleScenario(ProportionalDecline, declines),
    SimpleScenario(LargeDecline, declines),
    SimpleScenario(RareDecline, declines),
    SimpleScenario(CommonDecline, declines),
    SimpleScenario(Invasive, declines),
    SimpleScenario(RandHabitatLoss!, habloss)]
divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson,
    mean_abun, geom_mean_abun, pd]
q = 1.0

function runsim(times::Unitful.Time)
    burnin = 3year; interval = 1month; reps = 1000
    lensim = length(0month:interval:times)
    abun = SharedArray(zeros(1, length(divfuns), lensim, length(scenario), reps))
    @sync @parallel  for j in 1:reps
        for i in 1:length(scenario)
            abunvec = Multinomial(individuals, probs)
            native = fill(true, numSpecies + numInvasive)
            native[numSpecies+numInvasive] = false
            sppl = SpeciesList(numSpecies + numInvasive, 2, abunvec,
                               energy_vec, movement, param, native, [0.5, 0.5])
            abenv = simplenicheAE(numNiches, grid, totalK, area)
            rel = Match{eltype(abenv.habitat)}()
            eco = Ecosystem(trait_populate!, sppl, abenv, rel)
            thisabun = view(abun, :, :, :, i, j);
            simulate!(eco, burnin, timestep)
            simulate_record_diversity!(thisabun, eco, times, interval, timestep,
            scenario[i], divfuns, q)
        end
    end
    abun
end

function runsim(times::Unitful.Time)
    burnin = 3year; interval = 1month; reps = 1
    lensim = length(0month:interval:times)
    #abun = SharedArray(zeros(1, length(divfuns), lensim, length(scenario), reps))
    abun = zeros(1, length(divfuns), lensim, length(scenario), reps)
    for j in 1:reps
        for i in 1:length(scenario)
            print(scenario[i].fun, "\n")
            # Ecosystem set up
            abunvec = Multinomial(individuals, probs)
            native = fill(true, numSpecies + numInvasive)
            native[numSpecies+numInvasive] = false
            sppl = SpeciesList(numSpecies + numInvasive, 2, abunvec,
                               energy_vec, movement, param, native, [0.5, 0.5])
            abenv = simplenicheAE(numNiches, grid, totalK, area)
            rel = Match{eltype(abenv.habitat)}()
            eco = Ecosystem(trait_populate!, sppl, abenv, rel)

            # Run simulation
            thisabun = view(abun, :, :, :, i, j);
            simulate!(eco, burnin, timestep)
            simulate_record_diversity!(thisabun, eco, times, interval, timestep,
            scenario[i], divfuns, q)
        end
    end
    abun
end

times = 10year
div = runsim(times)
save("SantiniRun_trad150.jld", "div", div)
div = div[:, [1, 2, 5, 6, 4, 3, 7], :, :, :]
standardise(x) = (x .- mean(x))./std(x)
stanmat = copy(div)
for divfun in 1:7
    stanmat[:, divfun, :, :, :] = standardise(div[:, divfun, :, :, :])
end
function linmod(x)
    df = DataFrame(X = x, Y = 1:length(x))
    mod = GLM.lm(@formula(X ~ Y), df)
    return coef(mod)[2]
end
slopemat = mapslices(linmod, stanmat[1, :, :, :, :], 2)[:, 1, :, :]
meanslope = mapslices(mean, slopemat, 3)[:, :, 1] * 10
repslope = mapslices(x-> (sum(x .> 0)/length(x)) > 0.95 || (sum(x .< 0)/length(x)) > 0.95, slopemat, 3)[:, :, 1]
repslope = reshape(repslope, 49, 1)
@rput meanslope
@rput repslope
R"
library(fields);
library(viridis);library(RColorBrewer)
png('Meanslope_static150_trad.png', width = 1000, height = 800)
par(mfrow=c(1,1), mar=c(6,4,4,8));
image(meanslope, axes=FALSE, xlab='', ylab='', srt=45, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51),
    breaks =seq(-1, 1,length.out=52));
axis(1, at = seq(0,1, length.out=7),
labels = c('Sorenson', 'Richness', 'Mean abun',
'Geometric mean','Simpson', 'Shannon',  'PD'));
axis(2, at = seq(0,1, length.out=7),
labels = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive', 'Habitat \n Loss'));
image.plot(meanslope, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51), legend.only=TRUE,
    breaks =seq(-1, 1,length.out=52), legend.lab ='% change in diversity metric')
mat = expand.grid(seq(0,1, length.out=7), seq(0,1, length.out=7));
mat = mat[repslope, ]
points(mat[,1],mat[,2], pch=8, col ='grey20')
dev.off()
"
# Plot trends over time
@rput stanmat
standat = DataFrame(div = reshape(stanmat[1,:,:,1,:], 847000))
standat[:measure] = repmat(repmat(["Sorenson", "Richness", "Mean abun",
"Geometric mean","Simpson", "Shannon",  "PD"], 121), 1000)
standat[:time] = (repmat(vcat(map(x -> repmat([x], 7), 1:121)...), 1000) .- 1)./12
standat[:rep] = vcat(map(x -> repmat([x], 847), 1:1000)...)

@rput standat
R"library(ggplot2); library(cowplot); library(scales)
png('Temp_trends.png', width = 1000, height = 900)
g = ggplot(data = standat, aes(x=time, y = div, group = rep))+
geom_line(col = alpha('black', 0.1))+facet_wrap(~measure)+ylim(-7, 10)+
xlab('Time (years)') + ylab('Diversity value')
print(g);dev.off()
"

# Check for early and late warning signals
early = stanmat[:,:, 1:24,:,:]
late = stanmat[:,:, 25:end,:,:]

slopematE = mapslices(linmod, early[1, :, :, :, :], 2)[:, 1, :, :]
meanslopeE = mapslices(mean, slopematE, 3)[:, :, 1]
repslopeE = mapslices(x-> (sum(x .> 0)/length(x)) > 0.95 || (sum(x .< 0)/length(x)) > 0.95, slopematE, 3)[:, :, 1]
repslopeE = reshape(repslopeE, 49, 1)

slopematL = mapslices(linmod, late[1, :, :, :, :], 2)[:, 1, :, :]
meanslopeL = mapslices(mean, slopematL, 3)[:, :, 1]
repslopeL = mapslices(x-> (sum(x .> 0)/length(x)) > 0.95 || (sum(x .< 0)/length(x)) > 0.95, slopematL, 3)[:, :, 1]
repslopeL = reshape(repslopeL, 49, 1)

meanslopeT = Array{Float64, 2}(14, 7)
for i in 1:7
    meanslopeT[(2*i - 1), :] = meanslopeE[i, :]
    meanslopeT[(2*i), :] = meanslopeL[i, :]
end
meanslopeT *= 10
repslopeT = Array{Bool, 2}(49*2, 1)
for i in 1:49
    repslopeT[(2*i - 1), :] = repslopeE[i, :]
    repslopeT[(2*i), :] = repslopeL[i, :]
end
@rput meanslopeT
@rput repslopeT
R"
library(fields);
library(viridis);library(RColorBrewer)
png('EWS_static150_trad.png', width = 1000, height = 1200)
par(mfrow=c(1,1), mar=c(6,6,4,8));
metriclabels = c('Sorenson', 'Richness', 'Mean abun',
'Geometric mean','Simpson', 'Shannon',  'PD');

image(t(meanslopeT), axes=FALSE, xlab='', ylab='', srt=45, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51),
    breaks =seq(-1, 1,length.out=52));
axis(1, at = seq(0,1, length.out=7), metriclabels);
labels = c('Uniform \n \n', 'Proportional \n \n', 'Largest \n \n', 'Rarest \n \n', 'Common \n \n','Invasive \n \n', 'Habitat Loss \n \n')
sublabels = c('\n EWS', '\nLWS')
#labels = paste(expand.grid(sublabels, labels)[,1],
#    expand.grid(sublabels, labels)[,2])
axis(2, at = seq(0,1, length.out=14), rep(sublabels, 7));
axis(2, at = seq(0.05, 1, 0.15), labels, tick = F, font =2);
#text(x = -0.5, y = seq(0.05, 1, 0.15), labels, srt = 90)
image.plot(meanslopeT, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51), legend.only=TRUE,
    breaks =seq(-1, 1,length.out=52), legend.lab ='% change in diversity metric')
mat = expand.grid(seq(0,1, length.out=7), seq(0,1, length.out=14));
mat = mat[repslopeT, ]
points(mat[,1],mat[,2], pch=8, col ='grey20')
dev.off()
"


# Simple niche preference ecosystem with 150 species over 10,000 km2 area
#
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using RCall
using JLD

# Start parallel processes
addprocs(9)

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

div = Array(div)
standardise(x) = (x .- mean(x))./std(x)
#stanmat = mapslices(standardise, Array(div), (4,5))[1, :, :, :, :]
#stanmat = mapslices(standardise, div, 2)
stanmat = copy(div)
for divfun in 1:10
    stanmat[:, divfun, :, :, :] = standardise(div[:, divfun, :, :, :])
end
function linmod(x)
    df = DataFrame(X = x, Y = 1:length(x))
    mod = GLM.lm(@formula(X ~ Y), df)
    return coef(mod)[2]
end
slopemat = mapslices(linmod, stanmat[1, :, :, :, :], 2)[:, 1, :, :]
meanslope = mapslices(mean, slopemat, 3)[:, :, 1] .* 10
repslope = mapslices(x-> (sum(x .> 0)/length(x)) > 0.95 || (sum(x .< 0)/length(x)) > 0.95, slopemat, 3)[:, :, 1]
repslope = reshape(repslope, 80, 1)
meanslope = meanslope[:, end:-1:1]

#slopemat = slopes(Array(div)) .* 100
#meanslope = mapslices(x->mean(x[.!isnan.(x)]), slopemat, [1, 3, 5])[1,:,1,:, 1]
#repslope = mapslices(x->mean(x[.!isnan.(x)]), slopemat, [1, 3])[1,:,1,:, :]
#repslope = mapslices(x->all(x.>0) || all(x.<0), repslope, 3)[:,:,1]
#repslope = reshape(repslope, 70, 1)
@rput meanslope
@rput repslope
R"
library(fields);
library(viridis);library(RColorBrewer)
png('Meanslope_dynamic50_norm.png', width = 1000, height = 800)
par(mfrow=c(1,1), mar=c(4,4,6,8));
image(meanslope, axes=FALSE, xlab='', ylab='', srt=45, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51),
    breaks =seq(-2, 2,length.out=52));
axis(1, at = seq(0,1, length.out=10),
labels = c('Norm alpha q1', 'Raw alpha q1', 'Norm beta q1', 'Raw beta q1', 'Norm rho q1',
'Raw rho q1', 'Gamma q1', 'Richness', 'Shannon', 'Simpson'));
axis(2, at = seq(0,1, length.out=8),
labels = rev(c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive', 'Habitat \n Loss',
'Habitat \n rep')));
image.plot(meanslope, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51), legend.only=TRUE,
    breaks =seq(-2, 2,length.out=52), legend.lab ='% change in diversity metric')
mat = expand.grid(seq(0,1, length.out=10), seq(0,1, length.out=8));
mat = mat[repslope, ]
points(mat[,1],mat[,2], pch=8, col ='grey20')
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
