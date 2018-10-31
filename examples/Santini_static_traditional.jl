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

val = 1.0
var = 0.5

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
    SimpleScenario(RandHabitatLoss!, habloss),
    SimpleScenario(ClustHabitatLoss!, habloss),
    SimpleScenario(SusceptibleDecline, declines)]
divfuns = [sorenson, geom_mean_abun]
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
            Simulation.resettraits!(sppl.types.tree)
            sppl.susceptible = ContinuousEvolve(val, var, sppl.types.tree).mean
            #Simulation.resettraits!(sppl.types.tree)
            #sppl.requirement.energy = ContinuousEvolve(5.0, 0.3, sppl.types.tree).mean
            #prob = 1./sppl.requirement.energy
            #sppl.abun = rand(Multinomial(individuals, prob./sum(prob)))
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
            Simulation.resettraits!(sppl.types.tree)
            sppl.susceptible = ContinuousEvolve(val, var, sppl.types.tree).mean
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
#save("SantiniRun_trad150.jld", "div", div)
div = div[:, [1, 2, 5, 6, 4, 3, 7], :, :, :]
div[isnan.(div)] = 0
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
par(mfrow=c(1,1), mar=c(6,4,4,10));
image(meanslope, axes=FALSE, xlab='', ylab='', srt=45, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51),
   breaks =seq(-1, 1,length.out=52));
axis(1, at = seq(0,1, length.out=7), cex.axis=1.2,
labels = c('Sorenson', 'Richness', 'Mean abun',
'Geometric mean','Simpson', 'Shannon',  'PD'));
axis(2, at = seq(0,1, length.out=7), cex.axis=1.2,
labels = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive', 'Habitat \n Loss'));
image.plot(meanslope, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51), legend.only=TRUE,
   breaks =seq(-1, 1,length.out=52),legend.args=list( text='% change in metric',
   cex=1.2, side=4,padj = 0.8, line=2), legend.mar = 7)
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
png('Temp_trends_trad.png', width = 1000, height = 900)
g = ggplot(data = standat, aes(x=time, y = div, group = rep))+
geom_line(col = alpha('black', 0.1))+facet_wrap(~measure)+ylim(-10, 15)+
xlab('Time (years)') + ylab('Diversity value')
print(g);dev.off()
"

standat = DataFrame(div = reshape(div[1,2,:,:,:], 1089000))
standat[:time] = (repmat(collect(1:121), 9000).-1)./12
#standat[:time] = (repmat(vcat(map(x -> repmat([x], 9), 1:121)...), 1000) .- 1)./12
standat[:scenario] = repmat(vcat(map(x -> repmat([x], 121), ["Uniform", "Proportional", "Large", "Rare",
    "Common", "Invasive", "Rand hab", "Clust hab", "Susceptible"])...), 1000)
standat[:rep] = vcat(map(x -> repmat([x], 1089), 1:1000)...)
@rput standat
R"library(ggplot2); library(cowplot); library(scales)
standat$scenario = factor(standat$scenario, levels = c('Uniform', 'Proportional',
'Large', 'Rare', 'Common', 'Invasive', 'Rand hab', 'Clust hab', 'Susceptible'))
png('SR_trends_trad.png', width = 1000, height = 900)
g = ggplot(data = standat, aes(x=time, y = div, group = rep))+
geom_line(col = alpha('black', 0.1))+facet_wrap(~scenario)+
xlab('Time (years)') + ylab('Scenario')
print(g);dev.off()
"
standat = DataFrame(div = reshape(div[1,3,:,:,:], 1089000))
standat[:time] = (repmat(collect(1:121), 9000).-1)./12
#standat[:time] = (repmat(vcat(map(x -> repmat([x], 9), 1:121)...), 1000) .- 1)./12
standat[:scenario] = repmat(vcat(map(x -> repmat([x], 121), ["Uniform", "Proportional", "Large", "Rare",
    "Common", "Invasive", "Rand hab", "Clust hab", "Susceptible"])...), 1000)
standat[:rep] = vcat(map(x -> repmat([x], 1089), 1:1000)...)
@rput standat
R"library(ggplot2); library(cowplot); library(scales)
standat$scenario = factor(standat$scenario, levels = c('Uniform', 'Proportional',
'Large', 'Rare', 'Common', 'Invasive', 'Rand hab', 'Clust hab', 'Susceptible'))
png('Ab_trends_trad.png', width = 1000, height = 900)
g = ggplot(data = standat, aes(x=time, y = div, group = rep))+
geom_line(col = alpha('black', 0.1))+facet_wrap(~scenario)+
xlab('Time (years)') + ylab('Scenario')
print(g);dev.off()
"

sdmat = mapslices(std, slopemat, 3)[:, :, 1]
upper = meanslope .+ 1.96 .* (sdmat./sqrt(1000))
lower = meanslope .- 1.96 .* (sdmat./sqrt(1000))
@rput meanslope
@rput upper
@rput lower
R"library(ggplot2)
labels = c('Sorenson', 'Richness', 'Mean abun',
'Geometric mean','Simpson', 'Shannon',  'PD')
scenarios = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive',
    'Habitat Loss')
dat = data.frame()
for (i in 1:7){
    dat = rbind(dat, data.frame(mn = meanslope[i,], up = upper[i,], lo = lower[i,],
    measure = rep(labels[i], 7), scenario = scenarios))
    }
    dat$scenario = factor(dat$scenario, levels = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive',
        'Habitat Loss'))
    print(dat)
    g = ggplot(dat, aes(y= mn, x = scenario, fill = measure)) + geom_bar(stat = 'identity') +
        facet_wrap(~ measure, nrow = 2) + geom_hline(yintercept = 0)+
        geom_errorbar(aes(ymin=lo, ymax=up),
                width=.2, position=position_dodge(.9)) + xlab('Scenario of change')+
                ylab('Mean slope')+theme(legend.position = 'none',
                axis.text.x = element_text(size=12, angle = 45, vjust = 0.9, hjust = 1))
    png('barplot_trad.png', width = 1200, height = 800)
    print(g)
    dev.off()
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
    meanslopeT[(2*i - 1), :] = meanslopeE[:, i]
    meanslopeT[(2*i), :] = meanslopeL[:, i]
end
meanslopeT *= 10
repslopeT = Array{Bool, 2}(49*2, 1)
for i in 1:49
    repslopeT[(2*i - 1), :] = repslopeE[:, i]
    repslopeT[(2*i), :] = repslopeL[:, i]
end
@rput meanslopeT
@rput repslopeT
R"
library(fields);
library(viridis);library(RColorBrewer)
png('EWS_static150_reeve.png', width = 1000, height = 1200)
par(mfrow=c(1,1), mar=c(6,6,4,8));
metriclabels = c('Sorenson', 'Richness', 'Mean abun',
'Geometric mean','Simpson', 'Shannon',  'PD');

image(t(meanslopeT), axes=FALSE, xlab='', ylab='', srt=45, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51),
    breaks =seq(-1, 1,length.out=52));
axis(1, at = seq(0,1, length.out=7), metriclabels, cex = 1.2);
labels = c('Uniform \n \n', 'Proportional \n \n', 'Largest \n \n', 'Rarest \n \n', 'Common \n \n','Invasive \n \n', 'Habitat Loss \n \n')
sublabels = c('\n EWS', '\nLWS')
#labels = paste(expand.grid(sublabels, labels)[,1],
#    expand.grid(sublabels, labels)[,2])
axis(2, at = seq(0,1, length.out=14), rep(sublabels, 7));
axis(2, at = seq(0.05, 1, 0.15), labels, tick = F, font =2, cex = 1.2);
#text(x = -0.5, y = seq(0.05, 1, 0.15), labels, srt = 90)
image.plot(meanslopeT, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51), legend.only=TRUE,
    breaks =seq(-1, 1,length.out=52), legend.args=list( text='% change in metric',
    cex=1.2, side=4,padj = 0.8, line=2), legend.mar = 7)
mat = expand.grid(seq(0,1, length.out=7), seq(0,1, length.out=14));
mat = mat[repslopeT, ]
points(mat[,1],mat[,2], pch=8, col ='grey20')
dev.off()
"

abunvec = Multinomial(individuals, probs)
native = fill(true, numSpecies + numInvasive)
native[numSpecies+numInvasive] = false
sppl = SpeciesList(numSpecies + numInvasive, 2, abunvec,
                   energy_vec, movement, param, native, [0.5, 0.5])
abenv = simplenicheAE(numNiches, grid, totalK, area)
rel = Match{eltype(abenv.habitat)}()
eco = Ecosystem(trait_populate!, sppl, abenv, rel)
mat = Array{Float64, 3}(size(eco.abenv.budget.matrix, 1),size(eco.abenv.budget.matrix, 2), 121)
for i in 1:121
    ClustHabitatLoss!(eco, timestep, habloss)
    mat[:,:,i] = eco.abenv.budget.matrix
end
Simulation.reenergise!(eco, totalK, size(eco.abenv.budget.matrix))
mat2 = Array{Float64, 3}(size(eco.abenv.budget.matrix, 1),size(eco.abenv.budget.matrix, 2), 121)
for i in 1:121
    RandHabitatLoss!(eco, timestep, habloss)
    mat2[:,:,i] = eco.abenv.budget.matrix
end
val = reshape(mat[:,:, 121], 2500)
val2  = reshape(mat2[:,:, 121], 2500)
@rput val
@rput val2
R"library(viridis); library(ggplot2); library(cowplot)
dat = expand.grid(1:50, 1:50)
dat$val1 = val
dat$val2 = val2
g1 = ggplot(dat, aes(x= Var1, y = Var2)) + geom_tile(aes(fill = factor(val1)))+
    xlab('') + ylab('')+
    scale_fill_manual(values = c('tan4', 'palegreen4'),
        name='Habitat \n status',
        labels=c('Removed', 'Remaining'))+
theme(axis.text.x = element_blank(),
      axis.text.y=element_blank())+
theme(legend.position = 'right',
      axis.ticks = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank())
g2 = ggplot(dat, aes(x= Var1, y = Var2)) + geom_tile(aes(fill = factor(val2)))+
xlab('') + ylab('')+
scale_fill_manual(values = c('palegreen4', 'tan4'),
    name='Habitat \n status',
    labels=c('Removed', 'Remaining'))+
theme(axis.text.x = element_blank(),
  axis.text.y=element_blank())+
theme(legend.position = 'none',
  axis.ticks = element_blank(),
  axis.line.x = element_blank(),
  axis.line.y = element_blank())
g = ggdraw() +
  draw_plot(g1, 0, 0, 0.525, 0.95) +
  draw_plot(g2, 0.525, 0, 0.475, 0.95) +
  draw_plot_label(c('A', 'B'), c(0.01, 0.51), c(0.98, 0.98), size = 18)
png('Cluster_mat.png', width = 1200, height = 800)
print(g)
dev.off()"

rawmat = mapslices(linmod, div[1, :, :, :, :], 2)[:, 1, :, :]
rawslope = mapslices(mean, rawmat, 3)[:, :, 1] * 10
sdmat = mapslices(std, rawmat, 3)[:, :, 1]
upper = rawslope .+ 1.96 .* (sdmat./sqrt(1000))
lower = rawslope .- 1.96 .* (sdmat./sqrt(1000))
@rput rawslope
@rput upper
@rput lower
R"library(ggplot2);library(cowplot)
labels = c('Sorenson', 'Richness', 'Mean abun',
'Geometric mean','Simpson', 'Shannon',  'PD')
scenarios = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive',
    'Habitat Loss')
dat = data.frame()
for (i in 1:7){
    dat = rbind(dat, data.frame(mn = rawslope[i,], up = upper[i,], lo = lower[i,],
    measure = rep(labels[i], 7), scenario = scenarios))
    }
    dat$scenario = factor(dat$scenario, levels = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive',
        'Habitat Loss'))
    print(dat)
    g = ggplot(dat, aes(y= mn, x = scenario, fill = measure)) +
    geom_bar(stat = 'identity') +
        facet_wrap(~ measure, nrow = 2,  scales = 'free_y') +
        geom_hline(yintercept = 0)+
        geom_errorbar(aes(ymin=lo, ymax=up),
                width=.2) +
                xlab('Scenario of change')+
                ylab('Mean slope')+theme(legend.position = 'none',
                axis.text.x = element_text(size=12, angle = 45,
                 vjust = 0.9, hjust = 1))
    png('barplot_trad_raw.png', width = 1200, height = 800)
    print(g)
    dev.off()
"
