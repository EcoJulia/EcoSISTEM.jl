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
    SimpleScenario(LargeDecline, declines),
    SimpleScenario(RareDecline, declines),
    SimpleScenario(CommonDecline, declines),
    SimpleScenario(Invasive, declines),
    SimpleScenario(Invasive, declines),
    SimpleScenario(RandHabitatLoss!, habloss),
    SimpleScenario(ClustHabitatLoss!, habloss),
    SimpleScenario(SusceptibleDecline, declines)]
divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson,
    mean_abun, geom_mean_abun, pd]
q = 1.0

function runsim(times::Unitful.Time)
    interval = 1month; reps = 1000
    lensim = length(0month:interval:times)
    abun = SharedArray(zeros(1, length(divfuns), lensim, length(scenario), reps))
    @sync @parallel  for j in 1:reps
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


times = 10year
div = runsim(times)
save("SantiniRun_trad150.jld", "div", div)

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
geom_line(col = alpha('black', 0.1))+facet_wrap(~measure)+ylim(-5, 5)+
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
standat = DataFrame(div = reshape(div[1,4,:,:,:], 1089000))
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
xlab('Time (years)') + ylab('Geometric mean abundance')
print(g);dev.off()
"

sdmat = mapslices(std, slopemat, 3)[:, :, 1]
upper = meanslope .+ sdmat
lower = meanslope .- sdmat
rep = Array{String, 2}(7, 10)
for i in 1:7
    for j in 1:10
    if sum(slopemat[i,j,:] .> 0)/1000 >= 0.95 || sum(slopemat[i,j,:] .< 0)/1000 >= 0.95
        rep[i, j ] = "*"
    else
        rep[i, j] = ""
    end
end
end
@rput rep
@rput meanslope
@rput upper
@rput lower
R"library(ggplot2); library(cowplot)
labels = c('Sorenson', 'Richness', 'Mean abun',
'Geometric mean','Simpson', 'Shannon',  'PD')
scenarios = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive',
    'Phylo invasive','Rand hab loss', 'Clust hab loss', 'Susceptible')
dat = data.frame()
adj = c()
for (i in 1:7){
    dat = rbind(dat, data.frame(mn = meanslope[i,], up = upper[i,], lo = lower[i,],
    measure = rep(labels[i], 10), scenario = scenarios, rep = rep[i,]))
    }
    for (j in 1:nrow(dat)){
        dat$adj[j] =  ifelse(dat$mn[j] > 0, dat$mn[j] + 0.02, dat$mn[j] - 0.02)
        }
    dat$scenario = factor(dat$scenario, levels = c('Uniform', 'Proportional',
    'Largest', 'Rarest', 'Common','Invasive','Phylo invasive',
        'Rand hab loss', 'Clust hab loss', 'Susceptible'))
    g = ggplot(dat, aes(y= mn, x = scenario, fill = measure)) + geom_bar(stat = 'identity') +
        facet_wrap(~ measure, nrow = 2) + geom_hline(yintercept = 0)+
        geom_errorbar(aes(ymin=lo, ymax=up),
                width=.2, position=position_dodge(.9)) + xlab('Scenario of change')+
                ylab('Mean slope')+theme(legend.position = 'none',
                axis.text.x = element_text(size=12, angle = 45, vjust = 0.9, hjust = 1))+
                 geom_text(aes(y = adj, label=rep), size = 6)+
        scale_fill_manual(palette = colorRampPalette(brewer.pal(7, 'Paired')))
    png('barplot_trad.png', width = 1200, height = 800)
    print(g)
    dev.off()
"

# Check for early and late warning signals
early = stanmat[:,:, 1:24,:,:]

slopematE = mapslices(linmod, early[1, :, :, :, :], 2)[:, 1, :, :]
meanslopeE = mapslices(mean, slopematE, 3)[:, :, 1]
sdmatE = mapslices(std, slopematE, 3)[:, :, 1]
upperE = meanslopeE .+ sdmatE
lowerE = meanslopeE .- sdmatE
repE = Array{String, 2}(7, 10)
for i in 1:7
    for j in 1:10
    if sum(slopematE[i,j,:] .> 0)/1000 >= 0.95 || sum(slopematE[i,j,:] .< 0)/1000 >= 0.95
        repE[i, j ] = "*"
    else
        repE[i, j] = ""
    end
end
end
@rput repE
@rput meanslopeE
@rput upperE
@rput lowerE
R"library(ggplot2); library(cowplot)
labels = c('Sorenson', 'Richness', 'Mean abun',
'Geometric mean','Simpson', 'Shannon',  'PD')
scenarios = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive',
    'Phylo invasive','Rand hab loss', 'Clust hab loss', 'Susceptible')
dat = data.frame()
adj = c()
for (i in 1:7){
    dat = rbind(dat, data.frame(mn = meanslopeE[i,], up = upperE[i,], lo = lowerE[i,],
    measure = rep(labels[i], 10), scenario = scenarios, rep = repE[i,]))
    }
    for (j in 1:nrow(dat)){
        dat$adj[j] =  ifelse(dat$mn[j] > 0, dat$up[j] + 0.02, dat$lo[j] - 0.02)
        }
    dat$scenario = factor(dat$scenario, levels = c('Uniform', 'Proportional',
    'Largest', 'Rarest', 'Common','Invasive','Phylo invasive',
        'Rand hab loss', 'Clust hab loss', 'Susceptible'))
    g = ggplot(dat, aes(y= mn, x = scenario, fill = measure)) + geom_bar(stat = 'identity') +
        facet_wrap(~ measure, nrow = 2) + geom_hline(yintercept = 0)+
        geom_errorbar(aes(ymin=lo, ymax=up),
                width=.2, position=position_dodge(.9)) + xlab('Scenario of change')+
                ylab('Mean slope')+theme(legend.position = 'none',
                axis.text.x = element_text(size=12, angle = 45, vjust = 0.9, hjust = 1))+
                 geom_text(aes(y = adj, label=rep), size = 6)+
        scale_fill_manual(palette = colorRampPalette(brewer.pal(7, 'Paired')))
    png('Early_barplot_trad.png', width = 1200, height = 800)
    print(g)
    dev.off()
"

rawmat = mapslices(linmod, div[1, :, :, :, :], 2)[:, 1, :, :]
rawslope = mapslices(mean, rawmat, 3)[:, :, 1] * 10
sdmat = mapslices(std, rawmat, 3)[:, :, 1]
upper = rawslope .+ sdmat
lower = rawslope .- sdmat
maxrep = mapslices(maximum, upper, 2)
minrep = mapslices(minimum, lower, 2)
rng = abs(maxrep - minrep)
@rput rawslope
@rput upper
@rput lower
@rput rng
R"library(ggplot2); library(cowplot)
labels = c('Sorenson', 'Richness', 'Mean abun',
'Geometric mean','Simpson', 'Shannon',  'PD')
scenarios = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive',
    'Phylo invasive','Rand hab loss', 'Clust hab loss', 'Susceptible')
dat = data.frame()
adj = c()
for (i in 1:7){
    dat = rbind(dat, data.frame(mn = rawslope[i,], up = upper[i,], lo = lower[i,],
    measure = rep(labels[i], 10), scenario = scenarios, rep = rep[i,],
    rn = rep(rng[i], 10)))
    }
    for (j in 1:nrow(dat)){
        dat$adj[j] =  ifelse(dat$mn[j] > 0, dat$up[j] + 0.05 * dat$rn[j],
        dat$lo[j] - 0.05 * dat$rn[j])
        }
    dat$scenario = factor(dat$scenario, levels = c('Uniform', 'Proportional',
    'Largest', 'Rarest', 'Common','Invasive','Phylo invasive',
        'Rand hab loss', 'Clust hab loss', 'Susceptible'))
    g = ggplot(dat, aes(y= mn, x = scenario, fill = measure)) + geom_bar(stat = 'identity') +
        facet_wrap(~ measure, nrow = 2, scales = 'free_y') + geom_hline(yintercept = 0)+
        geom_errorbar(aes(ymin=lo, ymax=up),
                width=.2, position=position_dodge(.9)) + xlab('Scenario of change')+
                ylab('Mean slope')+theme(legend.position = 'none',
                axis.text.x = element_text(size=12, angle = 45, vjust = 0.9, hjust = 1))+
                 geom_text(aes(y = adj, label=rep), size = 6)+
        scale_fill_manual(palette = colorRampPalette(brewer.pal(7, 'Paired')))
    png('Raw_barplot_trad.png', width = 1200, height = 800)
    print(g)
    dev.off()
"
# Check for early and late warning signals
early = div[:,:, 1:24,:,:]

slopematE = mapslices(linmod, early[1, :, :, :, :], 2)[:, 1, :, :]
meanslopeE = mapslices(mean, slopematE, 3)[:, :, 1]
sdmatE = mapslices(std, slopematE, 3)[:, :, 1]
upperE = meanslopeE .+ sdmatE
lowerE = meanslopeE .- sdmatE
repE = Array{String, 2}(7, 10)
for i in 1:7
    for j in 1:10
    if sum(slopematE[i,j,:] .> 0)/1000 >= 0.95 || sum(slopematE[i,j,:] .< 0)/1000 >= 0.95
        repE[i, j ] = "*"
    else
        repE[i, j] = ""
    end
end
end
maxrepE = mapslices(maximum, upperE, 2)
minrepE = mapslices(minimum, lowerE, 2)
rngE = abs(maxrepE - minrepE)
@rput repE
@rput meanslopeE
@rput upperE
@rput lowerE
@rput rngE
R"library(ggplot2); library(cowplot)
labels = c('Sorenson', 'Richness', 'Mean abun',
'Geometric mean','Simpson', 'Shannon',  'PD')
scenarios = c('Uniform', 'Proportional', 'Largest', 'Rarest', 'Common','Invasive',
    'Phylo invasive','Rand hab loss', 'Clust hab loss', 'Susceptible')
dat = data.frame()
adj = c()
for (i in 1:7){
    dat = rbind(dat, data.frame(mn = meanslopeE[i,], up = upperE[i,], lo = lowerE[i,],
    measure = rep(labels[i], 10), scenario = scenarios, rep = repE[i,],
    rn = rep(rngE[i], 10)))
    }
    for (j in 1:nrow(dat)){
        dat$adj[j] =  ifelse(dat$mn[j] > 0, dat$up[j] + 0.05 * dat$rn[j],
        dat$lo[j] - 0.05 * dat$rn[j])
        }
    dat$scenario = factor(dat$scenario, levels = c('Uniform', 'Proportional',
    'Largest', 'Rarest', 'Common','Invasive','Phylo invasive',
        'Rand hab loss', 'Clust hab loss', 'Susceptible'))
    g = ggplot(dat, aes(y= mn, x = scenario, fill = measure)) + geom_bar(stat = 'identity') +
        facet_wrap(~ measure, nrow = 2, scales='free_y') + geom_hline(yintercept = 0)+
        geom_errorbar(aes(ymin=lo, ymax=up),
                width=.2, position=position_dodge(.9)) + xlab('Scenario of change')+
                ylab('Mean slope')+theme(legend.position = 'none',
                axis.text.x = element_text(size=12, angle = 45, vjust = 0.9, hjust = 1))+
                 geom_text(aes(y = adj, label=rep), size = 6)+
        scale_fill_manual(palette = colorRampPalette(brewer.pal(7, 'Paired')))
    png('Early_raw_barplot_trad.png', width = 1200, height = 800)
    print(g)
    dev.off()
"
