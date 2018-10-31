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
div = runsim(times)
save("SantiniRun_reeve150.jld", "div", div)

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
png('Meanslope_static150_reeve.png', width = 1000, height = 800)
par(mfrow=c(1,1), mar=c(6,4,4,10));
image(meanslope, axes=FALSE, xlab='', ylab='', srt=45, col = colorRampPalette(brewer.pal(11, 'RdBu'))(51),
   breaks =seq(-1, 1,length.out=52));
axis(1, at = seq(0,1, length.out=7), cex.axis=1.2,
labels = c('Raw alpha', 'Norm alpha', 'Raw beta',
'Norm beta','Raw rho', 'Norm rho',  'Gamma'));
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
png('Temp_trends_reeve.png', width = 1000, height = 900)
g = ggplot(data = standat, aes(x=time, y = div, group = rep))+
geom_line(col = alpha('black', 0.1))+facet_wrap(~measure)+ylim(-10, 15)+
xlab('Time (years)') + ylab('Diversity value')
print(g);dev.off()
"
sdmat = mapslices(std, slopemat, 3)[:, :, 1]
upper = meanslope .+ 1.96 .* (sdmat./sqrt(1000))
lower = meanslope .- 1.96 .* (sdmat./sqrt(1000))
@rput meanslope
@rput upper
@rput lower
R"library(ggplot2)
labels = c('Raw alpha', 'Norm alpha', 'Raw beta',
'Norm beta','Raw rho', 'Norm rho',  'Gamma')
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
    png('barplot_reeve.png', width = 1200, height = 800)
    print(g)
    dev.off()
"


function linmod(x)
    df = DataFrame(X = x, Y = 1:length(x))
    mod = GLM.lm(@formula(X ~ Y), df)
    return coef(mod)[2]
end
rawmat = mapslices(linmod, div[1, :, :, :, :], 2)[:, 1, :, :]
rawslope = mapslices(mean, rawmat, 3)[:, :, 1] * 10
sdmat = mapslices(std, rawmat, 3)[:, :, 1]
upper = rawslope .+ 1.96 .* (sdmat./sqrt(1000))
lower = rawslope .- 1.96 .* (sdmat./sqrt(1000))
@rput rawslope
@rput upper
@rput lower
R"library(ggplot2)
labels = c('Raw alpha', 'Norm alpha', 'Raw beta',
'Norm beta','Raw rho', 'Norm rho',  'Gamma')
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
    png('barplot_reeve_raw.png', width = 1200, height = 800)
    print(g)
    dev.off()
"
