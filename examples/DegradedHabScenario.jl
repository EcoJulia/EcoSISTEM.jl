using Simulation
using Distributions
#using RCall
using Unitful
using Unitful.DefaultSymbols
using Diversity
using RCall

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

grid = (10, 10)
area = 25.0km^2
totalK = 10000.0
individuals=1000

# Create ecosystem
kernel = GaussianKernel(0.5km, numSpecies, 10e-04)
movement = AlwaysMovement(kernel)

opts = repmat([5.0°C], numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies) * °C
traits = ContinuousTrait(opts, vars)
abun = Multinomial(individuals, numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
movement, param)
abenv = tempgradAE(0.0°C, 10.0°C,  grid, totalK, area, 0.0°C/month)
rel = Gauss{eltype(abenv.habitat)}()
eco = Ecosystem(sppl, abenv, rel)

# Set up scenario of total habitat loss at certain rate
loss = 0.01 /year
degradation = SimpleScenario(Simulation.ClustHabitatLoss!, loss)

# Set up scenario of total habitat loss at certain rate
loss = 0.01/year
level = 0.0
recovery = 0.1/month
lag = 0.0month
degradation = DisturbanceScenario(HabitatDisturbance!, loss, level,
recovery, lag)


function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 3month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, 1)
    simulate!(eco, burnin, timestep)
    simulate_record!(abun, eco, times, interval, timestep, degradation)
end

times = 10year;
abun = runsim(eco, times);
#cleanup!(abun)
plot_abun(abun, numSpecies, grid, 1)
budget = eco.abenv.budget.matrix
@rput budget
R"par(mfrow=c(1,1));library(fields);image.plot(t(budget))"
plot_mean(abun, numSpecies, grid)

# -------------------------------------------------------------------------#
#                       DIVERSITY SIMULATION                               #
# -------------------------------------------------------------------------#
abenv = simplehabitatAE(5.0, grid, totalK, area)
eco = Ecosystem(sppl, abenv, rel)
eco = repopulate!(eco)

# Set up scenario of total habitat loss at certain rate
loss = 0.01 /year
degradation = SimpleScenario(Simulation.HabitatDisturbance!, loss)
divfuns = [norm_sub_alpha, raw_sub_alpha, norm_sub_beta, raw_sub_beta,
norm_sub_rho,raw_sub_rho, sub_gamma]
qs = 1.0

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 3month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, length(divfuns), lensim, 10)
    simulate!(eco, burnin, timestep)
    simulate_record_diversity!(abun, eco, times, interval, timestep, degradation,
    divfuns, qs)
end

times = 10year;
diversity = runsim(eco, times);
#divclean = cleanup!(diversity)
plot_diversity(diversity, 1, (5,5), 1)

function slopes(divarray::Array{Float64, 1})
    map(x -> (divarray[x] - divarray[x - 1]), 2:length(divarray))
end

function slopes(divarray::Array{Float64, 4})
    mapslices(slopes, divarray, 3)
end
slopemat = slopes(diversity)
function sign_change(a::BitArray{1})
    sum(map(x-> if (a[x]!=a[x-1]) 1 else 0 end, 2:5))
end
function sign_change(slopemat::Array{Float64, 4})
    sign  = slopemat .>= 0
    mapslices(sign_change, sign, 3)
end
signmat = sign_change(slopemat)[:,:,1,1]
meanslope = mapslices(mean, slopemat, [1, 3, 4])[:,:,1,1]
@rput meanslope
R"library(fields);par(mfrow=c(1,1));image.plot(meanslope)"
meanslope = mapslices(mean, slopemat, [3, 4])[:,:,1,1]
meanslope = reshape(meanslope, grid[1], grid[2], length(divfuns))
names = ["Norm alpha","Raw alpha", "Norm beta", "Raw beta", "Norm rho","Raw rho",
 "Gamma"]
 abunmat = mapslices(sum, eco.abundances.grid, 1)[1,:,:]
@rput meanslope
@rput names
@rput abunmat
budget = eco.abenv.budget.matrix
@rput budget
R"par(mfrow=c(3,3));library(fields); library(viridis)
for (i in c(1:7)){
    image.plot(t(meanslope[,,i]), col=magma(9),
    main = names[i])
    }
image.plot(t(budget), main = 'Energy')
image.plot(t(abunmat), main = 'Abundances')
    "
#, breaks= seq(min(meanslope, na.rm=T),max(meanslope, na.rm = T), length.out=10)
