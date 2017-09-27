## Code for running model on spatial grid
#addprocs(5)
# Include simulation functions
using Diversity
using Simulation
using RCall
using Distributions
using StatsBase
using ProgressMeter
using DataStructures
#using ArrayViews
using DataFrames
using JLD
using BenchmarkTools
using Unitful.DefaultSymbols
## Run simulation on 2 by 1 grid - investigate spatial distributions
function runsim(times::Int64)
  # Set up initial parameters for ecosystem
  numSpecies = 50

  # Set up how much energy each species consumes
  energy_vec = SimpleRequirement(repmat([2.0], numSpecies))

  # Set probabilities
  birth = 0.6/year
  death = 0.6/year
  l = 1.0
  s = 0.5
  boost = 1.0
  timestep = 1month

  # Collect model parameters together (in this order!!)
  param = EqualPop(birth, death, l, s, boost)

  minT = 0.0°C
  maxT = 10.0°C

  grid = (20, 20)
  area = 4.0km^2
  totalK = 100000.0
  individuals=10000

  # Create ecosystem
  kernel = GaussianKernel(0.4km, numSpecies, 10e-4)
  movement = BirthOnlyMovement(kernel)

  #opts = repmat([5.0], numSpecies) #collect(linspace(minT, maxT, 8))
  #vars = collect(1.0:1.0:8.0)[randperm( 8)]
  temps = readtable(Pkg.dir("Simulation", "examples", "Tempdat.csv"))[1:50, :]
  prefs = [temps[:OrigValueStr] temps[:OrigValueStr]+1]
  conv = collect(linspace(-5, 20, 10))
  prefs = mapslices(x -> conv[x], prefs, 1)
  opts = vcat(mapslices(x -> rand(Uniform(x[1], x[2])), prefs, 2)...) * °C
  #vars = collect(linspace(1, 4, numSpecies))
  vars = rand(Uniform(0, 25/9), numSpecies) * °C
  traits = TempTrait(opts, vars)
  abun = Multinomial(individuals, numSpecies)
  names = map(x -> "$x", 1:numSpecies)
  sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
   movement, param)
  abenv = tempgradAE(minT, maxT, grid, totalK, area, 0.0°C/month)
  rel = TraitRelationship(GaussTemp)
  eco = Ecosystem(sppl, abenv, rel)

  #plotdiv(norm_meta_alpha, eco, collect(0:10))

  steady_state = 100month; burnin = 200month; interval = 1month; reps = 1;
  lensim = length(0month:interval:times);
  lensteady = length(0month:interval:steady_state);

  # Run simulations 10 times
  storage_steady = generate_storage(eco, lensteady, reps);
  storage_change = generate_storage(eco, lensim, reps);
  #storage_newsteady = generate_storage(eco, lensim, reps);
  for j in 1:reps
    if (j != 1) repopulate!(eco) end
    simulate!(eco, burnin, interval, timestep);
    thisstoresteady = view(storage_steady, :, :, :, j);
    simulate_record!(thisstoresteady, eco, steady_state, interval, timestep);
    resetrate!(eco, 0.01°C/month)
    thisstorechange = view(storage_change, :, :, :, j)
    simulate_record!(thisstorechange, eco, times, interval, timestep)
    #eco.abenv.habitat.change.rate = 0.0
    #thisstorenewsteady = view(storage_newsteady, :, :, :, j)
    #simulate_record!(thisstorenewsteady, eco, steady_state, interval, timestep)
  end
  #save("macrorun_smaller.jld", "storage_change", storage_change,
  #"storage_steady", storage_steady)
end

runsim(100year)

## New plots
storage_steady = load("examples/data/macrorun_smaller.jld", "storage_steady")
storage_change = load("examples/data/macrorun_smaller.jld", "storage_change")
storage = cat(3, storage_steady, storage_change)
store = reshape(storage, 50, 20, 20, 1302, 1)
im1 = mapslices(sum, store[:,:,:, 100, 1], [1, 3])
im2 = mapslices(sum, store[:,:,:, 1300, 1], [1, 3])
@rput im1
@rput im2

R"par(mfrow=c(2,1));library(plotrix);library(viridis)
plot(-10:9, im1, xlim=c(-15, 25), xaxs='i', yaxs='i',cex=0.5,
xlab ='', ylab='Total abundance', ylim=c(0, 550));
gradient.rect(-10, 0, 10, 10, col = magma(32)[1:20])
plot(2:21, im2, xlim=c(-15, 25), xaxs='i', yaxs='i',cex=0.5,
xlab ='Temperature', ylab='Total abundance', ylim=c(0, 550));
gradient.rect(2, 0, 22, 10, col = magma(32)[21:32])
"
"""
SPATIAL ALPHA PLOT
"""
divtimes = collect(1:10:1301)
alphas = zeros(400, 11, length(divtimes))
for i in 1:length(divtimes)
  met = Metacommunity(storage[:,:,divtimes[i],1], UniqueTypes(numSpecies))
  alphas[:,:,i] = norm_sub_alpha(met, 0:10)[:diversity]
end
alphas = reshape(alphas, 20, 20, 11, length(divtimes))
temps = unique(eco.abenv.habitat.matrix)
hab =eco.abenv.habitat.matrix
@rput hab

@rput alphas; @rput temps; @rput vars
for i in 1:length(divtimes)
  @rput i
R"library(fields);par(mfrow=c(1,1))
jpeg(paste('plots/newrun/alphaq2_smaller', i, '.jpg'), quality=100)
im = alphas[ , , 2, i]
im[is.na(im)] = 0
image.plot(1:20, 1:20, t(im),col=magma(30), breaks = c(0:30),
xlab='', ylab='')
dev.off()
"
end

"""
SPATIAL GAMMA PLOT
"""
gammas = zeros(400, 11, length(divtimes))
for i in 1:length(divtimes)
  met = Metacommunity(storage[:,:,divtimes[i],1], UniqueTypes(numSpecies))
  gammas[:,:,i] = sub_gamma(met, 0:10)[:diversity]
end
gammas = reshape(gammas, 20, 20, 11, length(divtimes))
temps = unique(eco.abenv.habitat.matrix)
@rput gammas; @rput temps; @rput vars
for i in 1:length(divtimes)
  @rput i
R"library(fields);par(mfrow=c(1,1))
jpeg(paste('plots/newrun/gammaq2_smaller', i, '.jpg'), quality=100)
im = gammas[ , , 2, i]
im[is.na(im)] = 0
image.plot(1:20, 1:20, t(im),col=magma(80), breaks = c(0:80),
xlab='', ylab='')
dev.off()
"
end

"""
TEMPORAL BETA PLOT
"""
tempbetas = zeros(11, 400)
for sub in 1:400
  if sum(storage[:, sub, divtimes[11:131], 1])== 0.0
    tempbetas[:, sub] = repmat([0.0], 11)
  else
    met = Metacommunity(storage[:, sub, divtimes[11:131], 1], UniqueTypes(numSpecies))
    tempbetas[:, sub] = raw_meta_beta(met, 0:10)[:diversity]
  end
end
tempbetas = reshape(tempbetas, 11, 20, 20)
@rput tempbetas;
R"library(fields);par(mfrow=c(1,1))
jpeg('plots/newrun/rawtempbetaq1_smaller.jpg', quality=100)
im = tempbetas[1 , , ]
im[im==0] = 1.0
image.plot(1:20, 1:20, t(im),col=magma(20), breaks = seq(0,0.05,length.out=21),
xlab='', ylab='')
dev.off()
"
"""
TEMPORAL GAMMA PLOT
"""
tempgammas = zeros(11, 400)
for sub in 1:400
  if sum(storage[:, sub, divtimes[11:131], 1])== 0.0
    tempgammas[:, sub] = repmat([0.0], 11)
  else
    met = Metacommunity(storage[:, sub, divtimes[11:131], 1], UniqueTypes(numSpecies))
    tempgammas[:, sub] = meta_gamma(met, 0:10)[:diversity]
  end
end
tempgammas = reshape(tempgammas, 11, 20, 20)
@rput tempgammas;
R"library(fields);par(mfrow=c(1,1))
jpeg('plots/newrun/tempgammaq1_smaller.jpg', quality=100)
im = tempgammas[1 , , ]
im[im==0] = 1.0
image.plot(1:20, 1:20, t(im),col=magma(35), breaks = c(0:35),
xlab='', ylab='')
dev.off()
"

"""
HABITAT PLOT
"""
eco = Ecosystem(sppl, abenv, rel)
eco.abenv.habitat.change.rate = 0.0
for i in 1:1203
  if (i>=101 && i<=1101) eco.abenv.habitat.change.rate = 0.01
  else eco.abenv.habitat.change.rate = 0.0 end
  temps = unique(eco.abenv.habitat.matrix)
  hab =eco.abenv.habitat.matrix
  @rput hab; @rput i
  if any(i.==divtimes)
    R"library(viridis);library(fields)
    jpeg(paste('plots/newrun/tempchange', i, '.jpg'))
    image.plot(1:20, 1:20, t(hab),col=viridis(35), breaks = c(-10:25),
    xlab='', ylab='')
    dev.off()"
  end
  TempChange(eco)
end
"""
SEPARATE INTO TEMPERATURE CLASSES
"""

cold = find(temps[:OrigValueStr] .<= 3)
warmer = find(temps[:OrigValueStr] .>= 4)

divtimes = collect(1:10:1301)
alphascold = zeros(400, 11, length(divtimes)); alphaswarm = zeros(400, 11, length(divtimes))
for i in 1:length(divtimes)
  metcold = Metacommunity(storage[cold, :, divtimes[i], 1], UniqueTypes(length(cold)))
  metwarm = Metacommunity(storage[warmer, :, divtimes[i], 1], UniqueTypes(length(warmer)))
  alphascold[:,:,i] = norm_sub_alpha(metcold, 0:10)[:diversity]
  alphaswarm[:,:,i] = norm_sub_alpha(metwarm, 0:10)[:diversity]
end
alphascold = reshape(alphascold, 20, 20, 11, length(divtimes))
alphaswarm = reshape(alphaswarm, 20, 20, 11, length(divtimes))

@rput alphaswarm; @rput alphascold
for i in 1:length(divtimes)
  @rput i
R"library(viridis);library(fields);
jpeg(paste('plots/newrun/alphasq1_warm_cold_smaller', i, '.jpg'), quality=100,
width = 1000)
par(mfrow=c(1,2), mar=c(4,3,4,5))
imc = alphascold[ , , 1, i]; imw = alphaswarm[ , , 1, i]
imc[is.na(imc)] = 0; imw[is.na(imw)] = 0
image.plot(1:20, 1:20, t(imc),col=magma(25), breaks = c(0:25),
xlab='', ylab='', main = list('Cold preference', cex=1.4));
image.plot(1:20, 1:20, t(imw),col=magma(25), breaks = c(0:25),
xlab='', ylab='', main = list('Warm preference', cex=1.4));
dev.off()
"
end
