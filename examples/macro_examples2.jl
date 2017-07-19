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
## Run simulation on 2 by 1 grid - investigate spatial distributions
function runsim(times::Int64)
  # Set up initial parameters for ecosystem
  numSpecies = 50
  numTraits = 2
  numNiches = 2

  # Set up how much energy each species consumes
  energy_vec = SimpleRequirement(repmat([2.0], numSpecies))

  # Set probabilities
  birth = 0.6/12
  death = 0.6/12
  l = 1.0
  s = 0.5
  timestep = 1/12

  # Collect model parameters together (in this order!!)
  param = EqualPop(birth, death, l, s)

  minT = -10.0
  maxT = 10.0

  grid = (20, 20)
  gridSize = 0.1
  totalK = 100000.0
  individuals=10000

  # Create ecosystem
  kernel = GaussianKernel(0.2, numSpecies, 10e-4)
  movement = BirthOnlyMovement(kernel)

  #opts = repmat([5.0], numSpecies) #collect(linspace(minT, maxT, 8))
  #vars = collect(1.0:1.0:8.0)[randperm( 8)]
  temps = readtable(Pkg.dir("Simulation", "examples", "Tempdat.csv"))[1:50, :]
  prefs = [temps[:OrigValueStr] temps[:OrigValueStr]+1]
  conv = collect(linspace(-5, 20, 10))
  prefs = mapslices(x -> conv[x], prefs, 1)
  opts = vcat(mapslices(x -> rand(Uniform(x[1], x[2])), prefs, 2)...)
  #vars = collect(linspace(1, 4, numSpecies))
  vars = rand(Uniform(0, 25/9), numSpecies)
  traits = TempTrait(opts, vars)
  abun = Multinomial(individuals, numSpecies)
  names = map(x -> "$x", 1:numSpecies)
  sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
   movement, param)
  abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 0.0)
  rel = TraitRelationship(GaussTemp)
  eco = Ecosystem(sppl, abenv, rel)

  #plotdiv(norm_meta_alpha, eco, collect(0:10))

  steady_state = 100; burnin = 200; interval = 1; reps = 1;
  lensim = length(0:interval:times);
  lensteady = length(0:interval:steady_state);

  # Run simulations 10 times
  storage_steady = generate_storage(eco, lensteady, reps);
  storage_change = generate_storage(eco, lensim, reps);
  #storage_newsteady = generate_storage(eco, lensim, reps);
  for j in 1:reps
    if (j != 1) repopulate!(eco) end
    simulate!(eco, burnin, interval, timestep);
    thisstoresteady = view(storage_steady, :, :, :, j);
    simulate_record!(thisstoresteady, eco, steady_state, interval, timestep);
    eco.abenv.habitat.change.rate = 0.01
    thisstorechange = view(storage_change, :, :, :, j)
    simulate_record!(thisstorechange, eco, times, interval, timestep)
    #eco.abenv.habitat.change.rate = 0.0
    #thisstorenewsteady = view(storage_newsteady, :, :, :, j)
    #simulate_record!(thisstorenewsteady, eco, steady_state, interval, timestep)
  end
  save("macrorun_1steady.jld", "storage_change", storage_change,
  "storage_steady", storage_steady)
end

runsim(1200)


## New plots
storage_steady = load("macrorun_1steady.jld", "storage_steady")
storage_change = load("macrorun_1steady.jld", "storage_change")
storage = cat(3, storage_steady, storage_change)

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
jpeg(paste('plots/newrun/alphaq2_1steady', i, '.jpg'), quality=100)
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
jpeg(paste('plots/newrun/gammaq2_1steady', i, '.jpg'), quality=100)
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
    tempbetas[:, sub] = norm_meta_beta(met, 0:10)[:diversity]
  end
end
tempbetas = reshape(tempbetas, 11, 20, 20)
@rput tempbetas;
R"library(fields);par(mfrow=c(1,1))
jpeg('plots/newrun/tempbetaq1_1steady.jpg', quality=100)
im = tempbetas[ 1 , , ]
im[im==0] = 1.0
image.plot(1:20, 1:20, t(im),col=magma(20), breaks = seq(0.99,2,length.out=21),
xlab='', ylab='')
dev.off()
"

"""
HABITAT PLOT
"""
abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 0.0)
eco = Ecosystem(sppl, abenv, rel)
eco.abenv.habitat.change.rate = 0.0
for i in 1:1302
  if (i>=101) eco.abenv.habitat.change.rate = 0.01
  else eco.abenv.habitat.change.rate = 0.0 end
  temps = unique(eco.abenv.habitat.matrix)
  hab =eco.abenv.habitat.matrix
  @rput hab; @rput i
  if any(i.==divtimes)
    R"library(viridis);library(fields)
    jpeg(paste('plots/newrun/tempchange_1steady', i, '.jpg'))
    image.plot(1:20, 1:20, t(hab),col=magma(35), breaks = c(-10:25),
    xlab='', ylab='')
    dev.off()"
  end
  TempChange(eco)
end
