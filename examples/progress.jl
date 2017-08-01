using RCall
using StatsBase
using ProgressMeter
#using ArrayViews
using JLD

addprocs(5)
@everywhere using Diversity
@everywhere using Simulation
@everywhere using DataFrames
@everywhere using Distributions
@everywhere using DataStructures

"""
ONE GRID SQUARE, VARYING SIZES
"""

@everywhere function runsim(times::Int64, reps::Int64, gridSize::Float64)
  # Set up initial parameters for ecosystem
  numSpecies = 25

  # Set up how much energy each species consumes
  energy_vec = SimpleRequirement(repmat([2.0], numSpecies))

  # Set probabilities
  birth = 0.6/12
  death = 0.6/12
  l = 1.0
  s = 0.0
  boost = 1000.0
  timestep = 1/12

  # Collect model parameters together (in this order!!)
  param = EqualPop(birth, death, l, s, boost)

  minT = 0.0
  maxT = 0.0

  grid = (1, 1)
  totalK = 1000.0 * (gridSize^2)
  individuals = 100

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
  #abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 0.0)
  abenv = simplehabitatAE(0.0, grid, totalK, gridSize)
  rel = TraitRelationship(GaussTemp)
  eco = Ecosystem(sppl, abenv, rel)

  burnin = 1000; interval = 1;
  lensim = length(0:interval:times);

  # Run simulations 10 times
  storage = SharedArray(generate_storage(eco, lensim, reps));
  #storage_newsteady = generate_storage(eco, lensim, reps);
  @sync @parallel for j in 1:reps
    if (j != 1) repopulate!(eco) end
    simulate!(eco, burnin, interval, timestep);
    thisstore = view(storage, :, :, :, j);
    simulate_record!(thisstore, eco, times, interval, timestep); end
    return storage
end
storage = runsim(200, 1, 1.0)
@rput storage
R"for (i in c(1:25)){
if (i==1) plot_fun = plot else plot_fun =lines
plot_fun(1:501, storage[i,1,,1], col=i, type='l',
ylim=c(0,25))}"



"""
Number of species surviving at equilibrium
"""
reps = 100
times = 500
sizes = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
meanSR =
  @sync @parallel (vcat) for k in sizes
    storage = runsim(times, reps, k)
    sum(storage[:, :, (times+1), 1:reps] .> 0)/reps
  end
meanSR
@rput sizes; @rput meanSR
R"plot(log(sizes), meanSR, type = 'l', xlab='Log size',
ylab='Average species richness', ylim=c(0,25),
xaxs='i',yaxs='i')"

"""
Species richness vs metabolism
"""

@everywhere function runsim(times::Int64, reps::Int64, maxE::Float64)
  # Set up initial parameters for ecosystem
  numSpecies = 25

  # Set up how much energy each species consumes
  energy_vec = SimpleRequirement(repmat([maxE], numSpecies))
#collect(linspace(2, maxE, 25))
  # Set probabilities
  birth = 0.6/12
  death = 0.6/12
  l = 1.0
  s = 0.0
  boost = 1000.0
  timestep = 1.0

  # Collect model parameters together (in this order!!)
  param = EqualPop(birth, death, l, s, boost)

  minT = 0.0
  maxT = 0.0

  grid = (1, 1)
  gridSize = 10.0
  totalK = 1000.0 * (gridSize^2)
  individuals = 100

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
  #abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 0.0)
  abenv = simplehabitatAE(0.0, grid, totalK, gridSize)
  rel = TraitRelationship(GaussTemp)
  eco = Ecosystem(sppl, abenv, rel)

  burnin = 1000; interval = 1;
  lensim = length(0:interval:times);

  # Run simulations 10 times
  storage = SharedArray(generate_storage(eco, lensim, reps));
  #storage_newsteady = generate_storage(eco, lensim, reps);
  @sync @parallel for j in 1:reps
    if (j != 1) repopulate!(eco) end
    simulate!(eco, burnin, interval, timestep);
    thisstore = view(storage, :, :, :, j);
    simulate_record!(thisstore, eco, times, interval, timestep);
  end
    storage
end
storage = runsim(500, 1, 20.0)
@rput storage

R"for (i in c(1:25)){
if (i==1) plot_fun = plot else plot_fun =lines
plot_fun(1:501, storage[i,1,,1], col=i, type='l',
ylim=c(0,10000))}"

"""
Number of species surviving at equilibrium
"""
reps = 100
times = 500
maxE = linspace(2.0,10.0, 10)
meanSR =
  @sync @parallel (vcat) for k in maxE
    storage = runsim(times, reps, k)
    sum(storage[:, :, (times+1), 1:reps] .> 0)/reps
  end
meanSR
meanAB =
  @sync @parallel (vcat) for k in maxE
    storage = runsim(times, reps, k)
    mapslices(mean, storage[:, :, (times+1), 1:reps], [1, 3])
  end
@rput maxE; @rput meanAB
R"plot(maxE, meanAB, type = 'l', xlab='Maximum energy requirement',
ylab='Average abundance')"


"""
Species richness vs grid
"""

@everywhere function runsim(times::Int64, reps::Int64, squares::Int64,
   move)
  # Set up initial parameters for ecosystem
  numSpecies = 25

  # Set up how much energy each species consumes
  energy_vec = SimpleRequirement(repmat([2.0], numSpecies))
#collect(linspace(2, maxE, 25))
  # Set probabilities
  birth = 0.6/12
  death = 0.6/12
  l = 1.0
  s = 0.5
  boost = 1000.0
  timestep = 1.0

  # Collect model parameters together (in this order!!)
  param = EqualPop(birth, death, l, s, boost)

  minT = 0.0
  maxT = 0.0

  grid = (squares, squares)
  gridSize = 10.0/squares
  totalK = 1000.0 * (gridSize^2)
  individuals = 100

  # Create ecosystem
  kernel = GaussianKernel(0.5, numSpecies, 10e-4)
  movement = move(kernel)

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
  #abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 0.0)
  abenv = simplehabitatAE(0.0, grid, totalK, gridSize)
  rel = TraitRelationship(GaussTemp)
  eco = Ecosystem(sppl, abenv, rel)

  burnin = 1000; interval = 1;
  lensim = length(0:interval:times);

  # Run simulations 10 times
  storage = SharedArray(generate_storage(eco, lensim, reps));
  #storage_newsteady = generate_storage(eco, lensim, reps);
  @sync @parallel for j in 1:reps
    if (j != 1) repopulate!(eco) end
    simulate!(eco, burnin, interval, timestep);
    thisstore = view(storage, :, :, :, j);
    simulate_record!(thisstore, eco, times, interval, timestep);
  end
    storage
end
storage = runsim(500, 1, 2)
@rput storage

R"for (i in c(1:25)){
if (i==1) plot_fun = plot else plot_fun =lines
plot_fun(1:501, storage[i,1,,1], col=i, type='l',
ylim=c(0,10000))}"

"""
Number of species surviving at equilibrium
"""
reps = 10
times = 100
squares = 2:10
meanAB =
  @sync @parallel (vcat) for k in squares
    storage1 = runsim(times, reps, k, NoMovement)
    storage2 = runsim(times, reps, k, BirthOnlyMovement)
    m1 = sum(storage1[:, :, (times+1), 1:reps] .> 0)/reps
    m2 = sum(storage2[:, :, (times+1), 1:reps] .> 0)/reps
    hcat(m1, m2)
  end
@rput squares; @rput meanAB
R"plot(squares, meanAB[,1]/10, type = 'l', xlab='Number of grid squares',
ylab='Average species richness', xaxs='i', yaxs='i', ylim=c(0, 25));
lines(squares, meanAB[,2]/10, col=2)"

m=copy(meanAB)


@everywhere function runsim(times::Int64, reps::Int64, move::Float64)
  # Set up initial parameters for ecosystem
  numSpecies = 25

  # Set up how much energy each species consumes
  energy_vec = SimpleRequirement(repmat([2.0], numSpecies))
#collect(linspace(2, maxE, 25))
  # Set probabilities
  birth = 0.6/12
  death = 0.6/12
  l = 1.0
  s = 0.0
  boost = 1000.0
  timestep = 1.0

  # Collect model parameters together (in this order!!)
  param = EqualPop(birth, death, l, s, boost)

  minT = 0.0
  maxT = 0.0

  grid = (10, 10)
  gridSize = 1.0
  totalK = 100000.0 * (gridSize^2)
  individuals = 1000

  # Create ecosystem
  kernel = GaussianKernel(move, numSpecies, 10e-4)
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
  #abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 0.0)
  abenv = simplehabitatAE(0.0, grid, totalK, gridSize)
  rel = TraitRelationship(GaussTemp)
  eco = Ecosystem(sppl, abenv, rel)

  burnin = 500; interval = 1;
  lensim = length(0:interval:times);

  # Run simulations 10 times
  storage = SharedArray(generate_storage(eco, lensim, reps));
  #storage_newsteady = generate_storage(eco, lensim, reps);
  @sync @parallel for j in 1:reps
    if (j != 1) repopulate!(eco) end
    simulate!(eco, burnin, interval, timestep);
    thisstore = view(storage, :, :, :, j);
    simulate_record!(thisstore, eco, times, interval, timestep);
  end
    storage
end
storage = runsim(100, 1, 0.2)
@rput storage

R"for (i in c(1:25)){
if (i==1) plot_fun = plot else plot_fun =lines
plot_fun(1:501, storage[i,1,,1], col=i, type='l',
ylim=c(0,10000))}"

"""
Number of species surviving at equilibrium
"""
reps = 1
times = 100
moves = linspace(0.001, 0.2, 10)
meanAB =
  @sync @parallel (vcat) for k in moves
    storage = runsim(times, reps, k)
    mapslices(mean, storage[:, :, (times+1), 1:reps], [1, 2, 3])
  end
@rput moves; @rput meanAB
R"plot(moves, meanAB, type = 'l', xlab='Maximum energy requirement',
ylab='Average abundance')"

@everywhere function runsim(times::Int64, reps::Int64, var::Float64)
  # Set up initial parameters for ecosystem
  numSpecies = 25

  # Set up how much energy each species consumes
  energy_vec = SimpleRequirement(repmat([2.0], numSpecies))
#collect(linspace(2, maxE, 25))
  # Set probabilities
  birth = 0.6/12
  death = 0.6/12
  l = 1.0
  s = 0.1
  boost = 1.0
  timestep = 1.0

  # Collect model parameters together (in this order!!)
  param = EqualPop(birth, death, l, s, boost)

  minT = -10.0
  maxT = 10.0

  grid = (10, 10)
  gridSize = 0.1
  totalK = 100000.0
  individuals = 1000

  # Create ecosystem
  kernel = GaussianKernel(0.2, numSpecies, 10e-4)
  movement = BirthOnlyMovement(kernel)

  opts = repmat([5.0], numSpecies) #collect(linspace(minT, maxT, 8))
  vars = repmat([var], numSpecies)
  traits = TempTrait(opts, vars)
  abun = Multinomial(individuals, numSpecies)
  names = map(x -> "$x", 1:numSpecies)
  sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
   movement, param)
  abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 0.0)
  rel = TraitRelationship(GaussTemp)
  eco = Ecosystem(sppl, abenv, rel)

  burnin = 200; interval = 1;
  lensim = length(0:interval:times);

  # Run simulations 10 times
  storage = generate_storage(eco, lensim, reps);
  #storage_newsteady = generate_storage(eco, lensim, reps);
  for j in 1:reps
    if (j != 1) repopulate!(eco) end
    simulate!(eco, burnin, interval, timestep);
    thisstore = view(storage, :, :, :, 1);
    simulate_record!(thisstore, eco, times, interval, timestep);
  end
  return storage
end
runsim(1, 1, 0.01)
"""
SPATIAL ALPHA PLOT
"""
times = 200
vars = [0.1, 0.5, 1.0, 5.0]
for k in vars
  storage = runsim(times, 1, k)
  met = Metacommunity(storage[:, :, 1, 1], UniqueTypes(25))
  alphas = norm_sub_alpha(met, 1)[:diversity]
  im = reshape(alphas, 10, 10)
  @rput k; @rput im
R"library(fields);library(viridis);par(mfrow=c(1,1))
jpeg(paste('plots/progress/var', k, '.jpg'), quality=100)
im[is.na(im)] = 0
image.plot(1:10, 1:10, t(im),col=magma(30), breaks = c(0:30),
xlab='', ylab='', axes=F)
dev.off()
"
end
times = 200
vars = [0.1, 0.5, 1.0, 5.0]
im = zeros(4, 10)
for k in vars
  storage = runsim(times, 1, k)
  store = reshape(storage, 25, 10, 10, (times+1), 1)
  im[findin(vars, k), :] = mapslices(sum, store[:,:,:, (times+1), 1], [1, 3])
end
@rput im; @rput vars
R"library(fields);library(viridis);par(mfrow=c(1,1))
for (i in c(1:4)){
  if (i ==1 ) plot_fun = plot else plot_fun=lines
 plot_fun(1:10, im[i, ], col=i, xlab='Northness', ylab = 'Total abundance',
 type='l', ylim=c(0, 3000), xaxs='i', yaxs='i')
  legend('topright', legend=vars, col= c(1:4), pch='-', cex=0.9)}
"
