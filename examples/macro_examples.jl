## Code for running model on spatial grid
#addprocs(5)
# Include simulation functions
using Diversity
using Simulation
#@everywhere using Simulation
using RCall
using Distributions
using StatsBase
using ProgressMeter
using DataStructures
using ArrayViews
## Run simulation on 2 by 1 grid - investigate spatial distributions

# Set up initial parameters for ecosystem
numSpecies=4
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2.0], numSpecies))

# Set probabilities
birth = 0.6
death = 0.6
l = 1.0
s = 0.9
timestep = 1.0

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, l, s)

minT = -10.0
maxT = 20.0

grid = (10, 10)
gridSize = 1.0
totalK = 100000.0
individuals=1000

# Create ecosystem
kernel = GaussianKernel(0.2, numSpecies, 10e-4)
movement = BirthOnlyMovement(kernel)

opts = repmat([5.0], numSpecies)#collect(linspace(minT, maxT, 8))
#vars = collect(1.0:1.0:8.0)[randperm( 8)]
vars = collect(linspace(1, 4, numSpecies))
traits = TempTrait(opts, vars)
abun = Multinomial(individuals, numSpecies)
names = map(x -> "$x", 1:numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
 movement, param)
abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 1.0)
rel = TraitRelationship(GaussTemp)
eco = Ecosystem(sppl, abenv, rel)

#plotdiv(norm_meta_alpha, eco, collect(0:10))

times = 100; burnin = 0; interval = 10; reps = 1
lensim = length(0:interval:times)

# Run simulations 10 times
storage = generate_storage(eco, lensim, reps)
@time for j in 1:reps
  if (j != 1) repopulate!(eco) end
  thisstore = @view storage[ :, :, :, j]
  simulate!(eco, burnin, interval, timestep)
  simulate_record!(thisstore, eco, times, interval, timestep)
end

temps = unique(eco.abenv.habitat.matrix)
hab =eco.abenv.habitat.matrix
@rput hab
im = eco.abundances.grid
@rput im; @rput temps; @rput vars
R"library(fields);par(mfrow=c(2,2))
for(i in c(1:4)){
  image.plot(temps, 1:10, im[i, , ], main =paste('Niche width =', vars[i]))
  }"

plotdiv(norm_sub_beta, eco, 1)

datf = norm_sub_beta(eco, 1)
size(datf, 1) == length(eco.abenv.habitat.matrix) ||
  error("Metacommunity measures cannot be plotted as grid")
im = reshape(datf[:diversity], size(eco.abenv.habitat.matrix))
hab = eco.abenv.habitat.matrix
@rput im; @rput hab
R"par(mfrow=c(1,2));library(fields);
im[im==0]=NA
image.plot(im)
image.plot(hab)"


abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 0.1)
rel = TraitRelationship(GaussTemp)
eco = Ecosystem(sppl, abenv, rel)

#plotdiv(norm_meta_alpha, eco, collect(0:10))

times = 10; burnin = 0; interval = 10; reps = 1
lensim = length(0:interval:times)

# Run simulations 10 times
storage = generate_storage(eco, lensim, reps)
for j in 1:reps
  if (j != 1) repopulate!(eco) end
  thisstore = @view storage[ :, :, :, j]
  simulate!(eco, burnin, interval, timestep)
  hab = gethabitat(eco).matrix
  heatmap(hab)
end
