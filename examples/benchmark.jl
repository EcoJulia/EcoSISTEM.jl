## Code for running model on spatial grid

# Include simulation functions and other modules
using Diversity
using Simulation
using RCall
using Distributions
using StatsBase
using ProgressMeter
using DataStructures
#using ArrayViews
using DataFrames
using BenchmarkTools
using ProfileView

# Set up initial parameters for ecosystem
numSpecies = 25
numTraits = 2
numNiches = 2

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2.0], numSpecies))

# Set probabilities
birth = 0.6/12
death = 0.6/12
l = 1.0
s = 0.5
boost = 1.05
timestep = 1/12

# Collect model parameters together
param = EqualPop(birth, death, l, s, boost)

minT = 0.0
maxT = 10.0

grid = (10, 10)
gridSize = 0.1
totalK = 100000.0
individuals=10000

# Create ecosystem
kernel = GaussianKernel(0.2, numSpecies, 10e-4)
movement = BirthOnlyMovement(kernel)

#opts = repmat([5.0], numSpecies) #collect(linspace(minT, maxT, 8))
#vars = collect(1.0:1.0:8.0)[randperm( 8)]
temps = readtable(Pkg.dir("Simulation", "examples", "Tempdat.csv"))[1:numSpecies, :]
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
abenv = tempgradAE(minT, maxT, grid, totalK, gridSize, 0.01)
rel = TraitRelationship(GaussTemp)
eco = Ecosystem(sppl, abenv, rel)

function runsim(ecosys::Ecosystem, times::Int64)

  interval = 1; reps = 1;
  lensim = length(0:interval:times);

  # Run simulations 10 times
  storage = generate_storage(ecosys, lensim, reps);
  #storage_newsteady = generate_storage(eco, lensim, reps);
  for j in 1:reps
    if (j != 1) repopulate!(ecosys) end
    thisstore = view(storage, :, :, :, j)
    simulate_record!(thisstore, ecosys, times, interval, timestep)
  end
end
runsim(eco, 1)  # run once to trigger compilation
Profile.clear()  # in case we have any previous profiling data
@profile runsim(eco, 100)
#Profile.print(format := flat)
ProfileView.view()
