using Simulation
using Base.Test
using Distributions
using Unitful.DefaultSymbols
using RCall

numSpecies=4
numTraits = 1
numNiches = 2
# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2.0], numSpecies))

# Set probabilities
birth = 6.0/year
death = 6.0/year
l = 1.0
s = 0.5
boost = 1.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, l, s, boost)

grid = (10,10)
area = 100.0km^2
totalK = 100000.0
individuals=10000

# Create ecosystem
kernel = GaussianKernel(0.1km, numSpecies, 10e-4)
movement = AlwaysMovement(kernel)
abun = Multinomial(individuals, numSpecies)

# Test out discrete trait update
sppl = SpeciesList(numSpecies, numTraits, abun,
                   energy_vec, movement, param)
abenv = simplenicheAE(numNiches, grid, totalK, area)
rel = TraitRelationship{eltype(abenv.habitat)}(SimpleNiche)
eco = Ecosystem(sppl, abenv, rel)

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 100month; interval = 1month; reps = 1;
    lensim = length(0month:interval:times);

    # Run simulations 10 times
    storage = generate_storage(eco, lensim, reps);
    for j in 1:reps
      if (j != 1) repopulate!(eco) end
      simulate!(eco, burnin, interval, timestep);
      thisstore = view(storage, :, :, :, j);
      simulate_record!(thisstore, eco, times, interval, timestep);
    end
    return storage
end

abun = runsim(eco, 10year)
plot_abun(abun, 4, grid)

meanabun = reshape(mapslices(mean, abun, [3,4])[:,:, 1,1], 4, 10, 10)
map(1:4) do spp
    [mean(meanabun[spp,:, :][hab.==getpref(DiscreteHab{Simulation.Niches}, eco, spp)]),
    mean(meanabun[spp,:, :][hab.!=getpref(DiscreteHab{Simulation.Niches}, eco, spp)])]
end
im = reshape(abun[:,:, end, 1], 4, 10, 10)
hab = eco.abenv.habitat.matrix
map(1:4) do spp
    [mean(im[spp,:,:][hab.==getpref(eco.spplist.traits, spp)]),
    mean(im[spp,:,:][hab.!=getpref(eco.spplist.traits, spp)])]
end
@rput meanabun
@rput hab
R"par(mfrow=c(1,2))
library(viridis); library(fields)
image.plot(meanabun[2,,], col=magma(50))
image.plot(hab, col=viridis(2))"
