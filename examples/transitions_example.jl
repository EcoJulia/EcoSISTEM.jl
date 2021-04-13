using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using Diversity

# Set up initial parameters for ecosystem
numSpecies = 10; grid = (10, 10); req= 10.0kJ; individuals=10_000; area = 1000.0*km^2; totalK = 1000.0kJ/km^2

# Set up how much energy each species consumes
energy_vec = SolarRequirement(fill(req, numSpecies))


# Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.0
boost = 1000.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
movement = AlwaysMovement(kernel, Torus())


# Create species list, including their temperature preferences, seed abundance and native status
opts = fill(274.0K, numSpecies)
vars = fill(0.5K, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
sppl.params.birth

# Create abiotic environment - even grid of one temperature
abenv = simplehabitatAE(274.0K, grid, totalK, area)


# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

#Create ecosystem
eco = Ecosystem(sppl, abenv, rel)

# Simulation Parameters
burnin = 5years; times = 50years; timestep = 1month; record_interval = 3months; repeats = 1
lensim = length(0years:record_interval:times)
abuns = zeros(Int64, numSpecies, prod(grid), lensim)
# Burnin
@time simulate!(eco, burnin, timestep);
@time simulate_record!(abuns, eco, times, record_interval, timestep);

eco = Ecosystem(sppl, abenv, rel)
@time biodiversity_simulate!(eco, burnin, timestep);
@time biodiversity_simulate_record!(abuns, eco, times, record_interval, timestep);

# Benchmark
using BenchmarkTools
eco = Ecosystem(sppl, abenv, rel);
@benchmark simulate!(eco, burnin, timestep)
eco = Ecosystem(sppl, abenv, rel);
@benchmark biodiversity_simulate!(eco, burnin, timestep)

using Plots
@gif for i in 1:lensim
    heatmap(abuns[:, :, i], clim = (50, 150))
end

using ProfileView
eco = Ecosystem(sppl, abenv, rel)
@profview simulate!(eco, burnin, timestep)
