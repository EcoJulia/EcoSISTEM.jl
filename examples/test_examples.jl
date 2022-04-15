## Test an example runs using the examples project folder

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Diversity
using OnlineStats
using Plots
using Distributions
using Diversity
using Test

## DIFFERENT OPTS ##

numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2);
individuals = 100_000_000; area = 100.0*km^2;
totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
bud = BudgetCollection2(abenv1.budget, abenv2.budget)
abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

vars = fill(2.0K, numSpecies)
opts = 298.0K .+ vars .* range(-3, stop = 3, length = numSpecies)

av_dist = fill(2.4, numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death = 0.15/ year
birth = death
long = 1.0
surv = 0.1
boost = 1.0

size_mean = 1.0m^2
# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

energy_vec = ReqCollection2(energy_vec1, energy_vec2)
param = EqualPop(birth, death, long, surv, boost)

# Create ecosystem

movement = BirthOnlyMovement(kernel, Torus())

traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
abun = rand(Multinomial(individuals, numSpecies))
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
rel = Gauss{typeof(first(opts))}()
eco = Ecosystem(sppl, abenv, rel)

times = 10years; timestep = 1month
lensim = length(0month:timestep:times)

@test_nowarn simulate!(eco, times, timestep)
