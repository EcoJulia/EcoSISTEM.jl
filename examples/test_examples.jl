# SPDX-License-Identifier: LGPL-3.0-or-later

## Test an example run using the examples project folder

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Diversity
using OnlineStats
using Plots
using Distributions

## DIFFERENT OPTS ##

numSpecies = 100;
grd = (10, 10);
demand = (450000.0kJ / m^2, 192.0nm / m^2);
individuals = 100_000_000;
area = 100.0 * km^2;
totalK = (4.5e11kJ / km^2, 192.0mm / km^2)

habitat1 = simplehabitat(298.0K, grd, totalK[1], area)
habitat2 = simplehabitat(298.0K, grd, totalK[2], area)
supply = SupplyCollection2(habitat1.supply, habitat2.supply)
habitat = GridHabitat{typeof(habitat1.regime),
                      typeof(supply)}(habitat1.regime,
                                      habitat1.active,
                                      supply,
                                      habitat1.names)

vars = fill(2.0K, numSpecies)
opts = 298.0K .+ vars .* range(-3, stop = 3, length = numSpecies)

av_dist = fill(2.4km, numSpecies)
kernel = GaussianKernel.(av_dist, 10e-10)

death = 0.15 / year
birth = death
long = 1.0
surv = 0.1
boost = 1.0

size_mean = 1.0m^2
# Set up how much resource each species consumes
resource_vec1 = SolarDemand(fill(demand[1] * size_mean, numSpecies))
resource_vec2 = WaterDemand(fill(demand[2] * size_mean, numSpecies))

resource_vec = DemandCollection2(resource_vec1, resource_vec2)
param = EqualPop(birth, death, long, surv, boost)

# Create ecosystem

movement = BirthOnlyMovement(kernel, Torus())

tolerance = NicheTolerance(Temperature, Normal, opts, vars)
native = fill(true, numSpecies)
abun = rand(Multinomial(individuals, numSpecies))
sppl = SpeciesList(numSpecies, tolerance, abun, resource_vec,
                   movement, param, native)
nichefit = NicheSuitability{typeof(first(opts))}()
eco = Ecosystem(sppl, habitat, nichefit)

times = 10years;
timestep = 1month;
lensim = length((0month):timestep:times)

simulate!(eco, times, timestep)
