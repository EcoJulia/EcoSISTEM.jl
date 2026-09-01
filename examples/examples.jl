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

# **The gate, added 2026-08-13 - this file had none and was the second-slowest thing in
# `extras_examples` at 36 s, a fifth of the whole run.** It is a *smoke test* of the two-resource
# path, so under the suite it needs enough individuals to exercise the arithmetic and no more:
# 100 million of them on 100 cells was buying nothing a million does not. Run it directly (or set
# `ECOSISTEM_SCALE=large`) for the published figures.
const SMALL = get(ENV, "ECOSISTEM_SCALE", "large") == "small"

numSpecies = 100;
grd = (10, 10);
demand = (450000.0kJ / m^2 / day, 192.0Unitful.L / m^2 / day);
individuals = SMALL ? 1_000_000 : 100_000_000;
area = 100.0 * km^2;
totalK = (4.5e11kJ / km^2 / day, 192.0mm / day)

# A `StudyArea` decides the grid, then `GridHabitat` builds on it. With no data-backed layer
# to shape it, `extent` and `cellsize` give a **synthetic** grid with no CRS - which is exactly what a
# synthetic `UniformSpec` needs, since it is generated at whatever shape it is handed. A square
# `area` split `grd[1]` ways gives the same 10x10 grid of 1 km^2 cells the old `simplehabitat(val, grd,
# maxsupply, area)` derived, so nothing about the simulation below changes.
#
# `verbosity = :silent` because this example runs under `@test_nowarn`, which fails on *any* stderr
# output - and the default `:normal` announces each guessed value. Worth knowing outside a test:
# leave it at the default and it tells you what it decided.
side = sqrt(area)
studyarea = StudyArea(extent = (side, side), cellsize = side / grd[1],
                      verbosity = :silent)
# Both resources in one call, so the two-resource environment is **built** rather than assembled.
# The `simplehabitat` version had to build one single-resource habitat per resource and then reach
# inside them to staple the supplies into a `LayerCollection` by hand.
habitat = GridHabitat(regime = UniformSpec(298.0K, axis = Temperature),
                      supply = (UniformSpec(totalK[1],
                                            axis = SolarRadiation),
                                UniformSpec(totalK[2],
                                            axis = Precipitation)),
                      area = studyarea,
                      topology = Torus())

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
resource_vec1 = Demand{SolarRadiation}(fill(demand[1] * size_mean, numSpecies))
resource_vec2 = Demand{Precipitation}(fill(demand[2] * size_mean, numSpecies))

resource_vec = SpeciesRequirementCollection((resource_vec1, resource_vec2))
param = EqualPop(birth, death, long, surv, boost)

# Create ecosystem

movement = BirthOnlyMovement(kernel)

tolerance = NicheTolerance(Temperature, Normal, opts, vars)
native = fill(true, numSpecies)
abun = rand(Multinomial(individuals, numSpecies))
sppl = SpeciesList(numSpecies, tolerance, abun, resource_vec,
                   movement, param, native)
nichefit = NicheSuitability{Temperature, typeof(first(opts))}()
eco = Ecosystem(sppl, habitat, nichefit)

times = SMALL ? 2year : 10year;
timestep = 1month_mean_duration;
lensim = length((0month_mean_duration):timestep:times)

simulate!(eco, times, timestep)
