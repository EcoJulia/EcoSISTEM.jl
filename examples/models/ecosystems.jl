# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The landscapes and communities the experiments run on, built with
# `GridHabitat`/`build_species`/`build_ecosystem`.
#
# **The published configuration, unchanged.** Supply is per unit area and demand per unit of
# plant, so a species' actual demand is `DEMAND .* size` — which is what lets the large-pool
# experiment give every species its own body size. These are the paper's numbers: 4.5e11 kJ per km²
# per day against 450,000 kJ per m² per day, so a 1 m² plant needs 450,000 kJ/day and a 1 km² cell
# supports a million of them.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

const SUPPLY = (sunlight = 4.5e11kJ / (km^2 * day), water = 192.0mm / day)
const DEMAND = (sunlight = 450000.0kJ / (m^2 * day),
                water = 192.0Unitful.L / (m^2 * day))

# A square study area of total `area`, divided into `cells × cells`. Cell size follows from the
# two rather than being fixed, which is what lets the grid-resolution experiment carve one landscape
# four ways while the area experiment grows the landscape at a fixed grid.
function _studyarea(cells, area)
    side = sqrt(area)
    return StudyArea(extent = (side, side), cellsize = side / cells,
                     verbosity = :silent)
end

# `topology` is an **environment** keyword since step 19 — whether the grid wraps is a property
# of the map, not of the species dispersing over it. `Torus()` here because these experiments
# were written against a wrapping grid; the ones that want an edge ask for `Island()`.
"""
    uniform_environment(cells, area; temperature = 298.0K)

Build a `cells × cells` environment over `area`, uniform in every respect: one temperature
everywhere, with uniform solar and water supply.

 The regime is a **single, unnamed** temperature layer, as the original was — a one-element named
regime is refused, deliberately, since there would be nothing to check it against. The supplies are
named, so a gradient can say which resource it grades.
"""
function uniform_environment(cells, area; temperature = 298.0K,
                             topology = Torus())
    return GridHabitat(regime = UniformSpec(temperature,
                                            axis = Temperature),
                       supply = (sunlight = UniformSpec(SUPPLY.sunlight,
                                                        axis = SolarRadiation),
                                 water = UniformSpec(SUPPLY.water,
                                                     axis = Precipitation)),
                       area = _studyarea(cells, area),
                       topology = topology)
end

"""
    community(env, opts, widths; kw...)

Put `length(opts)` species on `env`, with Gaussian temperature tolerances (`opts` the optima,
`widths` the niche widths) and Gaussian dispersal of mean distance `dispersal`.

Every per-species keyword takes a scalar or a length-`n` vector, so an experiment can vary body
size, dispersal or mortality across the pool without needing a different builder — which is what
replaced the hand-built `PopGrowth`/`SpeciesRequirementCollection` of the original. `individuals` is a
**total**, split across species from the seeded stream, so the starting state is identical on every
run.
"""
function community(env, opts, widths; dispersal = 2.4km,
                   individuals = 100_000_000,
                   birth = 0.15 / year, death = 0.15 / year, sizes = 1.0m^2,
                   seed = 1)
    species = build_species(length(opts), tolerance = (opts, widths),
                            toleranceaxis = Temperature,
                            demand = map(d -> d .* sizes, DEMAND),
                            demandaxis = (sunlight = SolarRadiation,
                                          water = Precipitation),
                            dispersal = dispersal,
                            birth = birth, death = death, survival = 0.1,
                            abundance = individuals, seed = seed)
    return build_ecosystem(species, env, seed = seed)
end

# Total abundance per cell, laid back out on the grid — `(y, x)`, as everywhere in the package.
percell(eco) = dropdims(sum(eco.abundances.grid, dims = 1), dims = 1)

# Total abundance per species, summed over the whole landscape.
perspecies(eco) = vec(sum(eco.abundances.matrix, dims = 2))
