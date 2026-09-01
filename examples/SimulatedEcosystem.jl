# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A fully synthetic ecosystem - no real climate/land-cover data, no downloads. A temperature
# gradient (the simulated "regime") and a flat solar resource (the simulated "resource") over a
# deliberately non-square grid, populated with random species (random niche preferences and
# dispersal distances) at random initial abundances, then simulated for several months.
#
# Everything is seeded, so re-running this script reproduces exactly the same result.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using Random

# --- the study area -----------------------------------------------------------
# The grid is decided first, on its own, and only then is anything built on it. A purely synthetic
# area needs no layers at all - its geometry *is* `extent` + `cellsize`, so there is no data whose
# CRS, resolution or footprint could have a say.
#
# A non-square physical extent (20km north-south × 30km east-west) divided into 1km cells, giving
# a 20 × 30 grid - deliberately non-square, so this example also exercises the y/x dimension-order
# fix directly, not just a square grid where an order mixup wouldn't show up. `within` restricts the
# active cells to a disc of radius 4km about the grid's centre, inactivating the rest.
area = StudyArea(extent = (20.0km, 30.0km), cellsize = 1.0km,
                 within = CircleMaskSpec(radius = 4.0km))

# --- the simulated environment ----------------------------------------------
# The simulated regime: a temperature gradient from 10°C (bottom) to 30°C (top), tagged with the
# `Temperature` niche axis so species niches (below) are matched to it by name.
regime = GradientSpec(10°C, 30°C, axis = Temperature, orientation = 45°)

# The simulated resource: a flat solar supply available in every active cell.
supply = UniformSpec(5000.0kJ / km^2 / day, axis = SolarRadiation)

habitat = GridHabitat(regime = regime, supply = supply, area = area)

# --- random species ----------------------------------------------------------
seed = 1
Random.seed!(seed)

numspecies = 12

# Each species gets its own random temperature niche (mean drawn across the gradient's own
# range, so every niche is genuinely reachable somewhere on the grid) and its own random
# dispersal distance - `build_species` accepts a per-species vector for any of these keywords.
# A tolerance's mean (a position) and width (an interval) must share one unit, so the °C means
# (drawn across the gradient's 10-30 °C range) are converted to the K frame the widths use.
niche_means = uconvert.(K, rand(Uniform(10, 30), numspecies) * °C)
niche_widths = rand(Uniform(1.0, 4.0), numspecies) * K
dispersal = rand(Uniform(1.0, 5.0), numspecies) * km

species = build_species(numspecies, tolerance = (niche_means, niche_widths),
                        toleranceaxis = Temperature, dispersal = dispersal,
                        demand = 5.0kJ / day, demandaxis = SolarRadiation,
                        # A total of 5000 individuals, split at random across species - genuinely
                        # random initial abundances, reproducible via `seed`.
                        abundance = 5000, seed = seed)

# --- assemble and simulate ----------------------------------------------------
eco = build_ecosystem(species, habitat, seed = seed)

times = 10year
timestep = 1month_mean_duration

startabun = sum(eco.abundances.matrix)
println("Per-species initial abundance: ",
        sum(eco.abundances.matrix, dims = 2)[:, 1])

simulate!(eco, times, timestep)
endabun = sum(eco.abundances.matrix)

println("Simulated $numspecies species over a $(size(eco.abundances.grid, 2)) × " *
        "$(size(eco.abundances.grid, 3)) grid ($(count(habitat.active)) active cells) " *
        "for $times.")
println("Total abundance: $startabun -> $endabun.")
println("Per-species final abundance: ",
        sum(eco.abundances.matrix, dims = 2)[:, 1])
