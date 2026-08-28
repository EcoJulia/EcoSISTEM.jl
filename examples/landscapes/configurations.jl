# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The five landscape configurations from the former `test/examples/` — `Island`, `Patch`, `Peaked`,
# `RegionEven` and `RegionGrad` — recreated on the modern `StudyArea`/spec interface.
#
# **They are not five landscapes. They are two independent choices**, and reading the originals
# side by side is what showed it:
#
# | | temperature | boundary | grid (full scale) | species |
# |---|---|---|---|---|
# | `Island`     | uniform 298 K   | `Island(false)` | 10 × 10 |    100 |
# | `Patch`      | uniform 298 K   | `Torus`             | 10 × 10 |    100 |
# | `RegionEven` | uniform 298 K   | `Cylinder`          | 40 × 25 |  1,000 |
# | `RegionGrad` | 293 → 303 K     | `Cylinder`          | 40 × 25 |  1,000 |
# | `Peaked`     | peaked 293–303 K| `Torus`             | 50 × 50 | 10,000 |
#
# `Island` and `Patch` are the **same habitat** — `simplehabitat(298.0K, (10,10), …)` in both. The
# only substantive difference between those two files is the boundary. So what the set really varies
# is **what the environment looks like** and **what happens at the edge**, and the names describe the
# second more than the first.
#
# That is worth stating because it is the comparison the set exists to make: the same ecology under
# a boundary that wraps in both directions (`Torus`), one (`Cylinder`), or neither (`Island`).
#
# **Note the `false` on the island**, which is `disperse_safely` and is what makes it an island at
# all. A bare `Island()` blocks moves that would leave the grid and then **renormalises** the
# survivors (`_drawmoves!`), so nothing is *lost* at the edge — an individual in a corner is
# redistributed inland rather than drowning. That is a *reflecting* boundary, and under it `Island`
# and `Patch` come out statistically indistinguishable (measured: ~1% apart, the sign flipping with
# the seed). `Island(false)` loses them instead, which is the wind-dispersed case and the one the
# name implies: measured, the island settles at **0.68×** the patch's abundance, consistently across
# seeds.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

# ---------------------------------------------------------------------------
# Scale
# ---------------------------------------------------------------------------
# **The suite runs these small; a direct run gets the originals' full size.** `test/runtests.jl`
# sets `ECOSISTEM_SCALE=small`, so anything mediated through it is a test — and nothing else has to
# work out whether it is under test. Unset means full scale, which is what
# `julia --project=examples examples/landscapes.jl` gets.
#
# Override in either direction: `ECOSISTEM_SCALE=large` runs the big configurations under the
# suite, `ECOSISTEM_SCALE=small` keeps a direct run quick.
runscale() = Symbol(get(ENV, "ECOSISTEM_SCALE", "large"))

"""
    CONFIGURATIONS

Grid, species count, total individuals and landscape area for each of the five landscapes at each
scale. The `large` row is the original published configuration; `small` keeps the same
shape at a size CI can afford.

 `area` is the **whole landscape**, so the cell size follows from it and the grid — which is how
the originals were written, and why `Peaked` has 9 km cells while `Island` has 1 km ones.
"""
const CONFIGURATIONS = (island = (small = (grid = (10, 10), numspecies = 20,
                                           individuals = 10_000,
                                           area = 100.0km^2),
                                  large = (grid = (10, 10), numspecies = 100,
                                           individuals = 100_000_000,
                                           area = 100.0km^2)),
                        patch = (small = (grid = (10, 10), numspecies = 20,
                                          individuals = 10_000,
                                          area = 100.0km^2),
                                 large = (grid = (10, 10), numspecies = 100,
                                          individuals = 100_000_000,
                                          area = 100.0km^2)),
                        region_even = (small = (grid = (20, 12),
                                                numspecies = 20,
                                                individuals = 10_000,
                                                area = 240.0km^2),
                                       large = (grid = (40, 25),
                                                numspecies = 1_000,
                                                individuals = 10^10,
                                                area = 10_000.0km^2)),
                        region_gradient = (small = (grid = (20, 12),
                                                    numspecies = 20,
                                                    individuals = 10_000,
                                                    area = 240.0km^2),
                                           large = (grid = (40, 25),
                                                    numspecies = 1_000,
                                                    individuals = 10^10,
                                                    area = 10_000.0km^2)),
                        peaked = (small = (grid = (15, 15), numspecies = 20,
                                           individuals = 10_000,
                                           area = 225.0km^2),
                                  large = (grid = (50, 50), numspecies = 10_000,
                                           individuals = 10^10,
                                           area = 200_000.0km^2)))

# The configuration for one landscape at the scale currently in force.
_configuration(name) = CONFIGURATIONS[name][runscale()]

# ---------------------------------------------------------------------------
# The environments
# ---------------------------------------------------------------------------
# **The energetics are the originals'**, and they are what make a scale change meaningful rather
# than merely bigger. The original `Island` supplied `4.5e11 kJ/km²/day` to 100,000,000
# individuals over 100 km², which is 450,000 kJ/day each — so that is the per-individual demand here,
# and the supply is derived from it so that **every configuration starts at its carrying capacity**.
#
# Getting this wrong is not subtle but it *is* invisible at one scale. Scaling `individuals` to
# the originals' 100,000,000 while leaving the supply at the small configuration's value left the
# carrying capacity ten thousand times too low: measured, the island crashed from 100,000,000 to
# **124** individuals and the island/patch ratio moved from 0.68 to 0.175 — enough to fail the
# runner's assertion, and only visible by actually running the large path.
# Named for the resource, not just `DEMAND`, because two other example files define a
# `DEMAND` that is a **two-resource named tuple** (`(sunlight = …, water = …)`) rather than a
# scalar. Nothing breaks — no two of them are ever included into one session — but the shared
# name genuinely misled a sweep that inferred each demand's axis from the constant it was
# written as, and gave this call a two-member `demandaxis` for a one-member demand.
const SOLARDEMAND = 4.5e5kJ / day
_supplydensity(individuals, area) = individuals * SOLARDEMAND / area

# A square study area of the given total area, divided into `grid` cells — so the cell size follows
# from the two rather than being fixed, which is what lets one configuration table drive both scales.
function _area(grid, area)
    cellsize = sqrt(area / prod(grid))
    return StudyArea(extent = (grid[1] * cellsize, grid[2] * cellsize),
                     cellsize = cellsize, verbosity = :silent)
end

"""Uniform temperature — the `Island`, `Patch` and `RegionEven` environment."""
function even_environment(; grid = (10, 10), area = 100.0km^2,
                          supply = _supplydensity(10_000, area),
                          temperature = 298.0K,
                          topology = Island())
    return GridHabitat(regime = UniformSpec(temperature,
                                            axis = Temperature),
                       supply = UniformSpec(supply,
                                            axis = SolarRadiation),
                       area = _area(grid, area),
                       topology = topology)
end

"""A monotone temperature gradient across the grid — the `RegionGrad` environment."""
function gradient_environment(; grid = (40, 25), area = 10_000.0km^2,
                              supply = _supplydensity(10_000, area),
                              from = 293.0K, to = 303.0K,
                              topology = Island())
    return GridHabitat(regime = GradientSpec(from, to,
                                             axis = Temperature),
                       supply = UniformSpec(supply,
                                            axis = SolarRadiation),
                       area = _area(grid, area),
                       topology = topology)
end

"""
A temperature *peak* in the middle of the grid, falling away to either side — the `Peaked`
environment.  Ecologically the interesting one: a gradient has one warm end, whereas a peak has a
warm centre with cool edges on **both** sides, so a warm-adapted species is squeezed rather than
pushed.
"""
function peaked_environment(; grid = (25, 25), area = 200_000.0km^2,
                            supply = _supplydensity(10_000, area),
                            from = 293.0K, to = 303.0K,
                            topology = Island())
    return GridHabitat(regime = PeakedSpec(from, to,
                                           axis = Temperature),
                       supply = UniformSpec(supply,
                                            axis = SolarRadiation),
                       area = _area(grid, area),
                       topology = topology)
end

# ---------------------------------------------------------------------------
# The boundaries — the axis the original file names were really about
# ---------------------------------------------------------------------------
# **The boundary is a property of both the environment and the disperser, and they are different
# things** — which is why they are separate. The **topology**
# (`Island`/`Torus`/`Cylinder`) says how the grid's edges join, so it belongs to the
# environment and is chosen there; **`disperse_safely`** says whether an individual aimed at a dead
# cell is reflected inland or lost, which is a fact about the *disperser* and so belongs here.
# This file is where the two were most visibly conflated: `Island(false)` meant "hard edges
# **and** dispersers are lost", one object carrying a map fact and a species fact at once.
function landscape_species(disperse_safely; numspecies = 20,
                           individuals = 10_000,
                           temperature = 298.0K, seed = 1)
    return build_species(numspecies,
                         tolerance = (temperature, 3.0K),
                         toleranceaxis = Temperature,
                         demand = SOLARDEMAND,
                         demandaxis = SolarRadiation,
                         movement = BirthOnlyMovement,
                         disperse_safely = disperse_safely,
                         abundance = individuals, seed = seed)
end

"""
    landscape(name::Symbol; kw...)

Build one of the five original configurations by name — `:island`, `:patch`, `:region_even`,
`:region_gradient` or `:peaked` — returning a ready-to-run `Ecosystem`.

Size comes from [`CONFIGURATIONS`](@ref) at the scale currently in force, so this is small under the
test suite and the original full size outside it; `numspecies` and `individuals` override it.
"""
function landscape(name::Symbol; scaling = _configuration(name),
                   numspecies = scaling.numspecies,
                   individuals = scaling.individuals, seed = 1)
    grid, area = scaling.grid, scaling.area
    # Supply follows the population, so every scale starts at its carrying capacity.
    supply = _supplydensity(individuals, area)
    # Each branch now names the two facts separately: the grid's **topology**, which goes to the
    # environment, and the species' **`disperse_safely`**, which goes to the species.
    env, safely = if name === :island
        # `false` is `disperse_safely` — dispersers aimed off the edge are lost, not reflected
        # back inland. Without it an "island" behaves exactly like a patch (see the header).
        (even_environment(grid = grid, area = area, supply = supply,
                          topology = Island()), false)
    elseif name === :patch
        (even_environment(grid = grid, area = area, supply = supply,
                          topology = Torus()), true)
    elseif name === :region_even
        (even_environment(grid = grid, area = area, supply = supply,
                          topology = Cylinder()), true)
    elseif name === :region_gradient
        (gradient_environment(grid = grid, area = area, supply = supply,
                              topology = Cylinder()), true)
    elseif name === :peaked
        (peaked_environment(grid = grid, area = area, supply = supply,
                            topology = Torus()), true)
    else
        error("unknown landscape `$name`; the five are :island, :patch, :region_even, " *
              ":region_gradient and :peaked.")
    end
    species = landscape_species(safely, numspecies = numspecies,
                                individuals = individuals, seed = seed)
    return build_ecosystem(species, env, seed = seed)
end

const LANDSCAPES = (:island, :patch, :region_even, :region_gradient, :peaked)
