# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Scale and the shared landscape for the intervention scenarios.
#
# **The suite runs these small; a direct run gets the published size.** `test/runtests.jl` sets
# `ECOSISTEM_SCALE=small`, so anything mediated through it is a test - nothing here has to detect
# whether it is under test. Unset means full scale, which is what
# `julia --project=examples examples/interventions.jl` gets. `ECOSISTEM_SCALE=large` overrides the
# suite's default in the other direction.
#
# **Included by the runners, not by the scenario files.** `examples/interventions.jl` and
# `examples/models.jl` each include this once, before the three scenario files that use it - the same
# arrangement as `examples/models/`. Including it from each scenario file instead would define these
# constants three times over.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

runscale() = Symbol(get(ENV, "ECOSISTEM_SCALE", "large"))

"""
    CONFIGURATIONS

Species count and total individuals at each scale. The `large` row is the published configuration
these scenarios were originally run at - `test/examples/Island.jl` and its siblings put **100
species and 100,000,000 individuals** on this same 100 km^2 grid.

 Only the population scales; the **grid does not**. The originals were 10 × 10 at 1 km too, so
there is nothing to shrink there, and holding it fixed keeps the cell-selecting scenarios
(`RandomCells`, `SpreadingCells`) drawing from the same number of cells at both scales.
"""
const CONFIGURATIONS = (small = (numspecies = 10, individuals = 10_000),
                        large = (numspecies = 100,
                                 individuals = 100_000_000))

# The configuration at the scale currently in force.
configuration() = CONFIGURATIONS[runscale()]

# **Demand is the published per-individual figure, and supply is derived from it** so that every
# configuration starts at its carrying capacity. This is not tidiness: scaling the population
# without scaling the supply leaves capacity unchanged, and the run collapses back to it. Measured
# while doing exactly that to `examples/landscapes/` - an island seeded with 100,000,000 individuals
# against the small configuration's supply crashed to **124**.
# Underscore-prefixed, and that is not decoration. `examples/models.jl` includes both
# `models/ecosystems.jl` and these files into **one** module, and that file has its own `DEMAND`
# (a named tuple over two resources) and its own `CELLS`/`AREA`. Julia 1.12 lets the second
# definition silently overwrite the first - no warning at all - so an unprefixed `DEMAND` here
# replaced the models' one and broke its `map(d -> d .* sizes, DEMAND)`. Measured, not guessed.
const _AREA = 100.0km^2
const _DEMAND = 4.5e5kJ / day
intervention_demand() = _DEMAND
supply_density(individuals) = individuals * _DEMAND / _AREA

# The shared 10 × 10 km grid at 1 km cells, which all three scenarios sit on.
function study_area()
    return StudyArea(extent = (10.0km, 10.0km), cellsize = 1.0km,
                     verbosity = :silent)
end
