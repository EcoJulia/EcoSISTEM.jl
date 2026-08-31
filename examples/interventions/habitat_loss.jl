# SPDX-License-Identifier: LGPL-3.0-or-later
#
# `RandHabitatLoss!` and `ClustHabitatLoss!` from `examples/models/scenarios.jl`, recreated.
#
# **The reproducibility difference is the headline.** Both originals drew from the **global** RNG
# (`sample(pos, howmany)`), so no run using one could be reproduced, and under MPI every rank chose
# *different* cells and silently diverged. `RandomCells`/`SpreadingCells` draw from a counter-based
# stream keyed on `(seed, :intervention, k, step)`, so a run replays exactly and every rank agrees
# without communicating. **New runs therefore will not reproduce old numbers** - the old ones were
# not reproducible at all, which is the point.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

function loss_ecosystem(; numspecies = configuration().numspecies,
                        individuals = configuration().individuals, seed = 1)
    env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                      supply = UniformSpec(supply_density(individuals),
                                           axis = SolarRadiation),
                      area = study_area())
    species = build_species(numspecies, tolerance = (285.0K, 5.0K),
                            toleranceaxis = Temperature,
                            demand = intervention_demand(),
                            demandaxis = SolarRadiation,
                            abundance = individuals,
                            seed = seed)
    return build_ecosystem(species, env, seed = seed)
end

# ---------------------------------------------------------------------------
# `RandHabitatLoss!` / `ClustHabitatLoss!` - cells lost at random, or spreading
# ---------------------------------------------------------------------------
# **One line each, and faithful.** Both originals zeroed the cells' supplies and abundances and
# then deactivated; `Deactivate` now *is* that - destroying a cell kills what lives there, because a
# deactivated cell is skipped by the hot loop entirely and anything left in it would be frozen rather
# than dead. (Reactivating does not bring them back: it makes the cell habitable again so dispersal
# can recolonise it, as vegetation returns to a slag heap.)
#
# **And the rate is the rate**, not a fixed count: a `Quantity` region count draws binomially over
# the step, which is the process the originals had (`jbinom(1, npos, rate)`) - but from a stream
# keyed on the step rather than the global RNG, so it replays exactly.
function random_loss(rate = 0.05 / year)
    return Intervention(EveryStep(), RandomCells(rate),
                        Deactivate())
end
function clustered_loss(rate = 0.05 / year)
    return Intervention(EveryStep(),
                        SpreadingCells(rate),
                        Deactivate())
end

# ---------------------------------------------------------------------------
# Land conversion - clear an area, then plant a crop on it
# ---------------------------------------------------------------------------
# Several operations over **one** resolved region, which is what the varargs are for: two separate
# interventions would each draw their own random cells, so the crop would be planted somewhere other
# than the cleared ground.
function convert_to_crop(crop, seeded = 500; rate = 0.02 / year)
    return Intervention(EveryStep(), RandomCells(rate), Deactivate(),
                        AddAbundance(crop, seeded))
end
