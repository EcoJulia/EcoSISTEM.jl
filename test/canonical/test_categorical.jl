# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Canonical results for the **categorical** branch: a `SimpleCategoricalTolerance` matched against a
# `CategoricalRegime`. Synthetic throughout, so - like `test_simulated.jl` - nothing here can move
# because a download did.
#
# **Why this file exists at all.** Until the two categorical tolerances merged, the difference
# between *soft* exclusion (a species does worse outside its classes) and *hard* exclusion (it cannot
# live there) was carried by the **type**: `CategoricalTolerance`/`Match` scored 0.5 outside, and
# `LandCoverTolerance`/`LCmatch` scored 0. It is now a `penalty` **number** on the tolerance, which is
# far more expressive and correspondingly easier to change by accident - a default, a comparison, or
# the exponent it enters the demographics through. Nothing in `canonical/` exercised the categorical
# branch at all, so every one of those changes would have left `reference.toml` untouched.
#
# Two runs, identical but for the penalty, are what makes that detectable: a change that moved both
# by the same amount would still be caught, and one that collapsed the distinction - the exact failure
# a wrong default produces - moves them *together*, which the `soft ≠ hard` assertions below refuse.

module CanonicalCategorical

using Test
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Distributions
using Random

include("canonical.jl")
using .Canonical

# One fixed seed for the whole file, as in `test_simulated.jl`.
const SEED = 20260817

# **Non-square, and coprime with every plausible split** - 7 species on a 7 × 11 grid of 77 cells.
# A square grid has hidden a `(y, x)` transposition in this repository three times, and an
# evenly-dividing one hid a real `gatherabundance` row/column mismatch.
const NY, NX = 7, 11
const NSPECIES = 7
const NNICHES = 3

const AREA = StudyArea(extent = (NY * 1.0km, NX * 1.0km), cellsize = 1.0km,
                       verbosity = :silent)

# **The niche map is drawn ONCE and then re-used by VALUE, and both halves of that matter.**
#
# `NicheSpec` lays its niches out by percolation and weighted sampling from the **global** stream,
# which `Ecosystem`'s own per-species seeding does not cover. Drawing it once, immediately after a
# `Random.seed!`, is what makes it a fixed input rather than a function of whatever else has run in
# the process first - and a canonical file is `include`d into a shared process when re-blessing but
# runs in its own when checking, so "whatever else has run first" genuinely differs between the two.
#
# **But each run still gets its OWN habitat, built from this same map.** An `Ecosystem` *shares* its
# habitat with the caller rather than copying it, so two ecosystems built on one habitat object are
# not independent: simulating the first leaves state the second then inherits. Measured - sharing
# one habitat made the soft and hard totals depend on the order the runs were built and simulated in,
# and on the thread count. Passing the built layer instead is safe, because `GridHabitat` **copies**
# a supplied layer in (`_ownlayer`) rather than aliasing it.
#
# With both in place the fixture is reproducible at 1, 2, 4 and 8 threads and across processes,
# which is the package's own stated requirement - a result must not depend on how work was divided.
const REGIME = begin
    Random.seed!(SEED)
    GridHabitat(regime = NicheSpec(NNICHES, axis = LandCoverTypology),
                supply = UniformSpec(4.5e11kJ / km^2 / day,
                                     axis = SolarRadiation),
                area = AREA,
                topology = Torus()).regime
end

# A categorical ecosystem on a fresh habitat carrying the fixed niche map, with each species
# tolerating one niche. `penalty` is the whole point of the fixture - the weight a species gets in a
# cell whose class it does not tolerate.
function categorical_ecosystem(penalty)
    habitat = GridHabitat(regime = REGIME,
                          supply = UniformSpec(4.5e11kJ / km^2 / day,
                                               axis = SolarRadiation),
                          area = AREA,
                          topology = Torus())
    # Species 1 tolerates niche 1, species 2 niche 2, ... cycling round, so every niche has occupants
    # and no species tolerates them all. Written as bare classes rather than sets, which is the
    # spelling the released `DiscreteTrait` took - the constructor wraps each in a one-element set,
    # and that equivalence is itself part of what this fixture pins.
    classes = [1 + (sp - 1) % NNICHES for sp in 1:NSPECIES]
    tolerance = SimpleCategoricalTolerance(classes, axis = LandCoverTypology,
                                           penalty = penalty)

    demand = Demand{SolarRadiation}(fill(450000.0kJ / m^2 / day * 1.0m^2,
                                         NSPECIES))
    movement = BirthOnlyMovement(GaussianKernel.(fill(2.0km, NSPECIES), 10e-10))
    # `survival` must be non-zero or this fixture is VACUOUS. Suitability enters the demographics
    # as `suitability^±survival`, and `x^0.0 == 1.0` for every `x` including `Inf` - so at
    # `survival = 0` the penalty is ignored entirely and the soft and hard runs would be identical.
    # Asserted below rather than merely commented.
    param = EqualPop(0.15 / year, 0.15 / year, 1.0, 0.1, 1.0)

    rng = Random.Xoshiro(SEED)
    abun = rand(rng, Multinomial(70_000, NSPECIES))
    sppl = SpeciesList(NSPECIES, tolerance, abun, demand, movement, param,
                       fill(true, NSPECIES))
    # Through `build_ecosystem` rather than `Ecosystem(..., nichefit)`, so the run also exercises
    # `_defaultsuitability` **inferring** `CategoricalSuitability` and checking the tolerance's axis
    # against the regime's - a check the categorical branch could not have at all until the merge
    # gave these tolerances an axis.
    return build_ecosystem(sppl, habitat, seed = SEED)
end

# The share of each species' abundance sitting in cells whose class it actually tolerates. This is
# what "hard exclusion" means operationally, and it is the property no blessed total can express.
function own_niche_share(eco)
    grid = eco.abundances.grid
    regime = eco.habitat.regime.matrix
    tolerance = eco.spplist.tolerance
    return map(1:NSPECIES) do sp
        total = sum(@view grid[sp, :, :])
        total == 0 && return 0.0
        inown = sum(grid[sp, y, x]
                    for y in axes(regime, 1), x in axes(regime, 2)
                    if regime[y, x] in EcoSISTEM.getpref(tolerance, sp))
        return inown / total
    end
end

# Build one run and simulate it, returning the finished ecosystem and the cell-1 suitability read
# *before* any demography.
#
# **One ecosystem is built, simulated and finished before the next is built, and that is not a
# style choice.** Measured: holding two live ecosystems and simulating them afterwards made the
# blessed totals depend on the thread count (they agreed at 2, 4 and 8 threads and differed at 1),
# while building and simulating one at a time is identical at all four. It is also exactly what
# `test_simulated.jl` does for its reproducibility check, so it is the supported shape rather than a
# workaround. Whether two concurrent ecosystems *should* be independent is a real question about
# the package, recorded separately - this file must not be the thing that answers it.
function categorical_run(penalty)
    eco = categorical_ecosystem(penalty)
    cell1 = suitability(eco, 1, 1)
    simulate!(eco, 2year, 1month_mean_duration)
    return (eco = eco, cell1 = cell1)
end

@testset "canonical: categorical tolerance" begin
    regimemap = copy(REGIME.matrix)
    # **The fixture's INPUT, blessed.** Every number below is conditional on this map, so pinning
    # the outputs without it would leave the experiment able to change underneath its own results -
    # and the self-consistency checks further down (each run got *this* map) cannot see that, because
    # they compare against the map rather than against a record of it.
    canonical("categorical/regime_map", vec(regimemap))
    soft = categorical_run(0.5)
    hard = categorical_run(0.0)

    # **The penalty arithmetic, pinned without any demography at all.** These two numbers are what
    # a change to how `penalty` reaches `suitability` moves, and they move *before* the simulation
    # does - so if the totals below shift and these do not, the cause is downstream of the tolerance.
    # Cell 1 and species 1: whether species 1 tolerates cell 1's class depends on the seeded niche
    # map, so which of the two values each is cannot be asserted here - but it is blessed, and the
    # `in {1.0, penalty}` property below says it can only ever be one of them.
    canonical("categorical/suitability_soft_cell1", soft.cell1)
    canonical("categorical/suitability_hard_cell1", hard.cell1)
    @test soft.cell1 in (1.0, 0.5)
    @test hard.cell1 in (1.0, 0.0)

    softabun, hardabun = soft.eco.abundances.matrix, hard.eco.abundances.matrix

    # Blessed totals for both runs. Per-species as well as grand, so a change that redistributes
    # abundance between species while conserving the total is still caught.
    canonical("categorical/soft_total_abundance", sum(softabun))
    canonical("categorical/soft_occupied_cells",
              count(>(0), sum(softabun, dims = 1)))
    canonical("categorical/soft_species_surviving",
              count(>(0), sum(softabun, dims = 2)))
    canonical("categorical/soft_abundance_by_species",
              vec(sum(softabun, dims = 2)))

    canonical("categorical/hard_total_abundance", sum(hardabun))
    canonical("categorical/hard_occupied_cells",
              count(>(0), sum(hardabun, dims = 1)))
    canonical("categorical/hard_species_surviving",
              count(>(0), sum(hardabun, dims = 2)))
    canonical("categorical/hard_abundance_by_species",
              vec(sum(hardabun, dims = 2)))

    # And the share of each species' abundance in cells it tolerates - the quantity that says the
    # penalty is doing its job spatially rather than merely scaling everything down.
    canonical("categorical/soft_own_niche_share", own_niche_share(soft.eco))
    canonical("categorical/hard_own_niche_share", own_niche_share(hard.eco))

    # --- the properties, which no re-blessing may silence -------------------------------------
    @test all(>=(0), softabun)
    @test all(>=(0), hardabun)
    @test sum(softabun) > 0
    # The fixture is only meaningful while `survival` is non-zero - see the constructor.
    @test soft.eco.spplist.params.survival > 0

    # **The assertion this whole file exists for.** Soft and hard exclusion must remain
    # *distinguishable*: if the penalty stopped reaching the demographics - a wrong default, a
    # comparison inverted, `survival` set to zero - these two runs would converge, and every blessed
    # number above would still move in lockstep and could be re-blessed without anyone noticing.
    @test sum(hardabun) != sum(softabun)
    # Hard exclusion kills what disperses into an intolerable cell, so it must cost abundance.
    @test sum(hardabun) < sum(softabun)
    # ...and must concentrate what survives into each species' own niche.
    @test sum(own_niche_share(hard.eco)) > sum(own_niche_share(soft.eco))

    # Each run must have got the *same* niche map - the premise of every comparison here.
    @test soft.eco.habitat.regime.matrix == regimemap
    @test hard.eco.habitat.regime.matrix == regimemap

    # Reproducibility is the premise of every blessed number here, so it is asserted rather than
    # assumed: a fresh ecosystem on the same habitat and seed must reproduce the run exactly.
    again = categorical_run(0.5)
    @test again.eco.abundances.matrix == softabun
end

end
