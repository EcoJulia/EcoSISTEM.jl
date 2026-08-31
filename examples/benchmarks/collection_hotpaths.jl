# SPDX-License-Identifier: LGPL-3.0-or-later

# Hot-path benchmark for multi-layer environments.
#
# Times the three paths whose cost depends on how many layers an environment has — `_suitability`
# (per cell, per species), `_resourceadjustment` (likewise) and `update_resource_usage!` (per cell)
# — plus one end-to-end `simulate!` for context, across four environment shapes:
#
#     1r1s   bare regime,       bare supply         the commonest case; no collection at all
#     2r1s   regime collection, bare supply         `_suitability` folds over the layers
#     1r2s   bare regime,       supply collection   `_resourceadjustment`/`update_resource_usage!` fold
#     2r2s   both                                   everything at once
#
# The point of the 1r1s row is that a bare layer answers the same interface as a collection
# (`values(layer)` is a one-tuple), so the commonest case now goes through the same fold as the rest
# and must not have paid for it.
#
# Run it directly; nothing here downloads anything:
#
#     julia --project=examples -t 1 examples/benchmarks/collection_hotpaths.jl
#
# Size and repetition count are read from the environment (all optional):
#
#     ECOSISTEM_BENCH_SPECIES   number of species                     (default 100)
#     ECOSISTEM_BENCH_GRID      square grid side length               (default 60)
#     ECOSISTEM_BENCH_STEPS     timesteps in the `simulate!` run      (default 5)
#     ECOSISTEM_BENCH_REPEATS   timed repetitions, minimum reported   (default 31)
#     ECOSISTEM_BENCH_LABEL     text printed in the header            (default "EcoSISTEM")
#
# ## Comparing two versions of the package
#
# The last column is a hash of the final abundance matrix. Every run is seeded, so two versions of
# EcoSISTEM that compute the same thing print the same hash — which is what makes a *speed*
# comparison meaningful rather than a comparison of two different simulations. To measure a change
# against the code it replaces, extract the committed (or staged) tree somewhere and point a second
# run at it, with the same manifest so both resolve identical dependency versions:
#
#     git checkout-index -a -f --prefix=/tmp/before/     # or `git worktree add` for a commit
#     cp Manifest.toml /tmp/before/
#     julia --project=/tmp/before -t 1 examples/benchmarks/collection_hotpaths.jl
#
# `checkout-index` only reads the index, so it cannot disturb the working tree. Run single-threaded
# (`-t 1`) unless threading is what is being measured: the species loop in `update!` is
# `Threads.@threads`, and the scheduling noise swamps the millisecond differences this is looking for.
#
# This script is written against the current API, so it will not run unaltered against a tree
# older than the one-type-per-family collections (v0.5.0). Four constructors need swapping for their
# fixed-arity predecessors there, and nothing else: `LayerCollection((a, b))` →
# `RegimeCollection2(a, b)` for a regime and `SupplyCollection2(a, b)` for a supply,
# `SpeciesRequirementCollection((a, b))` → `ToleranceCollection2(a, b)`, `SpeciesRequirementCollection((a, b))` →
# `DemandCollection2(a, b)`, and `MultiplicativeFit((a, b))` → `MultiplicativeFit2(a, b)`. That
# variant is what produced the before/after figures recorded in
# `~/.claude/plans/ecosistem-collection-types.md`.
#
# **Do not lift the parameters from `test/TestCases.jl`.** Its fixtures use `survival = 0.0`,
# which makes the suitability term `ϵ̄real^-0.0 == 1` — so the regime, tolerance and nichefit have no
# effect on the dynamics whatsoever, every shape below runs the identical trajectory, and the
# population explodes past 1e11 within a handful of steps. The values here are `build_species`' own
# defaults, which keep the run stable and the shapes genuinely different from one another.

module CollectionHotpaths

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using Random
using Printf

const NSPP = parse(Int, get(ENV, "ECOSISTEM_BENCH_SPECIES", "100"))
const SIDE = parse(Int, get(ENV, "ECOSISTEM_BENCH_GRID", "60"))
const STEPS = parse(Int, get(ENV, "ECOSISTEM_BENCH_STEPS", "5"))
const REPEATS = parse(Int, get(ENV, "ECOSISTEM_BENCH_REPEATS", "31"))
const LABEL = get(ENV, "ECOSISTEM_BENCH_LABEL", "EcoSISTEM")

const GRID = (SIDE, SIDE)
const CELLSIZE = 1.0km
const SEED = 1234

# A temperature regime and a solar supply, or a rainfall regime and a water supply, filled flat.
# The values are arbitrary; only the shapes matter.
#
# These were once built from the layer constructors directly, on the grounds that
# `GridHabitat` "needs a positioned `StudyArea`, and so a data source, for more than one
# layer". Measured 2026-08-07: that is not so — it takes a tuple on either side over a synthetic
# area. Going through it means the benchmark measures the shape users actually get.
#
# One real difference, and it is why this is worth stating: a bare `ContinuousRegime(…)` declares
# no axis, so building one directly gives `EcoSISTEM.NicheAxis` where `GridHabitat` gives it the
# axis of its spec (`Temperature`/`Precipitation`). That parameter is dispatched on, so a benchmark
# built the first way measures a shape no real environment has. Both arms are built through the
# habitat, so the 1-vs-2 comparison these benchmarks exist for is unaffected.
function _area()
    return StudyArea(extent = (GRID[1] * CELLSIZE, GRID[2] * CELLSIZE),
                     cellsize = CELLSIZE, verbosity = :silent)
end
_tempspec() = UniformSpec(10.0K, axis = Temperature)
_rainspec() = UniformSpec(10.0mm / day, axis = Precipitation)
_solarspec() = UniformSpec(1.0e4kJ / (km^2 * day), axis = SolarRadiation)
_waterspec() = UniformSpec(1.0e4mm / day, axis = Precipitation)

# One ecosystem per shape: `nregimes`/`nsupplies` of 1 give bare layers, 2 give collections.
function build(nregimes::Int, nsupplies::Int)
    regime = nregimes == 1 ? _tempspec() : (_tempspec(), _rainspec())
    supply = nsupplies == 1 ? _solarspec() : (_solarspec(), _waterspec())
    habitat = GridHabitat(regime = regime, supply = supply,
                          area = _area())

    temptol = NicheTolerance(Temperature, Normal, fill(10.0K, NSPP),
                             fill(2.0K, NSPP))
    raintol = NicheTolerance(Precipitation, Normal, fill(10.0mm / day, NSPP),
                             fill(2.0mm / day, NSPP))
    tolerance = nregimes == 1 ? temptol :
                SpeciesRequirementCollection((temptol, raintol))
    nichefit = nregimes == 1 ? NicheSuitability(temptol) :
               MultiplicativeFit((NicheSuitability(temptol),
                                  NicheSuitability(raintol)))

    solar = Demand{SolarRadiation}(fill(2.0kJ / day, NSPP))
    water = Demand{Precipitation}(fill(2.0Unitful.L / day, NSPP))
    demand = nsupplies == 1 ? solar :
             SpeciesRequirementCollection((solar, water))

    Random.seed!(SEED)
    abun = rand(Multinomial(2000 * NSPP, NSPP))
    movement = BirthOnlyMovement(GaussianKernel.(fill(1.0km, NSPP), 1.0e-3))
    # `build_species`' own defaults — see the warning in the header about `survival = 0.0`.
    param = EqualPop(0.6 / year, 0.6 / year, 1.0, 0.2, 1.0)
    species = SpeciesList(NSPP, tolerance, abun, demand, movement, param,
                          fill(true, NSPP))
    return Ecosystem(species, habitat, nichefit, seed = SEED)
end

# The fastest of `REPEATS` timings, after one untimed call to pay for compilation. The minimum
# rather than the mean: every source of noise here only ever makes a run slower.
function fastest(f)
    f()
    return minimum(@elapsed(f()) for _ in 1:REPEATS)
end

# Every cell × species evaluation of a path, which is the work one `update!` timestep does.
function sweep(f, eco)
    total = 0.0
    for sp in Base.axes(eco.abundances.matrix, 1),
        i in Base.axes(eco.abundances.matrix, 2)
        total += f(eco, i, sp)
    end
    return total
end

_suit(eco, i, sp) = EcoSISTEM.suitability(eco, i, sp)
function _adjust(eco, i, sp)
    birth, death = EcoSISTEM._resourceadjustment(eco, eco.habitat.supply, i,
                                                 sp)
    return birth + death
end

# `update_resource_usage!` returns early on a valid cache, so it has to be invalidated each time.
function _totaldemand(eco)
    eco.cache.valid = false
    return EcoSISTEM.update_resource_usage!(eco)
end

function run()
    @printf("# %s — %d species, %d×%d cells, %d thread(s), min of %d\n", LABEL,
            NSPP, SIDE, SIDE, Threads.nthreads(), REPEATS)
    @printf("%-6s %13s %13s %13s %13s %20s\n", "shape", "simulate!(s)",
            "suitability(ms)", "adjust(ms)", "totaldemand(ms)",
            "abundance hash")
    for (nregimes, nsupplies) in ((1, 1), (2, 1), (1, 2), (2, 2))
        eco = build(nregimes, nsupplies)
        sweep(_suit, eco)                      # compile every path before timing any of it
        sweep(_adjust, eco)
        _totaldemand(eco)

        # A separate ecosystem for the hash: `simulate!` mutates, and the timed runs each need a
        # fresh one anyway, so this one is simulated once and its final state hashed.
        checked = build(nregimes, nsupplies)
        simulate!(checked, STEPS * 1month_mean_duration, 1month_mean_duration)

        simulated = fastest() do
            return simulate!(build(nregimes, nsupplies),
                             STEPS * 1month_mean_duration, 1month_mean_duration)
        end
        @printf("%-6s %13.4f %13.4f %13.4f %13.4f %20s\n",
                string(nregimes, "r", nsupplies, "s"), simulated,
                1.0e3 * fastest(() -> sweep(_suit, eco)),
                1.0e3 * fastest(() -> sweep(_adjust, eco)),
                1.0e3 * fastest(() -> _totaldemand(eco)),
                string(hash(checked.abundances.matrix), base = 16))
    end
    return nothing
end

end

CollectionHotpaths.run()
