# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A deliberately **non-uniform, non-square, time-varying** ecosystem fixture, shared by the canonical
# blessed-result tests and the MPI cross-rank reproducibility test.
#
# **Why it exists.** A reproducible-results test - the canonical simulated run and the MPI
# 1/2/4-rank comparison alike - must not run on a *spatially uniform, temporally static, square*
# environment. Such a fixture is structurally incapable of detecting:
#
#   - anything **spatial**, including the `(y, x)` dimension-order class that only a **non-square**
#     grid exposes (see CLAUDE.md, which records several of exactly that shape);
#   - anything **temporal**: the whole layer-change subsystem (`SteadyLayerChange`,
#     `PatternedLayerChange`, `SeriesLayerChange`, `atend`, calendars) had no blessed coverage at all;
#   - a **partitioning** bug in MPI, since with identical cells every decomposition looks alike.
#
# This fixture is the counterpart to `TestCases.jl`'s uniform ones, not a replacement: the uniform
# blessed numbers stay exactly as they are, so the historical baseline is preserved and these are
# purely additional keys.
#
# Nothing here may depend on thread or rank count. One deterministic RNG stream per species
# (`seed`) is what makes that true, and it is the premise the MPI comparison tests.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Distributions
using Random

# **Not square, and - separately - not evenly divisible by any rank count.** Three distinct
# properties, each guarding something different:
#
#   - `NX != NY` so a **transposed index** is a shape error rather than a plausible wrong number;
#   - **77 cells** (11 × 7) so the *cell* partition is uneven at every rank count: 4 ranks take
#     20/19/19/19 and 2 ranks take 39/38, rather than dividing exactly;
#   - **7 species** so the *species* partition is uneven too: 2/2/2/1 over 4 ranks, 4/3 over 2.
#
# `MPIEcosystem` partitions by species (rows) **and** grid cells (columns), so both halves need a
# remainder to put the decomposition under any stress at all. An earlier version of this fixture
# used 12 × 7 = 84 cells and 8 species - both of which divide exactly by 2 and 4, so every rank got an
# identical share and the uneven-split paths were never taken.
const VARYING_NX = 11
const VARYING_NY = 7
const VARYING_CELL = 1.0km
const VARYING_SPECIES = 7

# The study area for the fixture: a rectangle, so `size(regime.matrix)` is asymmetric and any
# transposition shows up as a shape error rather than as a silently wrong number.
function varying_area()
    return StudyArea(extent = (VARYING_NY * VARYING_CELL,
                               VARYING_NX * VARYING_CELL),
                     cellsize = VARYING_CELL, verbosity = :silent)
end

# A regime that varies in **space** (a north-south temperature gradient) and in **time** (a steady
# warming). Both at once is the point: a change that got the spatial orientation right and the clock
# wrong, or vice versa, moves these numbers.
#
# `IncrementBy` is a *rate*, so the total warming over the run is a function of elapsed time - which
# is exactly the timestep-independence the model requires, and would break if a change were ever
# applied per *step* rather than per unit time.
function varying_regime()
    return Varying(GradientSpec(288.0K, 302.0K, axis = Temperature,
                                orientation = 180°),
                   IncrementBy(0.5K / year))
end

# Two supplies, one of them **also spatially varying**, so the resource ratio differs cell by cell
# and the birth/death arithmetic is exercised against a gradient rather than one constant.
#
# **The supply gradient runs across the grid (`orientation = 90°`) while the regime's runs down it
# (`180°`), and that orthogonality is the whole design.** Measured: the regime takes 7 distinct
# values (one per row, `Y`) and the supply 11 (one per column, `X`). A transposed index therefore
# swaps those counts - 11 temperatures and 7 supplies - which is loud, rather than producing a
# plausible wrong number. On the old square, uniform grid neither could be seen at all.
function varying_supply()
    return (GradientSpec(3.0e11kJ / km^2 / day, 6.0e11kJ / km^2 / day,
                         axis = SolarRadiation, orientation = 90°),
            UniformSpec(150.0mm / day, axis = Precipitation))
end

# The assembled environment on its own, so the MPI cross-rank test can decompose exactly the same
# non-uniform field the canonical runs bless. Sharing it is the point: a partitioning bug is
# invisible when every cell is identical.
function varying_environment()
    return GridHabitat(regime = varying_regime(),
                       supply = varying_supply(),
                       area = varying_area(),
                       topology = Torus())
end

# Species optima spread across the regime's own gradient, so species genuinely sort in space rather
# than all preferring the same cells. Exposed separately because any ecosystem built on this
# environment - serial or MPI - needs optima that actually sit inside it.
function varying_optima(numspecies = VARYING_SPECIES)
    return 288.0K .+
           (302.0 - 288.0)K .* range(0.1, stop = 0.9,
                                     length = numspecies)
end

function varying_species(; numspecies = VARYING_SPECIES, seed = 0)
    vars = fill(2.0K, numspecies)
    opts = varying_optima(numspecies)
    demand = SpeciesRequirementCollection((Demand{SolarRadiation}(fill(4.0e5kJ /
                                                                       m^2 /
                                                                       day *
                                                                       1.0m^2,
                                                                       numspecies)),
                                           Demand{Precipitation}(fill(150.0Unitful.L /
                                                                      m^2 /
                                                                      day *
                                                                      1.0m^2,
                                                                      numspecies))))
    movement = BirthOnlyMovement(GaussianKernel.(fill(1.5km, numspecies),
                                                 10e-10))
    param = EqualPop(0.15 / year, 0.15 / year, 1.0, 0.1, 1.0)
    # Drawn from the *seeded* stream, never the global RNG, or the starting state would differ
    # between runs and nothing downstream could be blessed.
    rng = Random.Xoshiro(seed)
    abun = rand(rng, Multinomial(100_000, numspecies))
    tolerance = NicheTolerance(Temperature, Normal, opts, vars)
    return SpeciesList(numspecies, tolerance, abun, demand, movement, param,
                       fill(true, numspecies)), tolerance
end

"""
    varying_ecosystem(; seed = 0, numspecies = VARYING_SPECIES)

Build the shared non-uniform, non-square, time-varying ecosystem. Returns a plain `Ecosystem`; the
MPI test builds an `MPIEcosystem` from the same `SpeciesList`/environment instead.
"""
function varying_ecosystem(; seed = 0, numspecies = VARYING_SPECIES)
    env = varying_environment()
    sppl, tolerance = varying_species(numspecies = numspecies, seed = seed)
    return Ecosystem(sppl, env, NicheSuitability(tolerance), seed = seed)
end

# -- The fixture the serial-versus-distributed pinning uses --------------------------------------
#
# **Defined once here on purpose.** The blessed `mpi/...` results are only evidence about the
# distributed run if the canonical bless and `SmallMPItest.jl` build the *identical* fixture, and two
# spelled-out copies would drift apart - which is precisely the failure this pinning exists to catch.
#
# Deliberately unlike `varying_species` above: faster rates, a wider kernel and a large `boost`, so a
# two-year run moves the numbers substantially instead of sitting near its starting state.
#
# **Seven species, and that is the point.** `MPIEcosystem` partitions by species *and* by grid cells,
# so both need a remainder to exercise the uneven-split paths: 7 species go 2/2/2/1 over four ranks
# and 4/3 over two, while the 77 cells go 20/19/19/19 and 39/38.
const MPIFIXTURE_FILL = 10
const MPIFIXTURE_BURNIN = 2year
const MPIFIXTURE_TIMESTEP = 1month_mean_duration

# The dispersal kernel both MPI fixtures share. Named once so that `mpifixture_species` and
# `mpifixture_ecosystem` cannot drift apart in it: the serial run is the reference the distributed
# one has to reproduce, and a difference here would look exactly like a partitioning bug.
mpifixture_kernels(numspecies) = fill(GaussianKernel(3.0km, 10e-10), numspecies)
function mpifixture_movement(numspecies = VARYING_SPECIES)
    return BirthOnlyMovement(mpifixture_kernels(numspecies))
end
# The animal-like counterpart, dispersing the standing population rather than only the newborns.
# It reaches the distributed loop by a different route - `EcoSISTEM._standingpopulation` rather than
# the `births` the loop already holds - so pinning only one of the two would leave that route unrun.
function mpifixture_always(numspecies = VARYING_SPECIES)
    return AlwaysMovement(mpifixture_kernels(numspecies))
end

function mpifixture_species(; numspecies = VARYING_SPECIES,
                            movement = mpifixture_movement(numspecies))
    # Two demands, matching the two supplies of the varying environment: `build_ecosystem` checks
    # that the species and environment sides align layer for layer.
    demand = SpeciesRequirementCollection((Demand{SolarRadiation}(collect(1:numspecies) .*
                                                                  1.0kJ /
                                                                  day),
                                           Demand{Precipitation}(fill(2.0Unitful.L /
                                                                      day,
                                                                      numspecies))))
    param = EqualPop(0.6 / year, 0.6 / year, 1.0, 0.2, 100.0)
    # Optima spread across the regime's own 288-302 K gradient, so species genuinely sort in space.
    # A shared optimum would put them all in the same cells, and optima outside the environment would
    # kill everything -- a run where nothing survives compares equal across ranks for the wrong reason.
    tolerance = NicheTolerance(Temperature, Normal, varying_optima(numspecies),
                               fill(2.0K, numspecies))
    return SpeciesList(numspecies, tolerance,
                       fill(div(1_000, numspecies), numspecies), demand,
                       movement,
                       param, fill(true, numspecies)), tolerance
end

mpifixture_nichefit() = NicheSuitability{Temperature, typeof(1.0K)}()

"""
    mpifixture_ecosystem(; seed = 0, numspecies = VARYING_SPECIES)

Build the **serial** ecosystem that a distributed run of the same fixture must reproduce, with its
abundances filled exactly as the MPI test fills its own.

Simulating this for `MPIFIXTURE_BURNIN` at `MPIFIXTURE_TIMESTEP` and gathering the distributed run
must give the same abundance matrix, cell for cell, at any rank count.
"""
function mpifixture_ecosystem(; seed = 0, numspecies = VARYING_SPECIES,
                              movement = mpifixture_movement())
    sppl, _ = mpifixture_species(numspecies = numspecies, movement = movement)
    eco = Ecosystem(sppl, varying_environment(), mpifixture_nichefit(),
                    seed = seed)
    # Filled artificially rather than populated, so the starting state is identical on every rank
    # without depending on how `populate!` divided the draw.
    eco.abundances.matrix .= MPIFIXTURE_FILL
    return eco
end
