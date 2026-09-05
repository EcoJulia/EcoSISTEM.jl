# SPDX-License-Identifier: LGPL-3.0-or-later

# Single-configuration benchmark worker.
#
# Builds one reasonably large model, warms up the simulation to pay JIT
# compilation, times a single `simulate!` run, and prints one machine-readable
# result line to stdout:
#
#     BENCH_RESULT,<mode>,<procs>,<threads>,<seconds>
#
# The model size and mode are read from the environment (all optional):
#
#     ECOSISTEM_BENCH_MODE         "threaded" (default) or "mpi"
#     ECOSISTEM_BENCH_GRID         square grid side length          (default 100)
#     ECOSISTEM_BENCH_SPECIES      number of species                (default 10000)
#     ECOSISTEM_BENCH_INDIVIDUALS  starting total individuals       (default 1e9)
#     ECOSISTEM_BENCH_YEARS        simulated years in the timed run (default 3)
#     ECOSISTEM_BENCH_SEED         RNG seed for reproducibility     (default 1234)
#     ECOSISTEM_BENCH_COUNT_TYPE   expected EcoSISTEM.Count type    (default Int64)
#
# ECOSISTEM_BENCH_COUNT_TYPE is only an assertion that the worker compiled with the
# expected count type; the type itself is chosen by run_benchmarks.jl through the
# EcoSISTEM `count_type` preference (a compile-time setting, not an env var).
# `EcoSISTEM.Count` exists only on the unmerged `rr/counttype` branch, so on any
# branch without it the assertion is skipped rather than failed - otherwise this
# worker cannot run at all there, which is what it did before this was guarded.
#
# The individual count defaults to ~1 billion so the simulation is compute-bound
# from the very first timestep (as the smaller default only becomes after tens of
# years of growth), which gives a stable speed estimate from a short run.
#
# This is normally launched by `run_benchmarks.jl`, not run directly, but it can
# be run on its own for a single measurement.
#
# The model this builds is unchanged from the one the deprecated `simplehabitat`/`GaussTrait`/
# `Gauss` construction produced, so timings stay comparable with earlier runs. Verified rather than
# assumed: same grid, cell size, per-cell supply and active mask; identical tolerance distributions,
# demand, movement kernels and growth parameters; and a bit-identical abundance matrix after twelve
# seeded timesteps.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols

# `EcoSISTEM.Count` is a compile-time type set through the `count_type` preference
# by the orchestrator. Fail loudly if this worker did not actually compile with
# the requested type (e.g. a stale precompile cache), so a benchmark can never
# silently measure the wrong storage type.
#
# Guarded by `isdefined`: the configurable count type is a feature of the unmerged
# `rr/counttype` branch, and an unguarded `EcoSISTEM.Count` is an outright
# `UndefVarError` on every branch without it. Checking rather than assuming keeps one
# worker correct on both - it still fails loudly where the type exists, and simply has
# nothing to check where it does not.
if isdefined(EcoSISTEM, :Count)
    const EXPECTED_COUNT = get(ENV, "ECOSISTEM_BENCH_COUNT_TYPE", "Int64")
    string(EcoSISTEM.Count) == EXPECTED_COUNT ||
        error("benchmark worker compiled with EcoSISTEM.Count=$(EcoSISTEM.Count) " *
              "but expected $EXPECTED_COUNT")
end

const MODE = get(ENV, "ECOSISTEM_BENCH_MODE", "threaded")
const GRID = parse(Int, get(ENV, "ECOSISTEM_BENCH_GRID", "100"))
const NUM_SPECIES = parse(Int, get(ENV, "ECOSISTEM_BENCH_SPECIES", "10000"))
const INDIVIDUALS = parse(Int,
                          get(ENV, "ECOSISTEM_BENCH_INDIVIDUALS",
                              string(1_000_000_000)))
const YEARS = parse(Int, get(ENV, "ECOSISTEM_BENCH_YEARS", "3"))
const SEED = parse(Int, get(ENV, "ECOSISTEM_BENCH_SEED", "1234"))

# Build the model: a single solar supply over a flat temperature grid, Gaussian thermal tolerance
# and birth-driven dispersal on a torus. The construction is identical for the plain and MPI
# ecosystems - `distributed` alone decides the container - so only that differs between modes.
#
# The grid is stated as an extent and a cell size, which is what a `StudyArea` decides a grid from:
# `AREA` is the total, so each side is its square root and each cell that divided by `GRID`.
const AREA = 1_000_000.0km^2
const SIDE = sqrt(AREA)

function build_model(mode::AbstractString)
    # `verbosity = :silent`: this prints one machine-readable result line, and a benchmark has no
    # use for the grid-decision commentary.
    area = StudyArea(extent = (SIDE, SIDE), cellsize = SIDE / GRID,
                     verbosity = :silent)
    environment = GridHabitat(regime = UniformSpec(274.0K,
                                                   axis = Temperature),
                              supply = UniformSpec(1000.0kJ / km^2 / day,
                                                   axis = SolarRadiation),
                              area = area)

    # `abundance` is given as an explicit **even** split rather than the integer total, which
    # `build_species` would instead divide up at random: every species doing the same work is what
    # makes a short run a stable speed estimate.
    species = build_species(NUM_SPECIES,
                            tolerance = (274.0K, 0.5K),
                            toleranceaxis = Temperature,
                            demand = 1.0kJ / day, demandaxis = SolarRadiation,
                            dispersal = 1.0km, pthresh = 10e-10,
                            movement = BirthOnlyMovement,
                            birth = 0.6 / year, death = 0.6 / year,
                            longevity = 1.0, survival = 0.2, boost = 1.0,
                            abundance = fill(div(INDIVIDUALS, NUM_SPECIES),
                                             NUM_SPECIES),
                            native = true)

    # `distributed` is forced rather than left at `:auto`, which would build a serial `Ecosystem`
    # on a single rank and so quietly benchmark the wrong container.
    return build_ecosystem(species, environment, seed = SEED,
                           distributed = mode == "mpi")
end

const TIMESTEP = 1month_mean_duration

if MODE == "mpi"
    using MPI
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    procs = MPI.Comm_size(comm)

    eco = build_model("mpi")

    # Warm-up run to trigger compilation, then the timed run.
    simulate!(eco, 2month_mean_duration, TIMESTEP)
    MPI.Barrier(comm)
    start = time_ns()
    simulate!(eco, YEARS * year, TIMESTEP)
    MPI.Barrier(comm)
    elapsed = (time_ns() - start) / 1e9

    if rank == 0
        println("BENCH_RESULT,mpi,$procs,$(Threads.nthreads()),$elapsed")
    end
    MPI.Finalize()
else
    eco = build_model("threaded")

    # Warm-up run to trigger compilation, then the timed run.
    simulate!(eco, 2month_mean_duration, TIMESTEP)
    start = time_ns()
    simulate!(eco, YEARS * year, TIMESTEP)
    elapsed = (time_ns() - start) / 1e9

    println("BENCH_RESULT,threaded,1,$(Threads.nthreads()),$elapsed")
end
