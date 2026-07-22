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
#
# The individual count defaults to ~1 billion so the simulation is compute-bound
# from the very first timestep (as the smaller default only becomes after tens of
# years of growth), which gives a stable speed estimate from a short run.
#
# This is normally launched by `run_benchmarks.jl`, not run directly, but it can
# be run on its own for a single measurement.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols

# `EcoSISTEM.Count` is a compile-time type set through the `count_type` preference
# by the orchestrator. Fail loudly if this worker did not actually compile with
# the requested type (e.g. a stale precompile cache), so a benchmark can never
# silently measure the wrong storage type.
const EXPECTED_COUNT = get(ENV, "ECOSISTEM_BENCH_COUNT_TYPE", "Int64")
string(EcoSISTEM.Count) == EXPECTED_COUNT ||
    error("benchmark worker compiled with EcoSISTEM.Count=$(EcoSISTEM.Count) " *
          "but expected $EXPECTED_COUNT")

const MODE = get(ENV, "ECOSISTEM_BENCH_MODE", "threaded")
const GRID = parse(Int, get(ENV, "ECOSISTEM_BENCH_GRID", "100"))
const NUM_SPECIES = parse(Int, get(ENV, "ECOSISTEM_BENCH_SPECIES", "10000"))
const INDIVIDUALS = parse(Int,
                          get(ENV, "ECOSISTEM_BENCH_INDIVIDUALS",
                              string(1_000_000_000)))
const YEARS = parse(Int, get(ENV, "ECOSISTEM_BENCH_YEARS", "3"))
const SEED = parse(Int, get(ENV, "ECOSISTEM_BENCH_SEED", "1234"))

# Build the model. Follows the pattern in examples/HPC/MPIRun.jl: a single solar
# supply with a flat temperature grid, Gaussian thermal tolerance and birth-driven
# dispersal on a torus. This keeps the construction identical for the plain and
# MPI ecosystems so only the container type differs between modes.
function build_ecosystem(mode::AbstractString)
    grid = (GRID, GRID)
    demand = 1.0kJ
    individuals = INDIVIDUALS
    area = 1_000_000.0km^2
    totalK = 1000.0kJ / km^2

    resource_vec = SolarDemand(fill(demand, NUM_SPECIES))

    birth = 0.6 / year
    death = 0.6 / year
    longevity = 1.0
    survival = 0.2
    boost = 1.0
    param = EqualPop(birth, death, longevity, survival, boost)

    kernel = fill(GaussianKernel(1.0km, 10e-10), NUM_SPECIES)
    movement = BirthOnlyMovement(kernel, Torus())

    opts = fill(274.0K, NUM_SPECIES)
    vars = fill(0.5K, NUM_SPECIES)
    tolerance = GaussTrait(opts, vars)
    native = fill(true, NUM_SPECIES)
    abun = fill(div(individuals, NUM_SPECIES), NUM_SPECIES)
    sppl = SpeciesList(NUM_SPECIES, tolerance, abun, resource_vec,
                       movement, param, native)

    habitat = simplehabitatAE(274.0K, grid, totalK, area)
    nichefit = Gauss{typeof(1.0K)}()

    if mode == "mpi"
        return MPIEcosystem(sppl, habitat, nichefit; seed = SEED)
    else
        return Ecosystem(sppl, habitat, nichefit; seed = SEED)
    end
end

const TIMESTEP = 1month

if MODE == "mpi"
    using MPI
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    procs = MPI.Comm_size(comm)

    eco = build_ecosystem("mpi")

    # Warm-up run to trigger compilation, then the timed run.
    simulate!(eco, 2months, TIMESTEP)
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
    eco = build_ecosystem("threaded")

    # Warm-up run to trigger compilation, then the timed run.
    simulate!(eco, 2months, TIMESTEP)
    start = time_ns()
    simulate!(eco, YEARS * year, TIMESTEP)
    elapsed = (time_ns() - start) / 1e9

    println("BENCH_RESULT,threaded,1,$(Threads.nthreads()),$elapsed")
end
