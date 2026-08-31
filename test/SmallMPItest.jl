# SPDX-License-Identifier: LGPL-3.0-or-later

## Small scale test for MPI Ecosystems
using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using MPI
using Random
using Diversity
using JLD2
using Test

# The shared non-uniform, non-square, time-varying environment (see `test/varyingcase.jl`).
# This test compares results across 1/2/4 rank+thread splits — the strongest reproducibility
# property in the repository — and so must not run on a *uniform, square* grid, where every
# decomposition looks alike and a partitioning bug cannot show. The field it decomposes varies down
# `Y` (regime) and across `X` (supply), on a 7 × 12 grid that no rank count divides evenly.
include(joinpath(@__DIR__, "varyingcase.jl"))
# For `canonical_reference` only -- the blessed values are READ here, never written.
include(joinpath(@__DIR__, "canonical", "canonical.jl"))
using .Canonical

nt = Threads.nthreads();
@info "Total Memory: $(Sys.total_memory() / 2^30)GB, threads: $nt"

# Set up MPI and print threads
if !MPI.Initialized()
    MPI.Init()
end

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# **The fixture is built by `mpifixture_species` in `varyingcase.jl`, not spelled out here.**
# The canonical `mpi/…` results are blessed from a SERIAL run of that same builder, and those
# numbers are only evidence about this run if both sides build the identical thing — two
# spelled-out copies would drift, which is the failure the pinning exists to catch.
numSpecies = VARYING_SPECIES;
sppl, tolerance = mpifixture_species()

# **This is the ecosystem the cross-rank comparison actually uses** — it is simulated, gathered
# and saved below, while the second one further down only exercises the synchronise paths. It must
# decompose the shared varying field rather than one uniform temperature on a **square** grid: every
# cell being alike is precisely what stops a partitioning bug showing.
habitat = varying_environment()

# Set nichefit between species and environment (gaussian)
nichefit = NicheSuitability{Temperature, typeof(1.0K)}()

# build_ecosystem auto-selects the type from the live MPI session: >1 rank ⇒ MPIEcosystem, a single
# rank ⇒ serial Ecosystem (this script runs under mpiexec -n 1, 2 and 4). `sppl`/`habitat`/`nichefit` are
# exactly what MPIEcosystem takes directly.
expected = MPI.Comm_size(comm) > 1 ? MPIEcosystem : Ecosystem
@test build_ecosystem(sppl, habitat, nichefit = nichefit, seed = 0) isa expected
@test build_ecosystem(sppl, habitat, nichefit = nichefit, seed = 0,
                      distributed = false) isa Ecosystem
@test build_ecosystem(sppl, habitat, nichefit = nichefit, seed = 0,
                      distributed = true) isa MPIEcosystem

# Create ecosystem
@test_nowarn MPIEcosystem(sppl, habitat, nichefit)
eco = MPIEcosystem(sppl, habitat, nichefit, seed = 0)

# Artifically fill ecosystem with individuals
eco.abundances.rows_matrix .= 10

# **Compared GLOBALLY, not per rank, and that distinction is the whole point of an uneven split.**
# `synchronise_from_rows!` moves data from the species-partitioned layout to the cell-partitioned
# one, so the two hold the same values *in total* — but a single rank's share of each is only the
# same size when both partitions divide evenly. With 7 species over 77 cells on 2 ranks, rank 0 holds
# 4 × 77 = 3080 by rows and 7 × 39 = 2730 by columns; both are correct.
#
# The per-rank form must not be used: on a fixture where every rank gets an identical share — 8
# species on a 4 × 4 grid, say — it passes by asserting an artefact of that choice, not the
# invariant.
allsum(x) = MPI.Allreduce(sum(x), +, comm)
expected = numSpecies * prod(size(eco.habitat.regime.matrix)) * 10

# Set columns vector to zero and check synchronise from rows
eco.abundances.cols_vector .= 0
@test_nowarn EcoSISTEM.synchronise_from_rows!(eco.abundances)
@test allsum(eco.abundances.cols_vector) == allsum(eco.abundances.rows_matrix)
@test allsum(eco.abundances.cols_vector) == expected

# Set rows matrix to zero and check synchronise from cols
eco.abundances.rows_matrix .= 0
@test_nowarn EcoSISTEM.synchronise_from_cols!(eco.abundances)
@test allsum(eco.abundances.cols_vector) == allsum(eco.abundances.rows_matrix)
@test allsum(eco.abundances.rows_matrix) == expected

## Reproducibility is provided by the per-species RNG streams seeded in the
## MPIEcosystem constructor (seed = 0 above), independent of the process/thread
## split; the global RNG is not used by the simulation.

# Simulation Parameters
burnin = 2year;
times = 10year;
timestep = 1month_mean_duration;
record_interval = 3month_mean_duration;
repeats = 1;
lensim = length((0year):record_interval:times)

# Burnin
MPI.Barrier(comm)
@test sum(getabundance(eco)) ≈ 1.0
@test sum(getmetaabundance(eco)) ≈ 1.0
@test_nowarn simulate!(eco, burnin, timestep)
@test sum(getabundance(eco)) ≈ 1.0
@test sum(getmetaabundance(eco)) ≈ 1.0

# Collect full abundance matrix together
true_abuns = gatherabundance(eco)
# On root node, print abundances and save out
if rank == 0
    @save joinpath(ARGS[1], "Test_abuns$nt.jld2") abuns=true_abuns
end

sppl, tolerance = mpifixture_species()

# The shared varying environment: a temperature gradient down the grid, a solar gradient across it,
# and a steady warming over the run. Non-square (7 × 12) on purpose.
habitat = varying_environment()

# Set nichefit between species and environment (gaussian)
nichefit = NicheSuitability{Temperature, typeof(1.0K)}()

# Create ecosystem
@test_nowarn MPIEcosystem(sppl, habitat, nichefit)
eco = MPIEcosystem(sppl, habitat, nichefit, seed = 0)

# Artifically fill ecosystem with individuals
eco.abundances.rows_matrix .= 10
sleep(rank)

# Global again — the second ecosystem's split is just as uneven as the first's (see the note above).
# Set columns vector to zero and check synchronise from rows
eco.abundances.cols_vector .= 0
@test_nowarn EcoSISTEM.synchronise_from_rows!(eco.abundances)
@test allsum(eco.abundances.cols_vector) == allsum(eco.abundances.rows_matrix)
@test allsum(eco.abundances.cols_vector) == expected

# Set rows matrix to zero and check synchronise from cols
eco.abundances.rows_matrix .= 0
@test_nowarn EcoSISTEM.synchronise_from_cols!(eco.abundances)
@test allsum(eco.abundances.cols_vector) == allsum(eco.abundances.rows_matrix)
@test allsum(eco.abundances.rows_matrix) == expected

## Reproducibility is provided by the per-species RNG streams seeded in the
## MPIEcosystem constructor (seed = 0 above), independent of the process/thread
## split; the global RNG is not used by the simulation.

# Simulation Parameters
burnin = 2year;
times = 10year;
timestep = 1month_mean_duration;
record_interval = 3month_mean_duration;
repeats = 1;
lensim = length((0year):record_interval:times)

# Burnin
MPI.Barrier(comm)
@test_nowarn simulate!(eco, burnin, timestep)

sleep(rank)

# Collect full abundance matrix together
true_abuns = gatherabundance(eco)
# On root node, print abundances and save out
if rank == 0
    @save joinpath(ARGS[1], "Test_abuns$nt.jld2") abuns=true_abuns
end

# **Serial is the reference the distributed loop has to reproduce, and nothing else checks it.**
# The 1/2/4-rank comparison in `ext_EcoSISTEMMPIExt.jl` is MPI against MPI: every one of those runs
# takes the duplicated hot loop in `ext/EcoSISTEMMPIExt/dynamics.jl`, so a change made there and not
# in `src/dynamics.jl` agrees with itself at every rank count and passes. That is exactly how the
# birth draw stayed pre-2021 in the distributed loop while the serial one was fixed -- the serial
# loop draws `Poisson(n * rate)`, the MPI loop drew `Poisson(n * (1 - exp(-rate)))`, which also
# breaks the timestep independence the model requires.
#
# Run at one rank only: the partition is trivial there, so a difference is the loop body alone
# rather than anything about how the work was divided, and the serial answer cannot depend on the
# rank count anyway. Everything is rebuilt rather than reused -- an `Ecosystem` shares the habitat
# it is built on, so simulating a second one on the same object would not be independent.
if MPI.Comm_size(comm) == 1
    serialeco = mpifixture_ecosystem()
    simulate!(serialeco, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP)
    @test serialeco.abundances.matrix == true_abuns
end

# **And pin BOTH paths to one blessed number, at EVERY rank count.** The comparison above needs a
# serial run alongside, so it can only happen at one rank; these keys are recorded once by
# `test/canonical/test_mpifixture.jl` from a serial run and asserted here however many ranks there
# are. That is what makes a distributed-only divergence visible at 2 and 4 ranks.
#
# **Read-only on purpose.** `canonical(...)` would *write* the reference file, and doing that from
# inside an `mpiexec` child — several of them at once — must never happen. `canonical_reference()`
# only reads.
if rank == 0
    reference = Canonical.canonical_reference()
    grid = reshape(true_abuns, numSpecies, VARYING_NY, VARYING_NX)
    for (key, value) in ("mpi/total_abundance" => sum(grid),
        "mpi/abundance_by_species" => vec(sum(grid, dims = (2, 3))),
        "mpi/abundance_by_row" => vec(sum(grid, dims = (1, 3))),
        "mpi/abundance_by_column" => vec(sum(grid, dims = (1, 2))))
        # A missing key means the canonical set has not been blessed here; say so rather than
        # passing silently, which would make this whole check vacuous.
        @test haskey(reference, key)
        haskey(reference, key) &&
            @test isapprox(float.(collect(value)), reference[key],
                           rtol = 1e-8)
    end
end

if !MPI.Finalized()
    MPI.Finalize()
end
