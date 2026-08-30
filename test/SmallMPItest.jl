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

nt = Threads.nthreads();
@info "Total Memory: $(Sys.total_memory() / 2^30)GB, threads: $nt"

# Set up MPI and print threads
if !MPI.Initialized()
    MPI.Init()
end

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Set up initial parameters for ecosystem
# **7, not 8, and that is the point of this test.** `MPIEcosystem` partitions by species *and* by
# grid cells, so both need a remainder to exercise the uneven-split paths: 7 species go 2/2/2/1 over
# four ranks and 4/3 over two, while the shared fixture's 77 cells go 20/19/19/19 and 39/38. With 8
# species on 84 cells every rank got an identical share and neither split was ever unbalanced.
numSpecies = VARYING_SPECIES;
demand = 10.0kJ / day;
individuals = 1_000;
area = 100.0 * km^2;
totalK = 10000.0kJ / km^2 / day;

# Set up how much resource each species consumes.
# **Two layers, matching the two supplies of the varying environment.** `build_ecosystem` checks
# that species demand and environment supply align layer for layer, so a solar-only demand against a
# solar+water environment is refused — correctly. Defined once here and reused by both species lists
# below, rather than being rebuilt differently half way down the file.
resource_vec = Demand{SolarRadiation}(collect(1:numSpecies) .* 1.0kJ / day)
water_vec = Demand{Precipitation}(fill(2.0Unitful.L / day, numSpecies))
total_use = SpeciesRequirementCollection((resource_vec, water_vec))
# Set probabilities
birth = 0.6 / year
death = 0.6 / year
longevity = 1.0
survival = 0.2
boost = 100.0

# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(3.0km, 10e-10), numSpecies)
movement = BirthOnlyMovement(kernel)

# Create species list, including their temperature preferences, seed abundance and native status
# Spread across the varying regime's own 288–302 K gradient. A single shared optimum would put
# every species in the same cells; `fill(274.0K, …)` would put them all outside the environment
# entirely, and a run where everything dies compares equal across ranks for the wrong reason.
opts = varying_optima(numSpecies)
vars = fill(2.0K, numSpecies)
tolerance = NicheTolerance(Temperature, Normal, opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, tolerance, abun, total_use,
                   movement, param, native)

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

sppl = SpeciesList(numSpecies, tolerance, abun, total_use, movement, param,
                   native)

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
# takes the duplicated hot loop in `ext/EcoSISTEMMPIExt/generate.jl`, so a change made there and not
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
    serialsppl = SpeciesList(numSpecies, tolerance, abun, total_use, movement,
                             param, native)
    serialeco = Ecosystem(serialsppl, varying_environment(), nichefit, seed = 0)
    # Matched to the distributed run above, which fills its own abundances the same way.
    serialeco.abundances.matrix .= 10
    simulate!(serialeco, burnin, timestep)
    @test serialeco.abundances.matrix == true_abuns
end

if !MPI.Finalized()
    MPI.Finalize()
end
