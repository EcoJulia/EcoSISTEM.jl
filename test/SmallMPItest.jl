# SPDX-License-Identifier: LGPL-3.0-or-later

## Small scale test for MPI Ecosystems
using EcoSISTEM
using EcoSISTEM: materialise
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using MPI
using Random
using Diversity
using DataFrames: nrow
using Phylo: getbranches, getlength
using JLD2
using Test
using ArchGDAL
using Dates: Date
using DimensionalData: DimensionalData, dims, refdims
using Rasters

# The shared non-uniform, non-square, time-varying environment (see `test/varyingcase.jl`).
# This test compares results across 1/2/4 rank+thread splits - the strongest reproducibility
# property in the repository - and so must not run on a *uniform, square* grid, where every
# decomposition looks alike and a partitioning bug cannot show. The field it decomposes varies down
# `Y` (regime) and across `X` (supply), on a 7 × 12 grid that no rank count divides evenly.
include(joinpath(@__DIR__, "varyingcase.jl"))
include(joinpath(@__DIR__, "layercompare.jl"))
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

# **Two testsets below run at some rank counts only, and this file is launched at 1, 2 and 4.**
# Each is about what ranks do *to each other*, so at one rank there is nothing to see - and they are
# the two most expensive things here, measured at 22 s and 9 s of a 53 s two-rank run, paid again in
# every launch. What every launch still runs is everything a rank count can change: the partition,
# the hot loop, the interventions, the recorders and the blessed comparisons.
const MANYRANKS = MPI.Comm_size(comm) > 1
# The shared build is the stricter case: a third and fourth rank receive exactly what the second
# does, so two ranks prove it and more only pay for it again.
const TWORANKS = MPI.Comm_size(comm) == 2

# Every rank asks for the same download at once: one fetches it, the others wait on its lock and
# then find the file, so the bytes travel once and nothing is corrupted. A `file://` source, so no
# network; the shared directory is made on rank 0 and its name sent round.
if !MANYRANKS
    @info "Skipping the shared-asset testset: it needs more than one rank."
else
    @testset "one rank fetches an asset the others wait for" begin
        shared = MPI.bcast(rank == 0 ? mktempdir() : "", 0, comm)
        source = joinpath(shared, "source.bin")
        rank == 0 && write(source, rand(UInt8, 4096))
        MPI.Barrier(comm)
        dest = joinpath(shared, "fetched.bin")
        p = EcoSISTEM.assetpath(EcoSISTEM.CachedAsset(TwentyCR,
                                                      "file://" * source,
                                                      path = dest))
        MPI.Barrier(comm)
        @test p == dest && read(dest) == read(source)
        @test isfile(dest * ".provenance.toml")
        @test !isfile(dest * ".part") && !isfile(dest * ".lock")
    end
end

# **A study area and its layers are built once, on the root, and every other rank receives them.**
# The data is a three-band
# WGS84 GeoTIFF, written by the root into a directory every rank can see and read as a dated series
# with a gap in its second slice, so a series, its units, its `NaN`, its CRS and its dates all cross
# between ranks.
#
# The ranks other than the root are given paths that do not exist. Their builds can only succeed
# if they read nothing - by the cache or around it - and take the root's result instead.
if !TWORANKS
    @info "Skipping the shared study area and layers testset: it runs at two ranks."
else
    @testset "a study area and its layers are built on the root and received by every rank" begin
        shared = MPI.bcast(rank == 0 ? mktempdir() : "", 0, comm)
        real = joinpath(shared, "rain.tif")
        if rank == 0
            ArchGDAL.create(real, driver = ArchGDAL.getdriver("GTiff"),
                            width = 8,
                            height = 7, nbands = 3, dtype = Float64) do ds
                for band in 1:3
                    values = [2.0 + i + 10j + 100band for i in 1:8, j in 1:7]
                    band == 2 && (values[3, 2] = NaN)
                    ArchGDAL.write!(ds, values, band)
                end
                ArchGDAL.setgeotransform!(ds, [10.0, 0.5, 0.0, 55.0, 0.0, -0.5])
                return ArchGDAL.setproj!(ds,
                                         ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
            end
        end
        MPI.Barrier(comm)
        dates = [Date(2000, 1, 15), Date(2000, 2, 15), Date(2000, 3, 15)]
        rain(path) = RasterFileSpec(path, axis = Precipitation, unit = mm / day,
                                    times = dates, atend = HoldAtEnd())
        mine = rain(rank == 0 ? real : joinpath(shared, "absent.tif"))

        # Gathered by serialisation - a route independent of the one under test - and compared on the
        # root against what the root built.
        function agreeswithroot(same, value, reference)
            everyone = MPI.gather(value, comm)
            rank == 0 || return true
            return all(same(reference, other) for other in everyone)
        end
        agreeswithroot(layer, reference) = agreeswithroot(samelayer, layer,
                                                          reference)

        # A shape file only the root can see, as the area's mask: a triangle, so the grid cropped to it
        # still has inactive cells.
        outline = joinpath(shared, "outline.geojson")
        if rank == 0
            triangle = ArchGDAL.createpolygon([[(10.5, 51.8), (13.6, 51.8),
                                                  (10.5, 54.6), (10.5, 51.8)]])
            ArchGDAL.create(outline, driver = ArchGDAL.getdriver("GeoJSON")
                            ) do ds
                makelayer = layer -> begin
                    ArchGDAL.addfeature(layer) do feature
                        return ArchGDAL.setgeom!(feature, triangle)
                    end
                    ArchGDAL.copy(layer, dataset = ds)
                end
                return ArchGDAL.createlayer(makelayer, name = "outline",
                                            geom = ArchGDAL.wkbPolygon)
            end
        end
        MPI.Barrier(comm)
        within = ShapeSpec(rank == 0 ? outline :
                           joinpath(shared, "absent.geojson"))

        area = StudyArea(supply = mine, verbosity = :silent)
        @test agreeswithroot(samereport, area.report,
                             rank == 0 ? area.report : nothing)
        @test count(area.report.active) == 56
        # The area's own records arrive with it, before any layer is built.
        @test agreeswithroot(==, provenance(area).inputs,
                             rank == 0 ? provenance(area).inputs : nothing)
        @test length(provenance(area).inputs) == 1
        investigated = investigate_study_area(supply = mine)
        @test agreeswithroot(samereport, investigated,
                             rank == 0 ? investigated : nothing)
        masked = StudyArea(supply = mine, within = within, verbosity = :silent)
        @test agreeswithroot(samereport, masked.report,
                             rank == 0 ? masked.report : nothing)
        @test 0 < count(masked.report.active) < length(masked.report.active)
        synthetic = StudyArea(extent = (70.0km, 80.0km), cellsize = 10.0km,
                              verbosity = :silent)
        @test agreeswithroot(samereport, synthetic.report,
                             rank == 0 ? synthetic.report : nothing)

        direct = materialise(mine, area, role = EcoSISTEM.Resource)
        @test direct.change isa EcoSISTEM.SeriesLayerChange
        @test any(isnan, direct.change.slices)
        reference = rank == 0 ?
                    EcoSISTEM._materialisespec(rain(real), area,
                                               EcoSISTEM.Resource) : nothing
        @test agreeswithroot(direct, reference)

        niche = NicheSpec(4, axis = EcoSISTEM.NicheAxis)
        habitat = GridHabitat(regime = (temperature = UniformSpec(290.0K,
                                                                  axis = Temperature),
                                        niche = niche),
                              supply = mine, area = area)
        # The habitat's own supply has had its gaps zeroed, so it is compared with the root's.
        @test agreeswithroot(habitat.supply,
                             rank == 0 ? habitat.supply : nothing)
        # Unseeded, so each rank would draw its own layout if it built its own.
        @test agreeswithroot(habitat.regime[:niche],
                             rank == 0 ? habitat.regime[:niche] : nothing)
        everyactive = MPI.gather(habitat.active, comm)
        rank == 0 && @test all(==(first(everyactive)), everyactive)

        # Only the root read, and every rank records what it read.
        @test isempty(area.report.cache.reads) == (rank != 0)
        everyinput = MPI.gather(provenance(habitat).inputs, comm)
        rank == 0 && @test !isempty(provenance(habitat).inputs)
        rank == 0 && @test all(==(first(everyinput)), everyinput)

        # A layer already built is copied on every rank, with nothing sent.
        rebuilt = GridHabitat(regime = habitat.regime[:temperature],
                              supply = habitat.supply, area = area)
        @test samelayer(rebuilt.supply, habitat.supply)

        # A failure on the root stops every rank, rather than leaving the others waiting.
        broken = rain(rank == 0 ? joinpath(shared, "absent.tif") : real)
        @test_throws Exception materialise(broken, area)
        MPI.Barrier(comm)

        # So does a shared build started inside another.
        @test_throws "started inside another one" EcoSISTEM._sharedlayer(area) do
            return EcoSISTEM._sharedlayer(() -> nothing, area)
        end
        MPI.Barrier(comm)
    end
end

# **The fixture is built by `mpifixture_species` in `varyingcase.jl`, not spelled out here.**
# The canonical `mpi/...` results are blessed from a SERIAL run of that same builder, and those
# numbers are only evidence about this run if both sides build the identical thing - two
# spelled-out copies would drift, which is the failure the pinning exists to catch.
numSpecies = VARYING_SPECIES;
sppl, tolerance = mpifixture_species()

# **This is the ecosystem the cross-rank comparison actually uses** - it is simulated, gathered
# and saved below, while the second one further down only exercises the synchronise paths. It must
# decompose the shared varying field rather than one uniform temperature on a **square** grid: every
# cell being alike is precisely what stops a partitioning bug showing.
habitat = varying_environment()

# Set nichefit between species and environment (gaussian)
nichefit = NicheSuitability{Temperature, typeof(1.0K)}()

# build_ecosystem auto-selects the type from the live MPI session: >1 rank => MPIEcosystem, a single
# rank => serial Ecosystem (this script runs under mpiexec -n 1, 2 and 4). `sppl`/`habitat`/`nichefit` are
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
# one, so the two hold the same values *in total* - but a single rank's share of each is only the
# same size when both partitions divide evenly. With 7 species over 77 cells on 2 ranks, rank 0 holds
# 4 * 77 = 3080 by rows and 7 * 39 = 2730 by columns; both are correct.
#
# The per-rank form must not be used: on a fixture where every rank gets an identical share - 8
# species on a 4 × 4 grid, say - it passes by asserting an artefact of that choice, not the
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

# **These replace two assertions that could not fail.** `sum(getabundance(eco)) ≈ 1.0` and the same
# of `getmetaabundance` were each true of any block normalised by its own total, at any rank count -
# so they passed while the metaabundance vector was silently only this rank's species (7 at one rank,
# 4 at two, 2 at four) and the weights summed to the rank count rather than to 1.
#
# What is asserted instead: that a whole-metacommunity question is **refused** rather than answered
# from one rank's block, and that the two quantities which are genuinely global agree with the serial
# answer whatever the rank count.
diversityrefusal = "no single rank holds them"
MPI.Barrier(comm)
@test_throws diversityrefusal getabundance(eco)
@test length(getmetaabundance(eco)) == numSpecies
@test sum(getmetaabundance(eco)) ≈ 1.0
@test sum(getweight(eco)) ≈ 1.0
@test length(getweight(eco)) == VARYING_NY * VARYING_NX

@test_nowarn simulate!(eco, burnin, timestep)

# The same after a run, and now with the values pinned rather than only their totals: both are
# metacommunity quantities, so they must not depend on how the work was divided.
@test_throws diversityrefusal getabundance(eco)
@test length(getmetaabundance(eco)) == numSpecies
@test sum(getmetaabundance(eco)) ≈ 1.0
@test sum(getweight(eco)) ≈ 1.0
@test isapprox(sum(Diversity.API._getscale(eco)), 1.0)

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

# Global again - the second ecosystem's split is just as uneven as the first's (see the note above).
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
# inside an `mpiexec` child - several of them at once - must never happen. `canonical_reference()`
# only reads.
function checkblessed(abuns, prefix)
    reference = Canonical.canonical_reference()
    grid = reshape(abuns, numSpecies, VARYING_NY, VARYING_NX)
    for (key, value) in ("$prefix/total_abundance" => sum(grid),
        "$prefix/abundance_by_species" => vec(sum(grid, dims = (2, 3))),
        "$prefix/abundance_by_row" => vec(sum(grid, dims = (1, 3))),
        "$prefix/abundance_by_column" => vec(sum(grid, dims = (1, 2))))
        # A missing key means the canonical set has not been blessed here; say so rather than
        # passing silently, which would make this whole check vacuous.
        @test haskey(reference, key)
        haskey(reference, key) &&
            @test isapprox(float.(collect(value)), reference[key],
                           rtol = 1e-8)
    end
end

rank == 0 && checkblessed(true_abuns, "mpi")

# **`AlwaysMovement` disperses the standing population, and only a multi-rank run can check it.**
# It reads that population back out of the landscape through `EcoSISTEM._standingpopulation`, which
# maps the global species index onto this rank's local row - and at **one** rank that map is the
# identity, so a run there passes whether or not the mapping is right. Measured: replacing the
# mapping with the raw global index leaves the 1-rank answer bit-identical and makes 2 ranks die in
# `Multinomial` on a negative count. That is why this is pinned to a blessed serial number and
# asserted at every rank count, rather than compared against a serial run at one.
alwayssppl, _ = mpifixture_species(movement = mpifixture_always())
alwayseco = MPIEcosystem(alwayssppl, varying_environment(), nichefit, seed = 0)
alwayseco.abundances.rows_matrix .= MPIFIXTURE_FILL
MPI.Barrier(comm)
@test_nowarn simulate!(alwayseco, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP)
always_abuns = gatherabundance(alwayseco)

rank == 0 && checkblessed(always_abuns, "mpi/always")

# The same one-rank equality the birth-only run gets above. It cannot catch a rank-mapping bug (see
# the note above), but it does catch the loop body diverging from the serial one, which is A22's
# failure and is invisible to any MPI-against-MPI comparison.
if MPI.Comm_size(comm) == 1
    alwaysserial = mpifixture_ecosystem(movement = mpifixture_always())
    simulate!(alwaysserial, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP)
    @test alwaysserial.abundances.matrix == always_abuns
end

# **An intervention is the only thing that can leave the two landscape layouts out of step**, and so
# the only thing that exercises where `synchronise_from_rows!` sits in the timestep. It writes
# `rows_matrix` after the dynamics; the next timestep opens with `update_resource_usage!`, which
# reads the *column* layout. A sync placed before the interventions leaves that read one intervention
# behind, and the distributed run then diverges from the serial one.
#
# Measured before the sync was moved: with this fixture, serial totalled 1752450322 against MPI's
# 1752458039 -- at a **single rank**, so it was the ordering rather than the partition. No other test
# in the suite runs an intervention under MPI at all, which is why nothing caught it.
ivsppl, _ = mpifixture_species()
iveco = MPIEcosystem(ivsppl, varying_environment(), nichefit, seed = 0)
iveco.abundances.rows_matrix .= MPIFIXTURE_FILL
MPI.Barrier(comm)
@test_nowarn simulate!(iveco, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP,
                       intervention = mpifixture_intervention())
iv_abuns = gatherabundance(iveco)

rank == 0 && checkblessed(iv_abuns, "mpi/intervention")

# **A callback on `simulate!` under MPI.** Every rank runs the loop and the callback, so the
# occurrences must be the same on every rank, and a callback reaching a collective - a gather - must
# complete, since every rank reaches it. The totals it records must equal a serial run of the same
# fixture at any rank count: the distributed loop is a duplicate of the serial one, and a duplicate is
# only ever checked against the original.
cbsppl, _ = mpifixture_species()
cbeco = MPIEcosystem(cbsppl, varying_environment(), nichefit, seed = 0)
cbeco.abundances.rows_matrix .= MPIFIXTURE_FILL
cbcounts = Int[]
cbtotals = Int[]
MPI.Barrier(comm)
simulate!(cbeco, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP,
          every = EveryInterval(3 * MPIFIXTURE_TIMESTEP)) do occurrence
    push!(cbcounts, occurrence.count)
    return push!(cbtotals, sum(gatherabundance(cbeco)))
end
cbeverywhere = MPI.Allgather(Int32(length(cbcounts)), comm)
@test all(==(first(cbeverywhere)), cbeverywhere)
# Building a habitat is collective under MPI, so every rank builds each serial reference and only
# the root runs it.
cbserial = mpifixture_ecosystem()
if rank == 0
    serialtotals = Int[]
    simulate!(cbserial, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP,
              every = EveryInterval(3 * MPIFIXTURE_TIMESTEP)) do _
        return push!(serialtotals, sum(cbserial.abundances.matrix))
    end
    @test cbcounts == eachindex(serialtotals)
    @test cbtotals == serialtotals
end

# **A table of abundances under MPI**, against the serial twin at any rank count. Each rank adds only
# its own species' rows, and a table read row by row is read whole on every rank, so a table offering
# columns and the same rows streamed must both reproduce the serial run. Rows cover every species,
# several cells and fractional counts, so a rank adding another's rows or rounding before summing
# would show.
tablerows = [(species = sp, cell = cell, count = 0.25 * k + sp,
              time = k * MPIFIXTURE_TIMESTEP)
             for k in 0:23
             for sp in 1:numSpecies
             for cell in (1, 20 + k, VARYING_NY * VARYING_NX - sp)]
# Reversed, so the columns are out of time order and searched through their permutation.
tablecolumns = (species = reverse([r.species for r in tablerows]),
                cell = reverse([r.cell for r in tablerows]),
                count = reverse([r.count for r in tablerows]),
                time = reverse([r.time for r in tablerows]))
function tableiv(source)
    return Intervention(EveryStep(), ActiveCells(),
                        AddAbundanceTable(source))
end
function tablerun(source)
    sppl, _ = mpifixture_species()
    built = MPIEcosystem(sppl, varying_environment(), nichefit, seed = 0)
    built.abundances.rows_matrix .= MPIFIXTURE_FILL
    MPI.Barrier(comm)
    simulate!(built, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP,
              intervention = tableiv(source))
    return gatherabundance(built)
end
table_abuns = tablerun(tablecolumns)
stream_abuns = tablerun(r for r in tablerows)
tableserial = mpifixture_ecosystem()
if rank == 0
    simulate!(tableserial, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP,
              intervention = tableiv(tablecolumns))
    @test tableserial.abundances.matrix != true_abuns
    @test table_abuns == tableserial.abundances.matrix
    @test stream_abuns == tableserial.abundances.matrix
end

# **The recorders under MPI**, against the serial twin at any rank count. The abundance recorder
# gathers to the root, the diversity recorder to every rank, so both must reproduce the serial
# recording, and each keeps the same record of the run on every rank.
recsppl, _ = mpifixture_species()
receco = MPIEcosystem(recsppl, varying_environment(), nichefit, seed = 0)
receco.abundances.rows_matrix .= MPIFIXTURE_FILL
recevery = EveryInterval(3 * MPIFIXTURE_TIMESTEP)
nrec = length((0year):(3 * MPIFIXTURE_TIMESTEP):MPIFIXTURE_BURNIN)
ncells = VARYING_NY * VARYING_NX
recab = RecordAbundance(zeros(Int, numSpecies, ncells, nrec))
recdiv = RecordDiversity(zeros(ncells, 2, nrec), norm_sub_alpha, [0.0, 1.0])
MPI.Barrier(comm)
simulate!(receco, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP,
          every = recevery) do occurrence
    recab(receco, occurrence)
    return recdiv(receco, occurrence)
end
# Each holds the run as it stood at its last write: the last multiple of the interval, a step before
# the run ends.
@test provenance(recdiv).run == provenance(recab).run
@test provenance(recdiv).run.elapsed ≈ uconvert(u"s", MPIFIXTURE_BURNIN)
recserial = mpifixture_ecosystem()
if rank == 0
    serialab = RecordAbundance(zeros(Int, numSpecies, ncells, nrec))
    serialdiv = RecordDiversity(zeros(ncells, 2, nrec), norm_sub_alpha,
                                [0.0, 1.0])
    simulate!(recserial, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP,
              every = recevery) do occurrence
        serialab(recserial, occurrence)
        return serialdiv(recserial, occurrence)
    end
    @test recab.storage == serialab.storage
    @test recdiv.storage ≈ serialdiv.storage
end
# A metacommunity measure has no value per cell to assemble. Every rank refuses it alike before the
# gather, so none is left waiting in the collective and the run carries on past it.
#
# **The refusal met here is the measure's own**, not `gatherdiversity`'s level check: both compute
# the measure first, and a metacommunity value is refused for a distributed ecosystem however it is
# asked for (below). The level check still guards the serial recorder - `test_Recorder.jl` holds
# that one - and still fires here for an individual-level measure.
@test_throws "differ between ranks" gatherdiversity(receco, meta_gamma,
                                                    [0.0, 1.0])
@test_throws "differ between ranks" RecordDiversity(zeros(1, 2, nrec),
                                                    meta_gamma,
                                                    [0.0, 1.0])(receco,
                                                                (count = 1,
                                                                 elapsed = 0.0u"s",
                                                                 date = nothing))

# **And asking for one directly is refused too, at every rank count.** A measure built on a
# distributed ecosystem covers this rank's cells, so a metacommunity value taken from it answers for
# that rank alone: measured at two ranks, `meta_gamma` gave 6.2237 against the serial 6.2512, the
# ranks disagreed, and nothing was said. Every `meta_*` shorthand reaches `metadiv`, so one of them
# stands for all.
@test_throws "differ between ranks" meta_gamma(receco, 1.0)
@test_throws "differ between ranks" norm_meta_alpha(receco, [0.0, 1.0])
# The subcommunity form still answers, for this rank's cells, which is what `gatherdiversity`
# assembles.
@test nrow(norm_sub_alpha(receco, 1.0)) ==
      receco.abundances.cols_tuple.last - receco.abundances.cols_tuple.first + 1
MPI.Barrier(comm)

# **Ordinariness is computed in the COLUMN partition**, where a rank owns every species for its own
# cells, so a cell's value is complete on the rank that owns it and the full species-by-species
# similarity matrix applies with no slice. Gathering the column blocks must therefore rebuild the
# serial matrix exactly.
#
# In the row partition the similarity matrix had to be cut to this rank's species on both axes,
# which silently discarded every similarity between a species here and one elsewhere.
#
# This check alone cannot see the similarity half: the species list here is `UniqueTypes`, whose
# similarity matrix is the identity, so the discarded off-diagonal was all zeros and the old code
# got the right answer anyway. The `GeneralTypes` check below is what covers it - measured, a
# mutation dropping every off-diagonal similarity leaves this assertion passing and fails that one.
ordblock = getordinariness!(eco)
@test size(ordblock, 1) == numSpecies
@test size(ordblock, 2) ==
      eco.abundances.cols_tuple.last - eco.abundances.cols_tuple.first + 1

ordcounts = Int32.(size(ordblock, 1) .* eco.sccounts)
ordfull = similar(vec(ordblock), Int64(sum(ordcounts)))
MPI.Allgatherv!(vec(ordblock), MPI.VBuffer(ordfull, ordcounts), comm)
ordgathered = reshape(ordfull, numSpecies, VARYING_NY * VARYING_NX)
@test sum(ordgathered) ≈ 1.0

# Compared against serial at **every** rank count, not just at one. The serial run is built with
# `Ecosystem` directly; every rank builds it, since building a habitat is collective, and rank 0 runs
# it. At a single rank the comparison is weak - a rank owning every species cannot notice a
# calculation restricted to its own species - so the multi-rank runs are the ones that carry the
# check.
ordserial = mpifixture_ecosystem()
if rank == 0
    simulate!(ordserial, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP)
    @test ordgathered == getordinariness!(ordserial)
end

# **The same again with a similarity matrix that is NOT the identity**, which is the only fixture
# here that can fail on a calculation dropping similarities between species held on different ranks.
# Ordinariness multiplies the similarity matrix by the abundances, so a species' value in a cell
# depends on every other species in that cell; restricting the matrix to one rank's species discards
# every off-diagonal entry crossing the partition. Under `UniqueTypes` those entries are all zero,
# so the check above passes either way - see `mpifixture_similarity`.
# **The diversity measures themselves, assembled across ranks.** Each rank computes the measure for
# its own cells - complete, because the column partition puts every species of a cell on one rank -
# and `gatherdiversity` concatenates rather than combining. So the answer must equal the serial one
# at every rank count, which is the reproducibility requirement the whole distributed design rests
# on.
#
# `sub_gamma` is included alongside a normalised measure because it divides by the metacommunity
# ordinariness, a sum over *every* cell: without the `Allreduce` in `_getmetaordinariness!` it would
# silently use this rank's share of the metacommunity as though it were all of it.
divserial = mpifixture_ecosystem()
rank == 0 && simulate!(divserial, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP)
for divmeasure in (norm_sub_alpha, sub_gamma)
    for order in (1.0, [0.0, 1.0, 2.0])
        assembled = gatherdiversity(eco, divmeasure, order)
        if rank == 0
            wanted = divmeasure(divserial, order)
            @test assembled[!, :diversity] ≈ wanted[!, :diversity]
            @test assembled[!, :partition_name] == wanted[!, :partition_name]
            @test assembled[!, :q] == wanted[!, :q]
        end
    end
end

gensppl, _ = mpifixture_generalspecies()
geneco = MPIEcosystem(gensppl, varying_environment(), nichefit, seed = 0)
geneco.abundances.rows_matrix .= MPIFIXTURE_FILL
MPI.Barrier(comm)
@test_nowarn simulate!(geneco, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP)

genord = getordinariness!(geneco)
gencounts = Int32.(size(genord, 1) .* geneco.sccounts)
genfull = similar(vec(genord), Int64(sum(gencounts)))
MPI.Allgatherv!(vec(genord), MPI.VBuffer(genfull, gencounts), comm)
gengathered = reshape(genfull, numSpecies, VARYING_NY * VARYING_NX)

genserial = mpifixture_generalecosystem()
if rank == 0
    simulate!(genserial, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP)
    @test gengathered ≈ getordinariness!(genserial)
end

# The two layouts must hold the same data once a timestep has finished. Asserted globally rather than
# per rank: they are the same data under two decompositions, and only the totals are comparable
# without reindexing.
@test MPI.Allreduce(sum(iveco.abundances.rows_matrix), +, comm) ==
      MPI.Allreduce(sum(iveco.abundances.cols_vector), +, comm)

if MPI.Comm_size(comm) == 1
    ivserial = mpifixture_ecosystem()
    simulate!(ivserial, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP,
              intervention = mpifixture_intervention())
    @test ivserial.abundances.matrix == iv_abuns
end

# **The distributed hot loop's own allocation and inference checks.** `test/test_dynamics.jl` has
# these for the serial loop, and until now nothing had them for this one - which matters because the
# distributed `update!` and `update_resource_usage!` are a **separate implementation**, so a serial
# check cannot speak for them.
#
# Every rank must run the measurement, not just rank 0: `update!` is collective, so a rank that
# skipped it would leave the others waiting inside an `Alltoallv`.
function _mpiupdatealloc(eco, timestep)
    EcoSISTEM.update!(eco, timestep)
    EcoSISTEM.update!(eco, timestep)
    return @allocated EcoSISTEM.update!(eco, timestep)
end

function _mpiecosystem(numspecies)
    sppl, _ = mpifixture_species(numspecies = numspecies)
    built = MPIEcosystem(sppl, varying_environment(), nichefit, seed = 0)
    built.abundances.rows_matrix .= MPIFIXTURE_FILL
    return built
end

fewalloc = _mpiupdatealloc(_mpiecosystem(VARYING_SPECIES), MPIFIXTURE_TIMESTEP)
manyalloc = _mpiupdatealloc(_mpiecosystem(8 * VARYING_SPECIES),
                            MPIFIXTURE_TIMESTEP)

if rank == 0
    # A slope, as in the serial checks, because the threading and the collectives both cost a
    # constant per call that cancels from a difference.
    #
    # The bound is looser than the serial one's byte-per-species: an `Allgatherv` buffer is sized
    # from the partition, so at more than one rank the figure wobbles by a couple of hundred bytes
    # between species counts without any trend. Measured at 8 bytes per cell per species before the
    # habitat's topology was parameterised, which over these 49 extra species is about 30 kB - so
    # this bound still catches that by a factor of about seventy.
    extraspecies = 8 * VARYING_SPECIES - VARYING_SPECIES
    @test manyalloc - fewalloc < 8 * extraspecies

    # The same sweep `test_dynamics.jl` runs, over the distributed types. `epoch`, `calendar` and
    # `active` are abstract by design and are explained there; anything else must be looked at.
    mpiabstract = Tuple{Symbol, Symbol}[]
    for obj in (eco, eco.abundances, eco.habitat, eco.cache)
        S = typeof(obj)
        for f in fieldnames(S)
            isconcretetype(fieldtype(S, f)) && continue
            (f === :epoch || f === :calendar || f === :active) && continue
            push!(mpiabstract, (nameof(S), f))
        end
    end
    @test mpiabstract == Tuple{Symbol, Symbol}[]

    # And the paths the distributed loop actually reads, which a field walk cannot see: these go
    # *through* the landscape, whose own container is abstract by design.
    E = typeof(eco)
    @test isconcretetype(Base.return_types(e -> e.abundances.rows_matrix, (E,))[1])
    @test isconcretetype(Base.return_types(e -> e.cache.totaldemand, (E,))[1])
    @test isconcretetype(Base.return_types(e -> e.cache.netmigration, (E,))[1])
    @test isconcretetype(Base.return_types(e -> e.habitat.topology, (E,))[1])
end

# **Every rank takes the first rank's seed and species list.** Anything drawn at random on each rank
# - a seed nobody gave, abundances split without a seed, a random phylogeny - would otherwise give
# each rank a different ecosystem, and a random intervention would select different cells on each.
@testset "every rank takes the first rank's seed and species list" begin
    everyone(x) = MPI.Allgather(x, comm)
    agree(x) = all(==(first(x)), x)
    gathered(x) = MPI.gather(x, comm)
    agreeing(x) = (all = gathered(x); rank == 0 ? agree(all) : true)

    # No seed: the first rank draws one, and a random selection is the same everywhere.
    sppl, _ = mpifixture_species()
    unseeded = MPIEcosystem(sppl, varying_environment(), nichefit)
    @test agree(everyone(unseeded.seed))
    unseeded.abundances.rows_matrix .= MPIFIXTURE_FILL
    MPI.Barrier(comm)
    simulate!(unseeded, 3 * MPIFIXTURE_TIMESTEP, MPIFIXTURE_TIMESTEP,
              intervention = Intervention(AtTime(1.0month_mean_duration),
                                          RandomCells(20), Deactivate()))
    @test agreeing(unseeded.habitat.active)
    @test count(.!unseeded.habitat.active) == 20

    # Abundances drawn differently on each rank: every rank simulates the first rank's.
    mine, _ = mpifixture_species()
    mine.abun .= rand(Xoshiro(rank), 1:100, length(mine.abun))
    shared = MPIEcosystem(mine, varying_environment(), nichefit, seed = 0)
    @test agreeing(shared.spplist.abun)
    @test agree(everyone(sum(shared.spplist.abun)))
    @test MPI.Allreduce(sum(shared.abundances.rows_matrix), +, comm) ==
          sum(shared.spplist.abun)

    # A random phylogeny, and the tolerances evolved on it, travel whole.
    kernels = GaussianKernel.(fill(1.0km, numSpecies), 1e-4)
    treed = SpeciesList(numSpecies, 2, fill(10, numSpecies),
                        Demand{SolarRadiation}(fill(1.0kJ / day, numSpecies)),
                        BirthOnlyMovement(kernels),
                        EqualPop(0.6 / year, 0.6 / year, 1.0, 0.2),
                        fill(true, numSpecies))
    ext = Base.get_extension(EcoSISTEM, :EcoSISTEMMPIExt)
    received = ext._rootspecies(treed, comm)
    @test typeof(received) === typeof(treed)
    tree = received.types.tree
    @test agreeing(sort([getlength(tree, b) for b in getbranches(tree)]))
    @test agreeing(received.tolerance.vals)

    # A seed must be the same on every rank, or given on none.
    if MPI.Comm_size(comm) > 1
        @test_throws "the same `seed` on every rank" MPIEcosystem(mpifixture_species()[1],
                                                                  varying_environment(),
                                                                  nichefit,
                                                                  seed = rank)
        @test_throws "rank 1: none" MPIEcosystem(mpifixture_species()[1],
                                                 varying_environment(),
                                                 nichefit,
                                                 seed = rank == 0 ? 1 :
                                                        nothing)
    end
    @test MPIEcosystem(mpifixture_species()[1], varying_environment(), nichefit,
                       seed = 7).seed == 7
    MPI.Barrier(comm)
end

if !MPI.Finalized()
    MPI.Finalize()
end
