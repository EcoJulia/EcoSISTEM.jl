# SPDX-License-Identifier: LGPL-3.0-or-later

module ExtEcoSISTEMMPIExt

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using Diversity
using MPI
using JLD2
using Test

if !MPI.Initialized()
    MPI.Init()
end

@testset "MPI" begin
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println(Threads.nthreads())
    numSpecies = 100
    grid = (10, 10)
    demand = 10.0kJ / day
    individuals = 1_000
    area = 100.0 * km^2
    totalK = 100.0kJ / km^2 / day
    # Set up initial parameters for ecosystem

    # Set up how much resource each species consumes
    resource_vec = Demand{SolarRadiation}(fill(demand, numSpecies))

    # Set probabilities
    birth = 0.6 / year
    death = 0.6 / year
    longevity = 1.0
    survival = 0.2

    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival)

    # Create kernel for movement
    kernel = fill(GaussianKernel(1.0km, 10e-10), numSpecies)
    movement = BirthOnlyMovement(kernel)

    # Create species list, including their temperature preferences, seed abundance and native status
    opts = fill(274.0K, numSpecies)
    vars = fill(0.5K, numSpecies)
    tolerance = NicheTolerance(Temperature, Normal, opts, vars)
    native = fill(true, numSpecies)
    # abun = rand(Multinomial(individuals, numSpecies))
    abun = fill(div(individuals, numSpecies), numSpecies)
    sppl = SpeciesList(numSpecies, tolerance, abun, resource_vec,
                       movement, param, native)

    # Create abiotic environment - an even grid of one temperature. The study area decides the
    # grid; `GridHabitat` only samples onto it.
    studyarea = StudyArea(extent = (sqrt(area), sqrt(area)),
                          cellsize = sqrt(area) / grid[1], verbosity = :silent)
    habitat = GridHabitat(regime = UniformSpec(274.0K,
                                               axis = Temperature),
                          supply = UniformSpec(totalK,
                                               axis = SolarRadiation),
                          area = studyarea)

    # Set nichefit between species and environment (gaussian)
    nichefit = NicheSuitability{Temperature, typeof(1.0K)}()

    # build_ecosystem auto-selection: with MPI initialised but a single rank, `:auto` stays serial
    # (Comm_size == 1) while `distributed = true` forces the distributed type - from the same
    # SpeciesList + GridHabitat that MPIEcosystem takes directly.
    @test build_ecosystem(sppl, habitat, nichefit = nichefit) isa Ecosystem
    @test build_ecosystem(sppl, habitat, nichefit = nichefit,
                          distributed = true) isa
          MPIEcosystem

    # Create ecosystem
    @test_nowarn eco = MPIEcosystem(sppl, habitat, nichefit)
    eco = MPIEcosystem(sppl, habitat, nichefit)
    @test sum(eco.sppcounts) == length(eco.spplist.names)
    @test eco.firstsp == 1
    @test sum(eco.sccounts) == prod(size(eco.habitat.regime.matrix))
    @test eco.firstsc == 1

    # Simulation Parameters
    burnin = 1month_mean_duration
    times = 3month_mean_duration
    timestep = 1month_mean_duration
    record_interval = 1month_mean_duration
    repeats = 1
    lensim = length((0year):record_interval:times)
    # Burnin
    MPI.Barrier(comm)
    @time simulate!(eco, burnin, timestep)

    # Set columns vector to zero and check synchronise from rows
    eco.abundances.cols_vector .= 0
    EcoSISTEM.synchronise_from_rows!(eco.abundances)
    @test sum(eco.abundances.cols_vector) == sum(eco.abundances.rows_matrix)

    # Set rows matrix to zero and check synchronise from cols
    eco.abundances.rows_matrix .= 0
    EcoSISTEM.synchronise_from_cols!(eco.abundances)
    @test sum(eco.abundances.cols_vector) == sum(eco.abundances.rows_matrix)

    # `getabundance` refuses on a distributed ecosystem whatever the rank count: it asks for the
    # whole metacommunity's abundances, which no rank holds, and Diversity's consumers each reduce
    # that matrix over a different axis. The refusal stands even here, where a single rank does
    # happen to hold everything - an answer that appeared at one rank and vanished at two would be
    # exactly the kind of rank-dependence the distributed code exists to avoid.
    @test_throws "no single rank holds them" getabundance(eco)

    # The genuinely global quantities, asserted on their length as well as their total: a sum of one
    # is true of any block normalised by its own total, which is how the per-species totals came to
    # cover only the calling rank's species without any test noticing.
    @test length(getmetaabundance(eco)) == numSpecies
    @test sum(getmetaabundance(eco)) ≈ 1.0
    @test sum(getordinariness!(eco)) ≈ 1.0
    @test sum(getweight(eco)) ≈ 1.0
    @test length(getweight(eco)) == countsubcommunities(eco)

    # Gather abundances and check against rows matrix - should be same for 1 process
    abuns = gatherabundance(eco)
    @test abuns == eco.abundances.rows_matrix
end

# Run one `mpiexec` launch, and fail rather than hang if its ranks never finish: a rank left waiting
# in a collective waits forever. A launch takes a few minutes, far inside the limit. Killing
# `mpiexec` stops its ranks too.
function runmpi(cmd; limit = 1800)
    process = run(cmd, wait = false)
    timedwait(() -> process_exited(process), limit, pollint = 1.0) === :ok &&
        return success(process)
    kill(process)
    wait(process)
    @error "`mpiexec` was stopped after $limit s without finishing: $cmd"
    return false
end

# The Julia command a launch runs, with coverage instrumentation only where it is wanted. It roughly
# doubles a launch - measured on a runner, this file took 861 s on the instrumented job against 402 s
# uninstrumented - and each launch runs the same package code, so one carrying it covers the same
# lines as three would. The parent process, which runs this file, keeps whatever it was given.
function launchcmd(; coverage::Bool)
    cmd = Base.julia_cmd()
    coverage && return cmd
    return Cmd(filter(a -> !startswith(a, "--code-coverage"), cmd.exec))
end

@testset "mpirun" begin
    # Keep the MPI outputs in a temp dir the OS cleans up (no manual `rm` needed for hygiene).
    # The child `mpiexec` processes read its path as their first command-line argument.
    datadir = mktempdir()
    # Compare 1 thread 4 processes vs. 4 threads 1 process vs. 2 threads 2 processes
    withenv("JULIA_NUM_THREADS" => "4") do
        nprocs = 1
        function cmd(n = nprocs)
            return `$(mpiexec()) -n $nprocs $(launchcmd(coverage = false)) --startup-file=no $(pkgdir(EcoSISTEM, "test", "SmallMPItest.jl")) $datadir`
        end
        @test runmpi(cmd())
    end
    withenv("JULIA_NUM_THREADS" => "2") do
        nprocs = 2
        function cmd(n = nprocs)
            return `$(mpiexec()) -n $nprocs $(launchcmd(coverage = true)) --startup-file=no $(pkgdir(EcoSISTEM, "test", "SmallMPItest.jl")) $datadir`
        end
        @test runmpi(cmd())
    end
    withenv("JULIA_NUM_THREADS" => "1") do
        nprocs = 4
        function cmd(n = nprocs)
            return `$(mpiexec()) -n $nprocs $(launchcmd(coverage = false)) --startup-file=no $(pkgdir(EcoSISTEM, "test", "SmallMPItest.jl")) $datadir`
        end
        @test runmpi(cmd())
    end

    ## All answers should be the same across process/thread splits. Load the
    ## saved values (note: `@load file var` binds `var` but returns the symbol
    ## list, so use `load(file, "abuns")` to get the data itself).
    abuns1thread = load(joinpath(datadir, "Test_abuns1.jld2"), "abuns")
    abuns2thread = load(joinpath(datadir, "Test_abuns2.jld2"), "abuns")
    abuns4thread = load(joinpath(datadir, "Test_abuns4.jld2"), "abuns")

    # Each launch also checks its own run against the blessed values, so a launch that does not
    # reproduce itself fails there, apart from any disagreement between the three here.
    @test abuns1thread == abuns2thread == abuns4thread
end

if !MPI.Finalized()
    MPI.Finalize()
end

end
