module TestMPI

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
    numSpecies = 100; grid = (10, 10); req= 10.0kJ; individuals=1_000; area = 100.0*km^2; totalK = 100.0kJ/km^2
    # Set up initial parameters for ecosystem

    # Set up how much energy each species consumes
    energy_vec = SolarRequirement(fill(req, numSpecies))

    # Set probabilities
    birth = 0.6/year
    death = 0.6/year
    longevity = 1.0
    survival = 0.2
    boost = 1.0

    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(1.0km, 10e-10), numSpecies)
    movement = BirthOnlyMovement(kernel, Torus())

    # Create species list, including their temperature preferences, seed abundance and native status
    opts = fill(274.0K, numSpecies)
    vars = fill(0.5K, numSpecies)
    traits = GaussTrait(opts, vars)
    native = fill(true, numSpecies)
    # abun = rand(Multinomial(individuals, numSpecies))
    abun = fill(div(individuals, numSpecies), numSpecies)
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)

    # Create abiotic environment - even grid of one temperature
    abenv = simplehabitatAE(274.0K, grid, totalK, area)

    # Set relationship between species and environment (gaussian)
    rel = Gauss{typeof(1.0K)}()

    # Create ecosystem
    @test_nowarn eco = MPIEcosystem(sppl, abenv, rel)
    eco = MPIEcosystem(sppl, abenv, rel)
    @test sum(eco.sppcounts) == length(eco.spplist.species.names)
    @test eco.firstsp == 1
    @test sum(eco.sccounts) == prod(size(eco.abenv.habitat.matrix))
    @test eco.firstsc == 1

    # Simulation Parameters
    burnin = 1month; times = 3months; timestep = 1month; record_interval = 1months; repeats = 1
    lensim = length(0years:record_interval:times)
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

    @test sum(getabundance(eco)) ≈ 1.0
    @test sum(getmetaabundance(eco)) ≈ 1.0
    @test sum(getordinariness!(eco)) ≈ 1.0
    @test sum(getweight(eco)) ≈ 1.0

    # Gather abundances and check against rows matrix - should be same for 1 process
    abuns = gather_abundance(eco)
    @test abuns == eco.abundances.rows_matrix
end

@testset "mpirun" begin
    # Keep outputs all one folder 
    isdir("data") || mkdir("data")
    testdir = @__DIR__
    # Compare 1 thread 4 processes vs. 4 threads 1 process vs. 2 threads 2 processes
    mpiexec() do mpirun
        withenv("JULIA_NUM_THREADS" => "4") do
            nprocs = 1
            cmd(n=nprocs) = `$mpirun -n $nprocs $(Base.julia_cmd()) --startup-file=no $(joinpath(testdir, "SmallMPItest.jl"))`
            @test success(run(cmd()))
        end
        withenv("JULIA_NUM_THREADS" => "2") do
            nprocs = 2
            cmd(n=nprocs) = `$mpirun -n $nprocs $(Base.julia_cmd()) --startup-file=no $(joinpath(testdir, "SmallMPItest.jl"))`
            @test success(run(cmd()))
        end
        withenv("JULIA_NUM_THREADS" => "1") do
            nprocs = 4
            cmd(n=nprocs) = `$mpirun -n $nprocs $(Base.julia_cmd()) --startup-file=no $(joinpath(testdir, "SmallMPItest.jl"))`
            @test success(run(cmd()))
        end
    end

    ## All answers should be the same
    abuns1thread = @load "data/Test_abuns1.jld2" abuns
    abuns2thread = @load "data/Test_abuns2.jld2" abuns
    abuns4thread = @load "data/Test_abuns4.jld2" abuns

    @test abuns1thread == abuns2thread == abuns4thread
end

# Clean up outputs
rm("data", recursive = true)

if !MPI.Finalized()
    MPI.Finalize()
end

end
