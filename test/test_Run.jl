using EcoSISTEM
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Diversity
using HDF5

include("TestCases.jl")

@testset "Simulate functions" begin
    @testset "Eco" begin
        eco = TestMultiEcosystem()

        times = 3.0months
        burnin = 1.0months
        interval = 1.0month
        timestep = 1month
        # Run simulation grid
        lensim = length((0.0month):interval:times)
        # Run simulations 10 times
        @test_nowarn generate_storage(eco, lensim, 1)
        abun = generate_storage(eco, lensim, 1)
        @test size(abun) ==
              (size(eco.abundances.matrix, 1), size(eco.abundances.matrix, 2),
               lensim, 1)
        @test_nowarn simulate!(eco, burnin, timestep)
        @test_nowarn simulate_record!(abun, eco, times, interval, timestep)
    end

    @testset "Epi" begin
        epi = Test1EpiSystem()
        times = 1.0month
        burnin = 1.0month
        interval = 1.0day
        timestep = 1.0day
        # Run simulation grid
        lensim = length((0.0month):interval:times)

        # Run simulations 10 times
        abun = zeros(Int64, size(epi.abundances.matrix, 1),
                     size(epi.abundances.matrix, 2), lensim)
        @test_nowarn simulate!(epi, burnin, timestep; save = true,
                               save_path = "TEMP")
        # Check if everything was saved
        @test isdir("TEMP")
        @test isfile(joinpath("TEMP", "configuration.jlso"))
        @test isfile(joinpath("TEMP", "initial_system.jlso"))
        @test isfile(joinpath("TEMP", "final_system.jlso"))
        # Delete temporary directory
        rm("TEMP", recursive = true)

        @test_nowarn simulate_record!(abun, epi, times, interval, timestep;
                                      save = true, save_path = "TEMP")
        # Check if everything was saved
        @test isdir("TEMP")
        @test isfile(joinpath("TEMP", "configuration.jlso"))
        @test isfile(joinpath("TEMP", "initial_system.jlso"))
        @test isfile(joinpath("TEMP", "final_system.jlso"))
        @test isfile(joinpath("TEMP", "abundances.h5"))
        # Check output abundance HDF5 file
        saved_abuns = h5read(joinpath("TEMP", "abundances.h5"),
                             "abundances/abuns")
        @test saved_abuns == abun
        # Delete temporary directory
        rm("TEMP", recursive = true)
        @test all(abun .>= 0)
        @test all(sum(abun, dims = (1, 2)) .> 0)
    end

    @testset "Transitions" begin
        eco = TestTransitions()
        timestep = 1month
        for su in eco.transitions.setup
            @test_nowarn run_rule!(eco, su, timestep)
        end
        for st in eco.transitions.state
            @test_nowarn run_rule!(eco, st, timestep)
        end
        for pl in eco.transitions.place
            @test_nowarn run_rule!(eco, pl, timestep)
        end
        for wd in eco.transitions.winddown
            @test_nowarn run_rule!(eco, wd, timestep)
        end
    end
end
