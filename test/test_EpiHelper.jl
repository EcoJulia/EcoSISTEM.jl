using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units

include("TestCases.jl")

epi = TestEpiSystem()

@testset "EpiHelper" begin
    times = 1.0month; burnin = 1.0month; interval = 1.0day
    timestep = 1.0day
    # Run simulation grid
    lensim = length(0.0month:interval:times)

    # Run simulations 10 times
    abun =  zeros(Int64, size(epi.abundances.matrix, 1), size(epi.abundances.matrix, 2), lensim)
    @test_nowarn simulate!(epi, burnin, timestep; save=true, save_path="TEMP")
    # Check if everything was saved
    @test isdir("TEMP")
    @test isfile(joinpath("TEMP", "configuration.jlso"))
    @test isfile(joinpath("TEMP", "initial_system.jlso"))
    @test isfile(joinpath("TEMP", "final_system.jlso"))
    # Delete temporary directory
    rm("TEMP", recursive=true)

    @test_nowarn simulate_record!(abun, epi, times, interval, timestep; save=true, save_path="TEMP")
    # Check if everything was saved
    @test isdir("TEMP")
    @test isfile(joinpath("TEMP", "configuration.jlso"))
    @test isfile(joinpath("TEMP", "initial_system.jlso"))
    @test isfile(joinpath("TEMP", "final_system.jlso"))
    @test isfile(joinpath("TEMP", "abundances.jlso"))
    # Delete temporary directory
    rm("TEMP", recursive=true)
    @test all(abun .>= 0)
    @test all(sum(abun, dims = (1,2)) .> 0)
end
