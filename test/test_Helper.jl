module TestHelper

using EcoSISTEM
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Diversity

include("TestCases.jl")

function TempIncrease!(eco::Ecosystem, timestep::Unitful.Time, rate::typeof(1.0K/year))
    resetrate!(eco, rate)
    eco.abenv.habitat.matrix[eco.abenv.habitat.matrix .< 0K] .= 0K
end
@testset "Simulate functions" begin
    eco = Test1Ecosystem()

    times = 3.0months; burnin = 1.0months; interval = 1.0month
    timestep = 1month
    # Run simulation grid
    lensim = length(0.0month:interval:times)
    @testset "simulate" begin
        # Run simulations 10 times
        @test_nowarn generate_storage(eco, lensim, 1)
        abun = generate_storage(eco, lensim, 1)
        @test size(abun) == (size(eco.abundances.matrix, 1), size(eco.abundances.matrix, 2), lensim, 1)
        @test_nowarn simulate!(eco, burnin, timestep)
        @test_nowarn simulate_record!(abun, eco, times, interval, timestep)
        isdir("data") || mkdir("data")
        @test_nowarn simulate!(eco, times, interval, timestep, "data", "testrun")
        rm("data", recursive = true)
    end
    @testset "scenarios" begin
        # Run with scenarios
        eco = TestMultiEcosystem()
        abun = generate_storage(eco, lensim, 1)
        scenario = SimpleScenario(TempIncrease!, 0.0K/year)
        @test_nowarn simulate_record!(abun, eco, times, interval, timestep, scenario)
    end
    @testset "diversity simulate" begin
        # Run diversity simulations 10 times
        qs = collect(1.0:3)
        @test_nowarn generate_storage(eco, length(qs), lensim, 1)
        abun = generate_storage(eco, length(qs), lensim, 1)
        @test size(abun) == (size(eco.abundances.matrix, 2), length(qs), lensim, 1)
        @test_nowarn simulate!(eco, burnin, timestep)
        @test_nowarn simulate_record_diversity!(abun, eco, times, interval, timestep, norm_sub_alpha, qs)

        abun = generate_storage(eco, length(qs), lensim, 1)
        scenario = SimpleScenario(TempIncrease!, 0.0K/year)
        @test_nowarn simulate_record_diversity!(abun, eco, times, interval, timestep, scenario, norm_sub_alpha, qs)

        divfuns = [norm_sub_alpha, norm_sub_beta]
        abun = generate_storage(eco, length(divfuns), lensim, 1)
        @test_nowarn simulate_record_diversity!(abun, eco, times, interval, timestep, divfuns, 1.0)

        divfuns = [norm_sub_alpha, norm_sub_beta]
        scenario = SimpleScenario(TempIncrease!, 0.0K/year)
        abun = generate_storage(eco, length(divfuns), lensim, 1)
        @test_nowarn simulate_record_diversity!(abun, eco, times, interval, timestep, scenario, divfuns, 1.0)

        qs = collect(1.0:3)
        abun1 = zeros(Float64, 100, 3, 3, 4)
        abun2 = zeros(Float64, 3, 3, 4)
        @test_nowarn simulate_record_diversity!(abun1, abun2, eco, times, interval, timestep, qs)
    end
end

end
