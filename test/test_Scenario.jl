using EcoSISTEM
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Diversity


include("TestCases.jl")

@testset "Scenarios" begin
    # Simple scenario
    eco = TestMultiEcosystem()
    function TempIncrease!(eco::Ecosystem, timestep::Unitful.Time, rate::typeof(1.0K/year))
        resetrate!(eco, rate)
        eco.abenv.habitat.matrix[eco.abenv.habitat.matrix .< 0K] .= 0K
    end
    scenario = SimpleScenario(TempIncrease!, 1.0K/year)
    @test_nowarn EcoSISTEM.runscenario!(eco, 1month, scenario, 1month)
    EcoSISTEM.update!(eco, 1month)
    @test mean(eco.abenv.habitat.matrix) > 10.0K

    # Fluctuating scenario
    eco = TestMultiEcosystem()
    function TempFluct!(eco::Ecosystem, timestep::Unitful.Time, rate::Quantity{Float64, 𝚯*𝐓^-1}, currentstep::Unitful.Time, startarray::Matrix{typeof(1.0K)})
        v = uconvert(K/year, rate)
        eco.abenv.habitat.matrix .= (v * year) .* sin(collect(-π:π/6:π)[mod(Int64(uconvert(NoUnits, currentstep/timestep)),12) + 1]) .+ startarray
    end
    scenario2 = FluctScenario(TempFluct!, 1.0K/year, eco.abenv.habitat.matrix)
    @test_nowarn EcoSISTEM.runscenario!(eco, 1month, scenario2, 1month)
    @test mean(eco.abenv.habitat.matrix) < 10.0K

    # Mutliple scenarios
    multiscenario = MultiScenario(scenario, scenario2)
    @test_nowarn EcoSISTEM.runscenario!(eco, 1month, multiscenario, 1month)

end
