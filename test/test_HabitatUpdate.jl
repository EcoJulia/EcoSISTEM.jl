using Simulation
using Distributions
using Test
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

include("TestCases.jl")

eco = TestEcosystem()
@test_nowarn Simulation.HabitatUpdate(
    Simulation.NoChange, 0.0/month, Unitful.Dimensions{()})
@test_nowarn Simulation.HabitatUpdate(
    Simulation.TempChange, 1.0K/month, typeof(dimension(1.0K)))
@test_nowarn Simulation.habitatupdate!(eco, 1month)
@test_nowarn Simulation.budgetupdate!(eco, 1month)
