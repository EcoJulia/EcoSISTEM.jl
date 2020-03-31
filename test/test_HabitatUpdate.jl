using Simulation
using Distributions
using Compat.Test
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

include("TestCases.jl")

eco = TestEcosystem()
@test_nowarn Simulation.HabitatUpdate{Unitful.Dimensions{()}}(Simulation.NoChange,
    0.0/month)
@test_nowarn Simulation.HabitatUpdate{typeof(dimension(1.0K))}(Simulation.TempChange,
    1.0K/month)
@test_nowarn Simulation.habitatupdate!(eco, 1month)
@test_nowarn Simulation.budgetupdate!(eco, 1month)
