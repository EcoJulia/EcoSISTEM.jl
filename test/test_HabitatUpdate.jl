using EcoSISTEM
using Distributions
using Compat.Test
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units

include("TestCases.jl")

eco = TestEcosystem()
@test_nowarn EcoSISTEM.HabitatUpdate{Unitful.Dimensions{()}}(EcoSISTEM.NoChange,
    0.0/month)
@test_nowarn EcoSISTEM.HabitatUpdate{typeof(dimension(1.0K))}(EcoSISTEM.TempChange,
    1.0K/month)
@test_nowarn EcoSISTEM.habitatupdate!(eco, 1month)
@test_nowarn EcoSISTEM.budgetupdate!(eco, 1month)
