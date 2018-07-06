using Simulation
using Base.Test
using Distributions
using Unitful.DefaultSymbols

opts = repmat([5.0°C], numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies)  * °C
@test_nowarn traits = GaussTrait(opts, vars)
