using Simulation
using Base.Test
using Distributions

opts = repmat([5.0], numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies)
@test_nowarn traits = TempTrait(opts, vars)
