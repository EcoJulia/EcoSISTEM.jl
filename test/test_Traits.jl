using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units

import Simulation: DiscreteTrait

opts = fill(5.0°C, numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies)  * °C
@test_nowarn traits = GaussTrait(opts, vars)
@test_nowarn DiscreteTrait(fill(1, 10))
@test_nowarn TempBin(fill(1, 10, 2))
@test_nowarn RainBin(fill(1, 10, 2))
@test_nowarn TraitCollection2(TempBin(fill(1, 10, 2)), RainBin(fill(1, 10, 2)))
