using Simulation
using Unitful.DefaultSymbols
using Compat.Test
using Simulation.Units
import Simulation.equalpop

birth = 0.6/month
death = 0.6/month
l = 1.0
s = 0.0
boost = 1000.0
numSpecies = 10

param = EqualPop(birth, death, l, s, boost)
@test_nowarn EqualPop(birth, death, l, s, boost)
equalparams = equalpop(param, numSpecies)
@test length(equalparams.birth) == numSpecies
@test all(equalparams.birth .== birth)
@test all(equalparams.death .== death)
@test_nowarn param = PopGrowth{typeof(unit(0.0/month))}(fill(birth, 5), fill(death, 5), l, s, boost)
@test_nowarn param = NoGrowth{typeof(unit(0.0/month))}(fill(birth, 5), fill(death, 5), l, s, boost)
