using Simulation
using Unitful.DefaultSymbols
using Test
using Simulation.Units
using Unitful
import Simulation.equalpop

birth = 0.6/month
death = 0.6/month
longevity = 1.0
survival = 0.0
boost = 1000.0
numSpecies = 10

param = EqualPop(birth, death, longevity, survival, boost)
@test_nowarn EqualPop(birth, death, longevity, survival, boost)
equalparams = equalpop(param, numSpecies)
@test length(equalparams.birth) == numSpecies
@test all(equalparams.birth .== birth)
@test all(equalparams.death .== death)
@test_nowarn param = PopGrowth{typeof(unit(0.0/month))}(fill(birth, 5), fill(death, 5), longevity, survival, boost)
@test_nowarn param = NoGrowth{typeof(unit(0.0/month))}(fill(birth, 5), fill(death, 5), longevity, survival, boost)

param = PopGrowth{typeof(unit(0.0/month))}(fill(birth, numSpecies), fill(death, numSpecies), longevity, survival, boost)
equalparams = equalpop(param, numSpecies)
@test length(equalparams.birth) == numSpecies
@test all(equalparams.birth .== birth)
@test all(equalparams.death .== death)

param = NoGrowth{typeof(unit(0.0/month))}(fill(birth, numSpecies), fill(death, numSpecies), longevity, survival, boost)
equalparams = equalpop(param, numSpecies)
@test length(equalparams.birth) == numSpecies
@test all(equalparams.birth .== birth)
@test all(equalparams.death .== death)
