using Simulation
using Unitful.DefaultSymbols
using Test
using Simulation.Units

birth = 0.6/month
death = 0.6/month
l = 1.0
s = 0.0
boost = 1000.0
timestep = 1.0month

@test_nowarn param = EqualPop(birth, death, l, s, boost)
@test_nowarn param = PopGrowth{typeof(unit(0.0/month))}(fill(birth, 5), fill(death, 5), l, s, boost)
@test_nowarn param = NoGrowth{typeof(unit(0.0/month))}(fill(birth, 5), fill(death, 5), l, s, boost)
