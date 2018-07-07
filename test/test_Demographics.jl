using Simulation
using Unitful.DefaultSymbols
using Base.Test
using myunitful

birth = 0.6/month
death = 0.6/month
l = 1.0
s = 0.0
boost = 1000.0
timestep = 1.0month

@test_nowarn param = EqualPop(birth, death, l, s, boost)
