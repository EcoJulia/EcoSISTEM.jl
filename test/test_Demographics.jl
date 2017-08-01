using Simulation
using Base.Test

birth = 0.6
death = 0.6
l = 1.0
s = 0.0
boost = 1000.0
timestep = 1.0

@test_nowarn param = EqualPop(birth, death, l, s, boost)
