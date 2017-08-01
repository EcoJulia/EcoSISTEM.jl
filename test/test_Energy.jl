using Simulation
using Base.Test

@test_nowarn energy_vec = SimpleRequirement(repmat([2], 10))
@test_nowarn energy_vec = SimpleRequirement(2:10)
@test_nowarn energy_vec = SimpleRequirement(2.0:10.0)
