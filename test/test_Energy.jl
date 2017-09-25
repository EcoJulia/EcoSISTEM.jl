using Simulation
using Base.Test

@test_nowarn energy_vec = SimpleRequirement(repmat([2], 10))
@test_nowarn energy_vec = SimpleRequirement(2:10)
@test_nowarn energy_vec = SimpleRequirement(2.0:10.0)
bud = Array{Float64, 2}(2, 2)
fill!(bud, 100.0)
@test_nowarn Simulation.SimpleBudget(bud)
