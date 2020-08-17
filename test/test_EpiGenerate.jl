using Simulation
using Test
using Unitful.DefaultSymbols
using Distributions
using Simulation.Units

import Simulation: virusupdate!, classupdate!, invalidatecaches!
include("TestCases.jl")

epi = TestEpiSystem()
totalpops = sum(epi.abundances.matrix, dims = 2)
@test_nowarn update!(epi, 1day)
update!(epi, 1day)
@test all(epi.abundances.matrix .>= 0)
@test_nowarn virusupdate!(epi, 1day)
@test_nowarn classupdate!(epi, 1day)
invalidatecaches!(epi)
@test sum(epi.cache.virusdecay) == 0
@test sum(epi.cache.virusmigration) == 0

epi = TestEpiLockdown()
@test_nowarn update!(epi, 1day)
update!(epi, 1day)
