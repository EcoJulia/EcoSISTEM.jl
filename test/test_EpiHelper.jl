using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units

include("TestCases.jl")

epi = TestEpiSystem()

times = 1.0month; burnin = 1.0month; interval = 1.0day
timestep = 1.0day
# Run simulation grid
lensim = length(0.0month:interval:times)

# Run simulations 10 times
abun =  zeros(Int64, 4, size(epi.abundances.matrix, 2), lensim)
@test_nowarn simulate!(epi, burnin, timestep)
@test_nowarn simulate_record!(abun, epi, times, interval, timestep)
@test all(abun .>= 0)
@test all(sum(abun, dims = (1,2)) .> 0)
