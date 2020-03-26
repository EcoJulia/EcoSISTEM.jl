using Simulation
using Test
using Distributions
using RCall
using Unitful.DefaultSymbols
using Simulation.Units

include("TestCases.jl")

eco = TestEcosystem()

times = 10.0months; burnin = 3.0months; interval = 1.0month
timestep = 1month
# Run simulation grid
lensim = length(0.0month:interval:times)

# Run simulations 10 times
@test_nowarn generate_storage(eco, lensim, 1)
abun = generate_storage(eco, lensim, 1)
@test_nowarn simulate!(eco, burnin, timestep)
@test_nowarn simulate_record!(abun, eco, times, interval, timestep)
