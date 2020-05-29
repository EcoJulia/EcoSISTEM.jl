using Simulation
using Test
using Unitful.DefaultSymbols
using Simulation.Units

# Allocating version:
param = (beta_env = 1.0/day, beta_force = 1.0/day, sigma = 0.02/day, virus_growth = 1e-3/day, virus_decay = 1e-3/day, mean_dispersal_dist = 10.0km)
runparams = (times = 2years, interval = 1day, timestep = 1day)
grid_size = (4,4)
area = 100.0km^2
@test_nowarn abuns = SIR_wrapper(grid_size, area, param, runparams)

# Non-allocating version:
times = runparams.times; interval = runparams.interval
Ncells = grid_size[1] * grid_size[2]
numclasses = 4
abuns = zeros(Int64, numclasses, Ncells, convert(Int64, floor(times / interval)) + 1)
@test_nowarn SIR_wrapper!(grid_size, area, param, runparams, abuns)

@test abuns == SIR_wrapper!(grid_size, area, param, runparams, abuns)
