using Simulation
using Simulation.Units
using Unitful
using Unitful.DefaultSymbols

# Set model parameters:
# beta_env = transmission from environmental reservoir
# beta_force = transmission from airborne force of infection
# sigma = recovery
# virus_growth = rate of virus produced per infected
# virus_decay = decay rate of environmental reservoir
# mean_dispersal_dist = average dispersal distance of virus per each infected
param = (beta_env = 1.0/day, beta_force = 1.0/day, sigma = 0.02/day, virus_growth = 1e-3/day, virus_decay = 1e-3/day, mean_dispersal_dist = 10.0km)
# Set grid size & area
grid_size = (4,4)
area = 100.0km^2
repeats = 10
abuns = SIR_wrapper(grid_size, area, param, repeats)
# This set up runs in ~0.14 seconds per repeat.
