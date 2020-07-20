## Load packages
using Simulation
using Simulation.Units
using Unitful
using Unitful.DefaultSymbols

# Simple SIRGrowth model

## Set model parameters
#=
beta_env = transmission from environmental reservoir
beta_force = transmission from airborne force of infection
sigma = recovery
virus_growth = rate of virus produced per infected
virus_decay = decay rate of environmental reservoir
mean_dispersal_dist = average dispersal distance of virus per each infected
=#
param = (beta_env = 2.0/day,
         beta_force = 3.0/day,
         sigma = 0.02/day,
         virus_growth = 1e-3/day,
         virus_decay = 1e-3/day,
         mean_dispersal_dist = 5.0km)

## Set run parameters:
#=
times = amount of time to simulate over
interval = how often to record output
timestep = timestep of simulation
=#
runparams = (times = 2years,
             interval = 1day,
             timestep = 1day)

## Set grid size & area
grid_size = (4,4)
area = 100.0km^2

##  Run model
# This set up runs in ~0.2 seconds per repeat.
@time abuns = SIR_wrapper(grid_size, area, param, runparams) #1st run: 15sec


## More complex example with current SEI3HRD structure
param = (beta_env = 1.0/day,
         beta_force = 1.0/day,
         virus_growth_symp = 1e-3/day,
         virus_growth_asymp = 1e-3/day,
         virus_growth_presymp = 1e-3/day,
         virus_decay = 1e-3/day,
         mean_dispersal_dist = 10.0km)

runparams = (times = 2years,
             interval = 1day,
             timestep = 1day)

grid_size = (4,4)
area = 100.0km^2
@time abuns = SEI3HRD_wrapper(grid_size, area, param, runparams)


## Again but with age categories - now supply betas and virus growths as vectors for the number of ages
age_cats = 3
param = (beta_env = fill(1.0/day, age_cats),
         beta_force = fill(1.0/day, age_cats),
         virus_growth_symp = fill(1e-3/day, age_cats),
         virus_growth_asymp = fill(1e-3/day, age_cats),
         virus_growth_presymp = fill(1e-3/day, age_cats),
         virus_decay = 1e-3/day,
         mean_dispersal_dist = 10.0km)
runparams = (times = 2years,
             interval = 1day,
             timestep = 1day)
grid_size = (4,4)
area = 100.0km^2
@time abuns = SEI3HRD_wrapper(grid_size, area, param, runparams, age_cats)
