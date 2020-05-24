using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref
using StatsBase
using Test
using Simulation: human, virus

# sort out settings to potentially save inputs/outputs of `simulate`
do_save = (@isdefined do_save) ? do_save : false
save_path = (@isdefined save_path) ? save_path : pwd()

## High transmission & 100% case fatality
##
# Set simulation parameters
numclasses = 7
birth = fill(0.0/day, numclasses)
death = fill(0.0/day, numclasses)
virus_growth_asymp = virus_growth_symp = 1e-3/day
virus_decay = 1/3days
beta_force = 1e3/day
beta_env = 1e3/day

# Prob of developing symptoms
p_s = 1.0
# Prob of hospitalisation
p_h = 0.2
# Case fatality ratio
cfr_home = cfr_hosp = 1.0
# Time exposed
T_lat = 5days
# Time asymptomatic
T_asym = 3days
# Time symptomatic
T_sym = 5days
# Time in hospital
T_hosp = 5days
# Time to recovery if symptomatic
T_rec = 11days

param = SEI2HRDGrowth(birth, death, virus_growth_asymp, virus_growth_symp, virus_decay, beta_force, beta_env, p_s, p_h, cfr_home, cfr_hosp, T_lat, T_asym, T_sym, T_hosp, T_rec)
param = transition(param)

# Set up simple gridded environment
grid = (4, 4)
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set population to initially have no individuals
initial_pops = (
  virus = 0,
  susceptible = 5_000_000,
  exposed = 0,
  asymptomatic = 0,
  symptomatic = 0,
  hospitalised = 0,
  recovered = 0,
  dead = 0,
)
abun = [initial_pops...]

# Dispersal kernels for virus dispersal from different disease classes
dispersal_dists = fill(1.0km, numclasses)
dispersal_dists[3:4] .= 500.0km

kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numclasses), fill(0.1K, numclasses))
epilist = Simulation.SEI2HRD(traits, abun, movement, param)
rel = Gauss{eltype(epienv.habitat)}()

# Create epi system with all information
new_symptomatic = new_asymptomatic = 100
new_virus = 1000
epi = EpiSystem(epilist, epienv, rel)
human(epi.abundances)[1, 1] = new_virus
human(epi.abundances)[4:5, 1] .= new_symptomatic

# Run simulation
abuns = zeros(Int64, size(epi.abundances.matrix, 1), size(epi.abundances.matrix, 2), 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "high_trans"))


# Test everyone becomes infected and dies
@test sum(abuns[2, :, 365]) == 0
@test sum(abuns[end, :, 365]) == (initial_pops[:susceptible] + new_symptomatic + new_asymptomatic)

## Low transmission & 100% case fatality
##

birth = fill(0.0/day, numclasses)
death = fill(0.0/day, numclasses)
virus_growth_asymp = virus_growth_symp = 1e-3/day
virus_decay = 1.0/day
beta_force = 1e-10/day
beta_env = 1e-10/day

# Prob of developing symptoms
p_s = 1.0
# Prob of hospitalisation
p_h = 0.2
# Case fatality ratio
cfr_home = cfr_hosp = 1.0
# Time exposed
T_lat = 5days
# Time asymptomatic
T_asym = 3days
# Time symptomatic
T_sym = 5days
# Time in hospital
T_hosp = 5days
# Time to recovery if symptomatic
T_rec = 11days

param = SEI2HRDGrowth(birth, death, virus_growth_asymp, virus_growth_symp, virus_decay, beta_force, beta_env, p_s, p_h, cfr_home, cfr_hosp, T_lat, T_asym, T_sym, T_hosp, T_rec)
param = transition(param)

# Set up simple gridded environment
grid = (8, 8)
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set population to initially have no individuals
initial_pops = (
  virus = 0,
  susceptible = 500_000 * prod(grid),
  exposed = 0,
  asymptomatic = 0,
  symptomatic = 0,
  hospitalised = 0,
  recovered = 0,
  dead = 5,
)

abun = [initial_pops...]

# Dispersal kernels for virus dispersal from different disease classes
dispersal_dists = fill(1.0km, numclasses)
dispersal_dists[3:4] .= 200.0km

kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numclasses), fill(0.1K, numclasses))
epilist = Simulation.SEI2HRD(traits, abun, movement, param)
rel = Gauss{eltype(epienv.habitat)}()

# Create epi system with all information
new_symptomatic = new_asymptomatic = 100
new_virus = 1000
epi = EpiSystem(epilist, epienv, rel)
virus(epi.abundances)[1, 1] = new_virus
human(epi.abundances)[4:5, 1] .= new_symptomatic

# Run simulation
times = 1year; interval = 1day; timestep = 1day
abuns = zeros(Int64, size(epi.abundances.matrix, 1), size(epi.abundances.matrix, 2),
              convert(Int64, floor(times / interval)) + 1)
@time simulate_record!(abuns, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "low_trans"))



# Test no one becomes infected & dies # TODO: Check this comment
@test sum(abuns[2, :, 365]) == initial_pops[:susceptible]
@test sum(abuns[end, :, 365]) == (new_symptomatic + new_asymptomatic) + initial_pops[:dead]

### TEST OUTPUTS
# TODO: When shifting out virus from Epilist, these indexes will need updating.

idx_sus = 2
idx_rec = 7
idx_dead = 8

# Test susceptible population decreasing or constant only [Source]
@test all(diff(vec(sum(abuns[idx_sus, :, :], dims = 1))) .<= 0)
@test sum(abuns[idx_sus, :, 1]) == initial_pops[:susceptible]

# Test recovered population increasing  or constant only [Sink]
@test all(diff(vec(sum(abuns[idx_rec, :, :], dims = 1))) .>= 0)
@test sum(abuns[idx_rec, :, 1]) == initial_pops[:recovered]

# Test dead population increasing or constant only [Sink]
@test all(diff(vec(sum(abuns[idx_dead, :, :], dims = 1))) .>= 0)
@test sum(abuns[idx_dead, :, 1]) == initial_pops[:dead]
