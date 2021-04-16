using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using StatsBase
using Test
using DataFrames
using EcoSISTEM: human, virus

@testset "SEI2HRD" begin

# sort out settings to potentially save inputs/outputs of `simulate`
do_save = (@isdefined do_save) ? do_save : false
save_path = (@isdefined save_path) ? save_path : pwd()

# Set up simple gridded environment
grid = (4, 4)
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all pathogen categories
abun_v = DataFrame([
    (name="Environment", initial=0),
    (name="Force", initial=0),
])
numvirus = nrow(abun_v)

# Set initial population sizes for all human categories
susceptible = 500_000 * prod(grid)
abun_h = DataFrame([
  (name="Susceptible", type=Susceptible, initial=susceptible),
  (name="Exposed", type=OtherDiseaseState, initial=0),
  (name="Asymptomatic", type=Infectious, initial=0),
  (name="Symptomatic", type=Infectious, initial=0),
  (name="Hospitalised", type=OtherDiseaseState, initial=0),
  (name="Recovered", type=Removed, initial=0),
  (name="Dead", type=Removed, initial=0),
])
numclasses = nrow(abun_h)

# Set non-pathogen mediated transitions
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
# Exposed -> asymptomatic
mu_1 = 1/T_lat
# Asymptomatic -> symptomatic
mu_2 = p_s * 1/T_asym
# Symptomatic -> hospital
hospitalisation = p_h * 1/T_sym
# Asymptomatic -> recovered
sigma_1 = (1 - p_s) * 1/T_asym
# Symptomatic -> recovered
sigma_2 = (1 - p_h) * (1 - cfr_home) * 1/T_rec
# Hospital -> recovered
sigma_hospital = (1 - cfr_hosp) * 1/T_hosp
# Symptomatic -> death
death_home = cfr_home * 2/T_hosp
# Hospital -> death
death_hospital = cfr_hosp * 1/T_hosp

transitions = DataFrame([
  (from="Exposed", to="Asymptomatic", prob=mu_1),
  (from="Asymptomatic", to="Symptomatic", prob=mu_2),
  (from="Symptomatic", to="Hospitalised", prob=hospitalisation),
  (from="Asymptomatic", to="Recovered", prob=sigma_1),
  (from="Symptomatic", to="Recovered", prob=sigma_2),
  (from="Hospitalised", to="Recovered", prob=sigma_hospital),
  (from="Symptomatic", to="Dead", prob=death_home),
  (from="Hospitalised", to="Dead", prob=death_hospital),
])

## High transmission & 100% case fatality
birth = fill(0.0/day, numclasses)
death = fill(0.0/day, numclasses)
virus_growth_asymp = virus_growth_symp = 1e-3/day
virus_decay = 1/3days
beta_force = 1e3/day
beta_env = 1e3/day

param = (birth = birth, death = death, virus_growth = [virus_growth_asymp virus_growth_symp], virus_decay = virus_decay, beta_force = beta_force, beta_env = beta_env)

# Dispersal kernels for virus dispersal from different disease classes
dispersal_dists = fill(500.0km, prod(grid))

kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = EpiMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param)
rel = Gauss{eltype(epienv.habitat)}()

# Create epi system with all information
new_symptomatic = new_asymptomatic = 100
new_virus = 1_000
epi = EpiSystem(epilist, epienv, rel)
virus(epi.abundances)[1, 1] = new_virus
human(epi.abundances)[3:4, 1] .= new_symptomatic

# Run simulation
abuns = zeros(Int64, size(epi.abundances.matrix, 1), size(epi.abundances.matrix, 2), 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "high_trans"))

# Test everyone becomes infected and dies
@test sum(abuns[1, :, 365]) == 0
@test sum(abuns[end, :, 365]) == (susceptible + new_symptomatic + new_asymptomatic)

## Low transmission & 100% case fatality
virus_decay = 1.0/day
beta_force = 1e-10/day
beta_env = 1e-10/day

param = (birth = birth, death = death, virus_growth = [virus_growth_asymp virus_growth_symp], virus_decay = virus_decay, beta_force = beta_force, beta_env = beta_env)

# Set up simple gridded environment
grid = (8, 8)
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all human categories
susceptible = 500_000 * prod(grid)
abun_h = DataFrame([
  (name="Susceptible", type=Susceptible, initial=susceptible),
  (name="Exposed", type=OtherDiseaseState, initial=0),
  (name="Asymptomatic", type=Infectious, initial=0),
  (name="Symptomatic", type=Infectious, initial=0),
  (name="Hospitalised", type=OtherDiseaseState, initial=0),
  (name="Recovered", type=Removed, initial=0),
  (name="Dead", type=Removed, initial=0),
])

# Dispersal kernels for virus dispersal from different disease classes
dispersal_dists = fill(200.0km, prod(grid))

kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = EpiMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param)

# Create epi system with all information
new_symptomatic = new_asymptomatic = 100
new_virus = 1000
epi = EpiSystem(epilist, epienv, rel)
virus(epi.abundances)[1, 1] = new_virus
human(epi.abundances)[3:4, 1] .= new_symptomatic

# Run simulation
times = 1year; interval = 1day; timestep = 1day
abuns = zeros(Int64, numclasses, prod(grid), div(times, interval) + 1)
@time simulate_record!(abuns, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "low_trans"))

# Find correct indices in arrays
idx_sus = findfirst(==("Susceptible"), abun_h.name)
idx_rec = findfirst(==("Recovered"), abun_h.name)
idx_dead = findfirst(==("Dead"), abun_h.name)

# Test no one becomes infected & dies # TODO: Check this comment
@test sum(abuns[idx_sus, :, end]) == abun_h.initial[idx_sus]
@test sum(abuns[idx_dead, :, end]) == (new_symptomatic + new_asymptomatic) + abun_h.initial[idx_dead]

### TEST OUTPUTS

# Test susceptible population decreasing or constant only [Source]
@test all(diff(vec(sum(abuns[idx_sus, :, :], dims = 1))) .<= 0)
@test sum(abuns[idx_sus, :, 1]) == abun_h.initial[idx_sus]

# Test recovered population increasing  or constant only [Sink]
@test all(diff(vec(sum(abuns[idx_rec, :, :], dims = 1))) .>= 0)
@test sum(abuns[idx_rec, :, 1]) == abun_h.initial[idx_rec]

# Test dead population increasing or constant only [Sink]
@test all(diff(vec(sum(abuns[idx_dead, :, :], dims = 1))) .>= 0)
@test sum(abuns[idx_dead, :, 1]) == abun_h.initial[idx_dead]
end
