using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref
using StatsBase

do_plot = false

numvirus = 1
numclasses = 7

# Set simulation parameters
birth_rates = 1e-5/day; death_rates = birth_rates
birth = [fill(birth_rates, 6); 0.0/day]
death = [fill(death_rates, 6); 0.0/day]
virus_growth_asymp = virus_growth_symp = 0.1/day
virus_decay = 1.0/day
beta_force = 10.0/day
beta_env = 10.0/day

# Prob of developing symptoms
p_s = 0.96
# Prob of hospitalisation
p_h = 0.2
# Case fatality ratio
cfr_home = cfr_hospital = 0.1
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

param = SEI2HRDGrowth(birth, death, virus_growth_asymp, virus_growth_symp, virus_decay, beta_force, beta_env, p_s, p_h, cfr_home, cfr_hospital, T_lat, T_asym, T_sym, T_hosp, T_rec)
param = transition(param)

# Read in population sizes for Scotland
scotpop = Array{Float64, 2}(readfile(Simulation.path("test", "examples", "ScotlandDensity2011.tif"), 0.0, 7e5, 5e5, 1.25e6))

# Set up simple gridded environment
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, area, NoControl(), scotpop)

# Set population to initially have no individuals
abun_h = fill(0, numclasses)
abun_v = fill(0, numvirus)

# Dispersal kernels for virus and disease classes
dispersal_dists = [fill(2.0km, 6); 1e-2km] # Virus disperses further than people for now
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = Simulation.SEI2HRD(traits, abun_v, abun_h, movement, param)
rel = Gauss{eltype(epienv.habitat)}()

# Create epi system with all information
epi = EpiSystem(epilist, epienv, rel)

# Add in initial infections randomly (samples weighted by population size)
N_cells = size(epi.abundances.matrix, 2)
samp = sample(1:N_cells, weights(1.0 .* human(epi.abundances)[1, :]), 100)
virus(epi.abundances)[1, samp] .= 100 # Virus pop
human(epi.abundances)[2, samp] .= 10 # Exposed pop

# Run simulation
abuns = zeros(Int64, numclasses, N_cells, 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

if do_plot
    using Plots
    # View summed SIR dynamics for whole area
    display(plot_epidynamics(epi, abuns))
end
