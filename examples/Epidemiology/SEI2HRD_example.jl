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
scotpop = Array{Float64, 2}(readfile("test/examples/ScotlandDensity2011.tif", 0.0, 7e5, 5e5, 1.25e6))

# Set up simple gridded environment
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, area, NoControl(), scotpop)

# Set population to initially have no individuals
abun_h = (
    Susceptible = 0,
    Exposed = 0,
    Asymptomatic = 0,
    Symptomatic = 0,
    Hospitalised = 0,
    Recovered = 0,
    Dead = 0,
)
abun_v = (Virus = 0,)

# Dispersal kernels for virus and disease classes
dispersal_dists = [fill(2.0km, 6); 1e-2km] # Virus disperses further than people for now
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = EpiList(traits, abun_v, abun_h, movement, param)
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
    plotlyjs()
    # View summed SIR dynamics for whole area
    display(plot(mapslices(sum, abuns[1, :, :], dims = 1)[1, :], color = :Blue, label = "Susceptible"))
    display(plot!(mapslices(sum, abuns[2, :, :], dims = 1)[1, :], color = :Orange, label = "Exposed"))
    display(plot!(mapslices(sum, abuns[3, :, :], dims = 1)[1, :], color = :Yellow, label = "Asymptomatic"))
    display(plot!(mapslices(sum, abuns[4, :, :], dims = 1)[1, :], color = :Red, label = "Symptomatic"))
    display(plot!(mapslices(sum, abuns[5, :, :], dims = 1)[1, :], color = :Darkred, label = "Hospital"))
    display(plot!(mapslices(sum, abuns[6, :, :], dims = 1)[1, :], color = :Green, label = "Recovered"))
    display(plot!(mapslices(sum, abuns[7, :, :], dims = 1)[1, :], color = :Black, label = "Deaths"))
end
