using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref
using StatsBase

do_plot = false

# Set simulation parameters
birth_rates = 1e-5/day; death_rates = birth_rates
birth = [0.0/day; fill(birth_rates, 6); 0.0/day]
death = [0.0/day; fill(death_rates, 6); 0.0/day]
virus_growth = fill(0.0/day, 8)
virus_growth[4:5] .= 0.1/day
virus_decay = fill(0.0/day, 8)
virus_decay[4] = 1.0/day
beta = fill(0.0/day, 8)
beta[3] = 1e-3/day

# Prob of developing symptoms
p_s = 0.96
# Prob of hospitalisation
p_h = 0.2
# Case fatality ratio
cfr = 0.1
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

param = SEI2HRDGrowth(birth, death, virus_growth, virus_decay, beta,
p_s, p_h, cfr, T_lat, T_asym, T_sym, T_hosp, T_rec)
param = transition(param)

# Read in population sizes for Scotland
scotpop = Array{Float64, 2}(readfile("test/examples/ScotlandDensity2011.tif", 0.0, 7e5, 5e5, 1.25e6))
active = Array{Bool, 2}(.!isnan.(scotpop))
scotpop = round.(scotpop)

# Set up simple gridded environment
grid = (700, 750)
area = 525_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, active, NoControl())

# Set population to initially have no individuals
abun = fill(0, 8)

# Dispersal kernels for virus and disease classes
dispersal_dists = [1e-2km; fill(2.0km, 6); 1e-2km] # Virus disperses further than people for now
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, 8), fill(0.1K, 8))
epilist = Simulation.SEI2HRD(traits, abun, movement, param)
rel = Gauss{eltype(epienv.habitat)}()

# Create epi system with all information
epi = EpiSystem(epilist, epienv, rel)
scotpop[isnan.(scotpop)] .= 0
epi.abundances.grid[2, :, :] .+= scotpop

# Add in initial infections randomly (samples weighted by population size)
samp = sample(1:525_000, weights(1.0 .* epi.abundances.matrix[2, :]), 100)
epi.abundances.matrix[1, samp] .= 100 # Virus pop
epi.abundances.matrix[4:5, samp] .= 10 # Inf pop

# Run simulation
abuns = zeros(Int64, 8, 525_000, 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

if do_plot
    using Plots
    plotlyjs()
    # View summed SIR dynamics for whole area
    display(plot(mapslices(sum, abuns[2, :, :], dims = 1)[1, :], color = :Blue, label = "Susceptible"))
    display(plot!(mapslices(sum, abuns[3, :, :], dims = 1)[1, :], color = :Orange, label = "Exposed"))
    display(plot!(mapslices(sum, abuns[4, :, :], dims = 1)[1, :], color = :Yellow, label = "Asymptomatic"))
    display(plot!(mapslices(sum, abuns[5, :, :], dims = 1)[1, :], color = :Red, label = "Symptomatic"))
    display(plot!(mapslices(sum, abuns[6, :, :], dims = 1)[1, :], color = :Darkred, label = "Hospital"))
    display(plot!(mapslices(sum, abuns[7, :, :], dims = 1)[1, :], color = :Green, label = "Recovered"))
    display(plot!(mapslices(sum, abuns[8, :, :], dims = 1)[1, :], color = :Black, label = "Deaths"))
end
