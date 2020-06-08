using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref
using StatsBase
using Distributions

do_plot = false

# Set simulation parameters
age_categories = 10
numclasses = 7
numvirus = 1
birth_rates = fill(0.0/day, numclasses, age_categories)
death_rates = fill(0.0/day, numclasses, age_categories)
birth_rates[:, 2:4] .= uconvert(day^-1, 1/20years); death_rates[1:end-1, :] .= uconvert(day^-1, 1/100years)
virus_growth_asymp = virus_growth_symp = fill(0.1/day, age_categories)
virus_decay = 1.0/day
beta_force = fill(10.0/day, age_categories)
beta_env = fill(10.0/day, age_categories)
ageing = fill(0.0/day, age_categories - 1)# no ageing for now

# Prob of developing symptoms
p_s = fill(0.96, age_categories)
# Prob of hospitalisation
p_h = fill(0.2, age_categories)
# Case fatality ratio
cfr_home = cfr_hospital = fill(0.1, age_categories)
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

param = SEI2HRDGrowth(birth_rates, death_rates, ageing,
                      virus_growth_asymp, virus_growth_symp, virus_decay,
                      beta_force, beta_env, p_s, p_h, cfr_home, cfr_hospital,
                      T_lat, T_asym, T_sym, T_hosp, T_rec)
param = transition(param, age_categories)

# Read in population sizes for Scotland
ukpop = Array{Float64, 2}(readfile(Simulation.path("test", "examples",
                                                   "UK.tif"),
                                   0.0, 7e5, 0, 1.25e6))
# Coarsen grid to 10km
ukpop = [sum(ukpop[i:i+9, j:j+9]) for i in 1:10:size(ukpop, 1), j in 1:10:size(ukpop, 2)]

# Set up simple gridded environment
area = 875_000.0km^2
epienv = simplehabitatAE(298.0K, area, NoControl(), ukpop)

# Set population to initially have no individuals
abun_h = (
    Susceptible = fill(0, age_categories),
    Exposed = fill(0, age_categories),
    Asymptomatic = fill(0, age_categories),
    Symptomatic = fill(0, age_categories),
    Hospitalised = fill(0, age_categories),
    Recovered = fill(0, age_categories),
    Dead = fill(0, age_categories)
)
disease_states = (
    susceptible = ["Susceptible"],
    infectious = ["Asymptomatic", "Symptomatic"]
)
abun_v = (Virus = 0,)

# Dispersal kernels for virus and disease classes
dispersal_dists = fill(1.0km, numclasses * age_categories)
cat_idx = reshape(1:(numclasses * age_categories), age_categories, numclasses)
dispersal_dists[vcat(cat_idx[:, 3:4]...)] .= 20.0km
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = EpiList(traits, abun_v, abun_h, disease_states, movement, param, age_categories)
rel = Gauss{eltype(epienv.habitat)}()

# Create epi system with all information
epi = EpiSystem(epilist, epienv, rel)

# Spread susceptibles randomly over age categories
split_pop = rand.(Multinomial.(Int.(epi.abundances.matrix[1, :]), 10))
for j in 1:size(epi.abundances.matrix, 2)
    epi.abundances.matrix[cat_idx[:, 1], j] .= split_pop[j]
end

# Add in initial infections randomly (samples weighted by population size)
N_cells = size(epi.abundances.matrix, 2)
samp = sample(1:N_cells, weights(1.0 .* epi.abundances.matrix[1, :]), 100)
epi.abundances.matrix[vcat(cat_idx[:, 2]...), samp] .= 10 # Exposed pop

# Run simulation
abuns = zeros(Int64, size(epi.abundances.matrix, 1), N_cells, 366)
times = 1year; interval = 1day; timestep = 1day
@time simulate_record!(abuns, epi, times, interval, timestep)

if do_plot
    using Plots
    # View summed SIR dynamics for whole area
    category_map = (
        "Susceptible" => cat_idx[:, 1],
        "Exposed" => cat_idx[:, 2],
        "Asymptomatic" => cat_idx[:, 3],
        "Symptomatic" => cat_idx[:, 4],
        "Hospital" => cat_idx[:, 5],
        "Recovered" => cat_idx[:, 6],
        "Deaths" => cat_idx[:, 7],
    )
    display(plot_epidynamics(epi, abuns; category_map=category_map))
end
