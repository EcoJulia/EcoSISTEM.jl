using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref
using StatsBase
using Plots
using Random
using Distributions

const stochasticmode = true
const seed = hash(time()) # seed used for Random.jl and therefore rngs used in Simulation.jl
const ic_rng = Random.MersenneTwister(0) # rng for initial conditions in this file

Random.seed!(seed)

function run_model(times::Unitful.Time, interval::Unitful.Time, timestep::Unitful.Time, do_plot::Bool = false)

    numvirus = 2
    numclasses = 8

    # Set simulation parameters
    birth_rates = 1e-5/day; death_rates = birth_rates
    birth = fill(birth_rates, numclasses)
    death = fill(death_rates, numclasses)
    virus_growth_asymp = virus_growth_presymp = virus_growth_symp = 0.1/day
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
    T_lat = 3days
    # Time asymptomatic
    T_asym = 5days
    # Time pre-symptomatic
    T_presym = 1.5days
    # Time symptomatic
    T_sym = 5days
    # Time in hospital
    T_hosp = 5days
    # Time to recovery if symptomatic
    T_rec = 11days

    param = SEI3HRDGrowth(birth, death, virus_growth_asymp, virus_growth_presymp, virus_growth_symp, virus_decay, beta_force, beta_env, p_s, p_h, cfr_home, cfr_hospital, T_lat, T_asym, T_presym, T_sym, T_hosp, T_rec)
    param = transition(param)

    # Read in population sizes for Scotland
    scotpop = Array{Float64, 2}(readfile(Simulation.path("test", "examples", "ScotlandDensity2011.tif"), 0.0, 7e5, 5e5, 1.25e6))

    # Set up simple gridded environment
    area = 525_000.0km^2
    epienv = simplehabitatAE(298.0K, size(scotpop), area, NoControl())

    # Set population to initially have no individuals
    abun_h = (
        Susceptible = 0,
        Exposed = 0,
        Asymptomatic = 0,
        Presymptomatic = 0,
        Symptomatic = 0,
        Hospitalised = 0,
        Recovered = 0,
        Dead = 0
    )
    disease_classes = (
        susceptible = ["Susceptible"],
        infectious = ["Asymptomatic", "Presymptomatic", "Symptomatic"]
    )
    abun_v = (Environment = 0, Force = 0)

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(2.0km, numclasses)
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, param)
    rel = Gauss{eltype(epienv.habitat)}()

    # multiple dispatch in action
    rngtype = stochasticmode ? Random.MersenneTwister : MedianGenerator

    # Create epi system with all information
    epi = EpiSystem(epilist, epienv, rel, scotpop; rngtype=rngtype)

    # Add in initial infections randomly (samples weighted by population size)
    N_cells = size(epi.abundances.matrix, 2)
    samp = sample(ic_rng, 1:N_cells, weights(1.0 .* human(epi.abundances)[1, :]), 100)
    virus(epi.abundances)[1, samp] .= 100 # Virus pop in Environment
    human(epi.abundances)[2, samp] .= 10 # Exposed pop

    # Run simulation
    abuns = zeros(Int64, numclasses, N_cells, 366)
    times = 2months; interval = 1day; timestep = 1day

    @time simulate_record!(abuns, epi, times, interval, timestep)

    if do_plot
        # View summed SIR dynamics for whole area
        display(plot_epidynamics(epi, abuns))
    end
end

times = 1year; interval = 1day; timestep = 1day
abuns = run_model(times, interval, timestep);
