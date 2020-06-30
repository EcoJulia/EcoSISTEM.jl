"""
    SIR_wrapper(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple)

Function to simulate simple SIR stochastic realisations, given a grid dimension, `grid_size`, area size, a set of model parameters, `params`, and running parameters, `runtimes`.

Outputs an abundance matrice of compartment by grid cell over time. Compartments for the SIR model are: Susceptible, Infected, Recovered, Dead.
"""
function SIR_wrapper(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple)
    times = runtimes.times; interval = runtimes.interval
    Ncells = grid_size[1] * grid_size[2]
    numclasses = 4
    abuns = zeros(Int64, numclasses, Ncells, convert(Int64, floor(times / interval)) + 1)
    return SIR_wrapper!(grid_size, area, params, runtimes, abuns)
end
"""
    SIR_wrapper!(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple, abuns::Array{Int64, 3})

Function to simulate simple SIR stochastic realisations, given a grid dimension, `grid_size`, area size, a set of model parameters, `params`, and running parameters, `runtimes`.

Fills an abundance matrice of compartment by grid cell over time. Compartments for the SIR model are: Susceptible, Infected, Recovered, Dead.
"""
function SIR_wrapper!(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple, abuns::Array{Int64, 3})
    # Set up
    numclasses = 4
    numvirus = 2
    Ncells = grid_size[1] * grid_size[2]

    # Extract model params from tuple
    beta_force = params.beta_force
    beta_env = params.beta_env
    sigma = params.sigma
    virus_growth = params.virus_growth
    virus_decay = params.virus_decay
    mean_dispersal_dist = params.mean_dispersal_dist

    # Set simulation parameters & create transition matrices
    birth = fill(0.0/day, numclasses)
    death = fill(0.0/day, numclasses)
    param = SIRGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
    param = transition(param)

    # Set up simple gridded environment
    epienv = simplehabitatAE(298.0K, grid_size, area, NoControl())

    # Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
    abun_h = (
        Susceptible = 500_000 * Ncells,
        Infected = 0,
        Recovered = 0,
        Dead = 0
    )
    disease_classes = (
        susceptible = ["Susceptible"],
        infectious = ["Infected"]
    )
    abun_v = (Environment = 0, Force = 0)

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(sqrt(area/Ncells)/5, numclasses)
    dispersal_dists[2] = mean_dispersal_dist
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = AlwaysMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, param)

    # Create epi system with all information
    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel)

    # Seed infected category at a single location
    human(epi.abundances)[2, 1] = 100 * Ncells

    # Run simulation
    times = runtimes.times; interval = runtimes.interval; timestep = runtimes.timestep
    simulate_record!(abuns, epi, times, interval, timestep; save=false)
    return abuns
end


"""
    SEI3HRD_wrapper(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple)

Function to simulate SEI3HRD stochastic realisations, given a grid dimension, `grid_size`, area size, a set of model parameters, `params`, and running parameters, `runtimes`.

Outputs an abundance matrice of compartment by grid cell over time. Compartments for the SEI3HRD model are: Susceptible, Exposed, Asymptomatic Infected, Pre-symptomatic Infected, Symptomatic Infected, Hospitalised, Recovered, Dead.
"""
function SEI3HRD_wrapper(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple)
    times = runtimes.times; interval = runtimes.interval
    Ncells = grid_size[1] * grid_size[2]
    numclasses = 8
    abuns = zeros(Int64, numclasses, Ncells, convert(Int64, floor(times / interval)) + 1)
    return SEI3HRD_wrapper!(grid_size, area, params, runtimes, abuns)
end
"""
    SEI3HRD_wrapper!(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple, abuns::Array{Int64, 3})

Function to simulate SEI3HRD stochastic realisations, given a grid dimension, `grid_size`, area size, a set of model parameters, `params`, and running parameters, `runtimes`.

Fills an abundance matrice of compartment by grid cell over time. Compartments for the SIR model are: Susceptible, Exposed, Asymptomatic Infected, Pre-symptomatic Infected, Symptomatic Infected, Hospitalised, Recovered, Dead.
"""
function SEI3HRD_wrapper!(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple, abuns::Array{Int64, 3})
    # Set up
    numclasses = 8
    numvirus = 2
    Ncells = grid_size[1] * grid_size[2]

    # Extract model params from tuple
    beta_force = params.beta_force
    beta_env = params.beta_env
    virus_growth_asymp = params.virus_growth_asymp
    virus_growth_presymp = params.virus_growth_presymp
    virus_growth_symp = params.virus_growth_symp
    virus_decay = params.virus_decay
    mean_dispersal_dist = params.mean_dispersal_dist

    # Set simulation parameters & create transition matrices
    birth = fill(0.0/day, numclasses)
    death = fill(0.0/day, numclasses)
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

    # Set up simple gridded environment
    epienv = simplehabitatAE(298.0K, grid_size, area, NoControl())

    # Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
    abun_h = (
        Susceptible = 500_000 * Ncells,
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
    dispersal_dists = fill(sqrt(area/Ncells)/5, numclasses)
    dispersal_dists[2] = mean_dispersal_dist
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = AlwaysMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, param)

    # Create epi system with all information
    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel)

    # Seed infected category at a single location
    human(epi.abundances)[2, 1] = 100 * Ncells

    # Run simulation
    times = runtimes.times; interval = runtimes.interval; timestep = runtimes.timestep
    simulate_record!(abuns, epi, times, interval, timestep; save=false)
    return abuns
end
