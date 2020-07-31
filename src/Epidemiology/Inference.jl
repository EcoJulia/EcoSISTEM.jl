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
    param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force)
    paramDat = DataFrame([["Infected"], [ "Recovered"], [sigma]], [:from, :to, :prob])

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
    dispersal_dists = fill(mean_dispersal_dist, Ncells)
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, paramDat, param)

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
function SEI3HRD_wrapper(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple, num_ages::Int64 = 1)
    times = runtimes.times; interval = runtimes.interval
    Ncells = grid_size[1] * grid_size[2]
    numclasses = 8
    abuns = zeros(Int64, numclasses * num_ages, Ncells, convert(Int64, floor(times / interval)) + 1)
    return SEI3HRD_wrapper!(grid_size, area, params, runtimes, abuns, num_ages)
end
"""
    SEI3HRD_wrapper!(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple, abuns::Array{Int64, 3})

Function to simulate SEI3HRD stochastic realisations, given a grid dimension, `grid_size`, area size, a set of model parameters, `params`, and running parameters, `runtimes`.

Fills an abundance matrice of compartment by grid cell over time. Compartments for the SIR model are: Susceptible, Exposed, Asymptomatic Infected, Pre-symptomatic Infected, Symptomatic Infected, Hospitalised, Recovered, Dead.
"""
function SEI3HRD_wrapper!(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, runtimes::NamedTuple, abuns::Array{Int64, 3}, num_ages::Int64 = 1)
    # Set up
    numclasses = 8
    numvirus = num_ages + 1
    Ncells = grid_size[1] * grid_size[2]
    cat_idx = reshape(1:(numclasses * num_ages), num_ages, numclasses)

    # Extract model params from tuple
    beta_force = params.beta_force
    beta_env = params.beta_env
    virus_growth_asymp = params.virus_growth_asymp
    virus_growth_presymp = params.virus_growth_presymp
    virus_growth_symp = params.virus_growth_symp
    virus_decay = params.virus_decay
    mean_dispersal_dist = params.mean_dispersal_dist

    (length(beta_force) == num_ages) && (length(beta_env) == num_ages) || error("Length of beta parameters must equal number of age classes")

    (length(virus_growth_asymp) == num_ages) && (length(virus_growth_presymp) == num_ages) && (length(virus_growth_symp) == num_ages) || error("Length of virus growth parameters must equal number of age classes")

    # Set simulation parameters & create transition matrices
    birth = fill(0.0/day, numclasses, num_ages)
    death = fill(0.0/day, numclasses, num_ages)
    age_mixing = fill(1.0, num_ages, num_ages)

    # Prob of developing symptoms
    p_s = fill(0.96, num_ages)
    # Prob of hospitalisation
    p_h = fill(0.2, num_ages)
    # Case fatality ratio
    cfr_home = cfr_hospital = fill(0.1, num_ages)
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

    # Exposed -> asymptomatic
    mu_1 = (1 .- p_s) .* 1/T_lat
    # Exposed -> Pre-symptomatic
    mu_2 = p_s .* 1/T_lat
    # Pre-symptomatic -> symptomatic
    mu_3 = fill(1 / T_presym, length(beta_force))
    # Symptomatic -> hospital
    hospitalisation = p_h .* 1/T_sym
    # Asymptomatic -> recovered
    sigma_1 = (1 .- p_s) .* 1/T_asym
    # Symptomatic -> recovered
    sigma_2 = (1 .- p_h) .* (1 .- cfr_home) .* 1/T_rec
    # Hospital -> recovered
    sigma_hospital = (1 .- cfr_hospital) .* 1/T_hosp
    # Symptomatic -> death
    death_home = cfr_home .* 2/T_hosp
    # Hospital -> death
    death_hospital = cfr_hospital .* 1/T_hosp

    paramDat = DataFrame([["Exposed", "Exposed", "Presymptomatic", "Symptomatic", "Asymptomatic", "Symptomatic", "Hospitalised", "Symptomatic", "Hospitalised"], ["Asymptomatic", "Presymptomatic", "Symptomatic", "Hospitalised", "Recovered", "Recovered", "Recovered", "Dead", "Dead"], [mu_1, mu_2, mu_3, hospitalisation, sigma_1, sigma_2, sigma_hospital, death_home, death_hospital]], [:from, :to, :prob])

    param = (birth = birth, death = death, virus_growth = [virus_growth_asymp virus_growth_presymp virus_growth_symp], virus_decay = virus_decay, beta_force = beta_force, beta_env = beta_env)

    # Set up simple gridded environment
    epienv = simplehabitatAE(298.0K, grid_size, area, NoControl())

    # Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
    abun_h = (
        Susceptible = fill(500_000 * Ncells, num_ages),
        Exposed = fill(0, num_ages),
        Asymptomatic = fill(0, num_ages),
        Presymptomatic = fill(0, num_ages),
        Symptomatic = fill(0, num_ages),
        Hospitalised = fill(0, num_ages),
        Recovered = fill(0, num_ages),
        Dead = fill(0, num_ages)
    )
    disease_classes = (
        susceptible = ["Susceptible"],
        infectious = ["Asymptomatic", "Presymptomatic", "Symptomatic"]
    )
    abun_v = (Environment = 0, Force = fill(0, num_ages))

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(mean_dispersal_dist, Ncells)
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, paramDat, param, num_ages)

    # Create epi system with all information
    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel)

    # Seed infected category at a single location
    human(epi.abundances)[cat_idx[:, 2], 1] .= 100 * Ncells

    # Run simulation
    times = runtimes.times; interval = runtimes.interval; timestep = runtimes.timestep
    simulate_record!(abuns, epi, times, interval, timestep; save=false)
    return abuns
end
