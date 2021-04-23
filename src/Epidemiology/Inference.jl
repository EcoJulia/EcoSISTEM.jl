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
    Ncells = prod(grid_size)

    # Set initial population sizes for all pathogen categories
    virus = 0
    abun_v = DataFrame([
        (name="Environment", initial=virus),
        (name="Force", initial=0),
    ])
    numvirus = nrow(abun_v)

    # Set initial population sizes for all human categories
    susceptible = 500_000 * Ncells
    abun_h = DataFrame([
        (name="Susceptible", type=Susceptible, initial=susceptible),
        (name="Infected", type=Infectious, initial=0),
        (name="Recovered", type=Removed, initial=0),
        (name="Dead", type=Removed, initial=0),
    ])
    numclasses = nrow(abun_h)

    # Set non-pathogen mediated transitions
    sigma = params.sigma
    transitions = DataFrame([
        (from="Infected", to="Recovered", prob=sigma),
    ])

    # Set simulation parameters & create transition matrices
    birth = fill(0.0/day, numclasses)
    death = fill(0.0/day, numclasses)
    param = (birth = birth, death = death, virus_growth = params.virus_growth, virus_decay = params.virus_decay, beta_env = params.beta_env, beta_force = params.beta_force)

    # Set up simple gridded environment
    epienv = simplehabitatAE(298.0K, grid_size, area, NoControl())

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(params.mean_dispersal_dist, Ncells)
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = SpeciesList(traits, abun_v, abun_h, movement, transitions, param)

    # Create epi system with all information
    rel = Gauss{eltype(epienv.habitat)}()
    epi = Ecosystem(epilist, epienv, rel)

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

    Ncells = prod(grid_size)

    # Set initial population sizes for all pathogen categories
    virus = 0
    abun_v = DataFrame([
        (name="Environment", initial=virus),
        (name="Force", initial=fill(0, num_ages)),
    ])
    numvirus = sum(length.(abun_v.initial))

    # Set initial population
    abun_h = DataFrame([
        (name="Susceptible", type=Susceptible, initial=fill(500_000 * Ncells, num_ages)),
        (name="Exposed", type=OtherDiseaseState, initial=fill(0, num_ages)),
        (name="Asymptomatic", type=Infectious, initial=fill(0, num_ages)),
        (name="Presymptomatic", type=Infectious, initial=fill(0, num_ages)),
        (name="Symptomatic", type=Infectious, initial=fill(0, num_ages)),
        (name="Hospitalised", type=OtherDiseaseState, initial=fill(0, num_ages)),
        (name="Recovered", type=Removed, initial=fill(0, num_ages)),
        (name="Dead", type=Removed, initial=fill(0, num_ages)),
    ])
    numclasses = nrow(abun_h)
    numstates = sum(length.(abun_h.initial))
    cat_idx = reshape(1:numstates, num_ages, numclasses)

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
    mu_3 = fill(1 / T_presym, length(params.beta_force))
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

    transitions = DataFrame([
        (from="Exposed", to="Asymptomatic", prob=mu_1),
        (from="Exposed", to="Presymptomatic", prob=mu_2),
        (from="Presymptomatic", to="Symptomatic", prob=mu_3),
        (from="Symptomatic", to="Hospitalised", prob=hospitalisation),
        (from="Asymptomatic", to="Recovered", prob=sigma_1),
        (from="Symptomatic", to="Recovered", prob=sigma_2),
        (from="Hospitalised", to="Recovered", prob=sigma_hospital),
        (from="Symptomatic", to="Dead", prob=death_home),
        (from="Hospitalised", to="Dead", prob=death_hospital),
    ])

    # Set simulation parameters
    birth = fill(0.0/day, numclasses, num_ages)
    death = fill(0.0/day, numclasses, num_ages)
    age_mixing = fill(1.0, num_ages, num_ages)

    @assert length(params.beta_force) == length(params.beta_env) == length(params.virus_growth_asymp) == length(params.virus_growth_presymp) == length(params.virus_growth_symp) == num_ages

    param = (birth = birth, death = death, virus_growth = [params.virus_growth_asymp params.virus_growth_presymp params.virus_growth_symp], virus_decay = params.virus_decay, beta_force = params.beta_force, beta_env = params.beta_env, age_mixing = age_mixing)

    # Set up simple gridded environment
    epienv = simplehabitatAE(298.0K, grid_size, area, NoControl())

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(params.mean_dispersal_dist, Ncells)
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = SpeciesList(traits, abun_v, abun_h, movement, transitions, param, num_ages)

    # Create epi system with all information
    rel = Gauss{eltype(epienv.habitat)}()
    epi = Ecosystem(epilist, epienv, rel)

    # Seed infected category at a single location
    human(epi.abundances)[cat_idx[:, 2], 1] .= 100 * Ncells

    # Run simulation
    times = runtimes.times; interval = runtimes.interval; timestep = runtimes.timestep
    simulate_record!(abuns, epi, times, interval, timestep; save=false)
    return abuns
end
