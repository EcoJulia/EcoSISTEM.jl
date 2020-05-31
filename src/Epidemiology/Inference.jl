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
    numvirus = 1
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
    initial_pops = (
        virus = 0,
        susceptible = 500_000 * Ncells,
        infected = 0,
        recovered = 0,
        dead = 0,
    )
    abun = [initial_pops...]
    abun_h = abun[2:end]
    abun_v = [abun[1]]

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(sqrt(area/Ncells)/5, numclasses)
    dispersal_dists[2] = mean_dispersal_dist
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = AlwaysMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = SIR(traits, abun_v, abun_h, movement, param)

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
