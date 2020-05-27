"""
    SIR_wrapper(grid_size::Tuple{Int64, Int64}, reps::Int64=1)

Function to simulate simple SIR stochastic realisations, given a grid dimension, `grid_size`, area size, set of parameters, `params` and number of repeats, `reps`.

Outputs an vector of abundance matrices, one for each repeat, of compartment by grid cell. Compartments for the SIR model are: Susceptible, Infected, Recovered, Dead.
"""
function SIR_wrapper(grid_size::Tuple{Int64, Int64}, area::Unitful.Area{Float64}, params::NamedTuple, reps::Int64=1)
    # Set up
    do_save = (@isdefined do_save) ? do_save : false
    save_path = (@isdefined save_path) ? save_path : pwd()
    abuns = Vector{Array{Int64, 3}}(undef, reps)
    numclasses = 4
    Ncells = grid_size[1] * grid_size[2]

    # Extract model params from tuple
    beta_force = params[:beta_force]
    beta_env = params[:beta_env]
    sigma = params[:sigma]
    virus_growth = params[:virus_growth]
    virus_decay = params[:virus_decay]
    mean_dispersal_dist = params[:mean_dispersal_dist]
    # Loop through repeats
    for i in 1:reps
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
            infected = 100 * Ncells,
            recovered = 0,
            dead = 0,
        )
        abun = [initial_pops...]

        # Dispersal kernels for virus and disease classes
        dispersal_dists = fill(sqrt(area/Ncells)/5, numclasses)
        dispersal_dists[2] = mean_dispersal_dist
        kernel = GaussianKernel.(dispersal_dists, 1e-10)
        movement = AlwaysMovement(kernel)

        # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
        traits = GaussTrait(fill(298.0K, numclasses), fill(0.1K, numclasses))
        epilist = SIR(traits, abun, movement, param)

        # Create epi system with all information
        rel = Gauss{eltype(epienv.habitat)}()
        epi = EpiSystem(epilist, epienv, rel)

        # Run simulation
        times = 2years; interval = 1day; timestep = 1day
        abuns[i] = zeros(Int64, size(epi.abundances.matrix, 1), Ncells, convert(Int64, floor(times / interval)) + 1)
        thisabun = abuns[i]
        @time simulate_record!(thisabun, epi, times, interval, timestep;save=do_save, save_path=joinpath(save_path, "rep_$(i)"))
    end
    return abuns
end
