using Simulation
using SimulationData
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref
using StatsBase
using Distributions
using AxisArrays
using HTTP
using Random
using DataFrames
using Plots

function run_model(api::DataPipelineAPI, times::Unitful.Time, interval::Unitful.Time, timestep::Unitful.Time; do_plot::Bool = false, do_download::Bool = true, save::Bool = false, savepath::String = pwd())
    # Download and read in population sizes for Scotland
    dir = Simulation.path("test", "TEMP")
    file = joinpath(dir, "demographics.h5")
    if do_download
        !isdir(Simulation.path("test", "TEMP")) && mkdir(Simulation.path("test", "TEMP"))
        io = open(Simulation.path("test", "TEMP", "demographics.h5"), "w")
        r = HTTP.request("GET", "https://raw.githubusercontent.com/ScottishCovidResponse/temporary_data/master/human/demographics/scotland/data/demographics.h5")
        write(io, r.body)
        close(io)
    end
    scotpop = parse_hdf5(file, grid = "1km", component = "grid1km/10year/persons")

    # Read number of age categories
    age_categories = size(scotpop, 3)

    # Set up simple gridded environment
    area = (AxisArrays.axes(scotpop, 1)[end] + AxisArrays.axes(scotpop, 1)[2] -
        2 * AxisArrays.axes(scotpop, 1)[1]) *
        (AxisArrays.axes(scotpop, 2)[end] + AxisArrays.axes(scotpop, 2)[2] -
        2 * AxisArrays.axes(scotpop, 2)[1]) * 1.0

    # Sum up age categories and turn into simple matrix
    total_pop = dropdims(sum(Float64.(scotpop), dims=3), dims=3)
    total_pop = AxisArray(total_pop, AxisArrays.axes(scotpop)[1], AxisArrays.axes(scotpop)[2])
    total_pop.data[total_pop .≈ 0.0] .= NaN
    # Shrink to smallest bounding box. The NaNs are inactive.
    total_pop = shrink_to_active(total_pop);

    # Set simulation parameters
    numclasses = 8
    numvirus = age_categories + 1
    birth_rates = fill(0.0/day, numclasses, age_categories)
    death_rates = fill(0.0/day, numclasses, age_categories)
    birth_rates[:, 2:4] .= uconvert(day^-1, 1/20years)
    death_rates[1:end-1, :] .= uconvert(day^-1, 1/100years)
    virus_growth_asymp = virus_growth_presymp = virus_growth_symp = fill(0.1/day, age_categories)
    virus_decay = 1.0/day
    beta_force = fill(10.0/day, age_categories)
    beta_env = fill(10.0/day, age_categories)
    age_mixing = fill(1.0, age_categories, age_categories)

    # Prob of developing symptoms
    p_s = fill(0.96, age_categories)
    # Prob of hospitalisation
    p_h = [0.143, 0.143, 0.1141, 0.117, 0.102, 0.125, 0.2, 0.303, 0.303, 0.303]
    # Case fatality ratio
    cfr_home = cfr_hospital = [0.0, 0.002, 0.002, 0.002, 0.004, 0.013, 0.036, 0.08, 0.148, 0.148]
    # Time exposed
    T_lat = 3days

    # Time asymptomatic
    T_asym = days(read_estimate(
        api,
        "human/infection/SARS-CoV-2/asymptomatic-period",
        "asymptomatic-period"
    )Unitful.hr)
    @show T_asym

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

    param = (birth = birth, death = death, virus_growth = [virus_growth_asymp virus_growth_presymp virus_growth_symp], virus_decay = virus_decay, beta_force = beta_force, beta_env = beta_env, age_mixing = age_mixing)

    epienv = simplehabitatAE(298.0K, size(total_pop), area, Lockdown(20days))

    # Set population to initially have no individuals
    abun_h = (
        Susceptible = fill(0, age_categories),
        Exposed = fill(0, age_categories),
        Asymptomatic = fill(0, age_categories),
        Presymptomatic = fill(0, age_categories),
        Symptomatic = fill(0, age_categories),
        Hospitalised = fill(0, age_categories),
        Recovered = fill(0, age_categories),
        Dead = fill(0, age_categories)
    )

    disease_classes = (
        susceptible = ["Susceptible"],
        infectious = ["Asymptomatic", "Presymptomatic", "Symptomatic"]
    )

    abun_v = (Environment = 0, Force = fill(0, age_categories))

    movement_balance = (home = fill(0.5, numclasses * age_categories), work = fill(0.5, numclasses * age_categories))

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(1.0km, length(total_pop))
    thresholds = fill(1e-3, length(total_pop))
    kernel = GaussianKernel.(dispersal_dists, thresholds)
    home = AlwaysMovement(kernel)

    # Import commuter data (for now, fake table)
    active_cells = findall(.!isnan.(total_pop[1:end]))
    from = active_cells
    to = sample(active_cells, weights(total_pop[active_cells]), length(active_cells))
    count = round.(total_pop[to]/10)
    home_to_work = DataFrame([from, to, count], [:from, :to, :count])
    work = Commuting(home_to_work)
    movement = EpiMovement(home, work)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, paramDat, param, age_categories, movement_balance)
    rel = Gauss{eltype(epienv.habitat)}()

    initial_infecteds = 100
    # Create epi system with all information
    @time epi = EpiSystem(epilist, epienv, rel, total_pop, UInt32(1), initial_infected = initial_infecteds)

    # Populate susceptibles according to actual population spread
    cat_idx = reshape(1:(numclasses * age_categories), age_categories, numclasses)
    reshaped_pop =
        reshape(scotpop[1:size(epienv.active, 1), 1:size(epienv.active, 2), :],
                size(epienv.active, 1) * size(epienv.active, 2), size(scotpop, 3))'
    epi.abundances.matrix[cat_idx[:, 1], :] = reshaped_pop
    N_cells = size(epi.abundances.matrix, 2)

    # Turn off work moves for <20s and >70s
    epi.epilist.human.home_balance[cat_idx[1:2, :]] .= 1.0
    epi.epilist.human.home_balance[cat_idx[7:10, :]] .= 1.0
    epi.epilist.human.work_balance[cat_idx[1:2, :]] .= 0.0
    epi.epilist.human.work_balance[cat_idx[7:10, :]] .= 0.0

    # Run simulation
    abuns = zeros(UInt32, size(epi.abundances.matrix, 1), N_cells, floor(Int, times/timestep) + 1)
    @time simulate_record!(abuns, epi, times, interval, timestep, save = save, save_path = savepath)

    # Write to pipeline
    write_array(api, "simulation-outputs", "final-abundances", DataPipelineArray(abuns))

    if do_plot
        # View summed SIR dynamics for whole area
        category_map = (
            "Susceptible" => cat_idx[:, 1],
            "Exposed" => cat_idx[:, 2],
            "Asymptomatic" => cat_idx[:, 3],
            "Presymptomatic" => cat_idx[:, 4],
            "Symptomatic" => cat_idx[:, 5],
            "Hospital" => cat_idx[:, 6],
            "Recovered" => cat_idx[:, 7],
            "Deaths" => cat_idx[:, 8],
        )
        display(plot_epidynamics(epi, abuns, category_map = category_map))
        display(plot_epiheatmaps(epi, abuns, steps = [30]))
    end
    return abuns
end

config = "data_config.yaml"
download_data_registry(config)
times = 2months; interval = 1day; timestep = 1day
abuns = StandardAPI(config, "test_uri", "test_git_sha") do api
    run_model(api, times, interval, timestep)
end;
