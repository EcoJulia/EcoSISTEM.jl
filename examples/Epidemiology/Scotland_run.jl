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
    scotpop = parse_scottish_population(api)

    # Read number of age categories
    age_categories = size(scotpop, 3)

    # Set initial population sizes for all pathogen categories
    abun_v = DataFrame([
        (name="Environment", initial=0),
        (name="Force", initial=fill(0, age_categories)),
    ])
    numvirus = sum(length.(abun_v.initial))

    # Set population to initially have no individuals
    abun_h = DataFrame([
        (name="Susceptible", type=Susceptible, initial=fill(0, age_categories)),
        (name="Exposed", type=OtherDiseaseState, initial=fill(0, age_categories)),
        (name="Asymptomatic", type=Infectious, initial=fill(0, age_categories)),
        (name="Presymptomatic", type=Infectious, initial=fill(0, age_categories)),
        (name="Symptomatic", type=Infectious, initial=fill(0, age_categories)),
        (name="Hospitalised", type=OtherDiseaseState, initial=fill(0, age_categories)),
        (name="Recovered", type=Removed, initial=fill(0, age_categories)),
        (name="Dead", type=Removed, initial=fill(0, age_categories)),
    ])
    numclasses = nrow(abun_h)
    numstates = sum(length.(abun_h.initial))

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

    # Prob of developing symptoms
    p_s = fill(read_estimate(
               api,
               "human/infection/SARS-CoV-2/symptom-probability",
               "symptom-probability"
           ), age_categories)

    param_tab = read_table(api, "prob_hosp_and_cfr/data_for_scotland", "cfr_byage")
    # Prob of hospitalisation
    p_h = param_tab.p_h[1:end-1] # remove HCW
    pushfirst!(p_h, p_h[1]) # extend age categories
    append!(p_h, fill(p_h[end], 2)) # extend age categories
    # Case fatality ratio
    cfr_home = param_tab.cfr[1:end-1]
    pushfirst!(cfr_home, cfr_home[1])
    append!(cfr_home, fill(cfr_home[end], 2))
    cfr_hospital = param_tab.p_d[1:end-1]
    pushfirst!(cfr_hospital, cfr_hospital[1])
    append!(cfr_hospital, fill(cfr_hospital[end], 2))

    @assert length(p_s) == length(p_h) == length(cfr_home)

    # Time exposed
    T_lat = days(read_estimate(
        api,
        "human/infection/SARS-CoV-2/latent-period",
        "latent-period"
    )Unitful.hr)

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
    T_sym = days(read_estimate(
        api,
        "human/infection/SARS-CoV-2/infectious-duration",
        "infectious-duration"
    )Unitful.hr) - T_presym
    # Time in hospital
    T_hosp = read_estimate(
        api,
        "fixed-parameters/T_hos",
        "T_hos"
    )days
    # Time to recovery if symptomatic
    T_rec = read_estimate(
        api,
        "fixed-parameters/T_rec",
        "T_rec"
    )days

    # Exposed -> asymptomatic
    mu_1 = (1 .- p_s) .* 1/T_lat
    # Exposed -> Pre-symptomatic
    mu_2 = p_s .* 1/T_lat
    # Pre-symptomatic -> symptomatic
    mu_3 = fill(1 / T_presym, age_categories)
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
        (from="Hospitalised", to="Dead", prob=death_hospital)
    ])

    # Set simulation parameters
    birth_rates = fill(0.0/day, numclasses, age_categories)
    death_rates = fill(0.0/day, numclasses, age_categories)
    birth_rates[:, 2:4] .= uconvert(day^-1, 1/20years)
    death_rates[1:end-1, :] .= uconvert(day^-1, 1/100years)
    virus_growth_asymp = virus_growth_presymp = virus_growth_symp = fill(0.1/day, age_categories)
    virus_decay = 1.0/day
    beta_force = fill(10.0/day, age_categories)
    beta_env = fill(10.0/day, age_categories)
    age_mixing = fill(1.0, age_categories, age_categories)

    param = (birth = birth_rates, death = death_rates, virus_growth = [virus_growth_asymp virus_growth_presymp virus_growth_symp], virus_decay = virus_decay, beta_force = beta_force, beta_env = beta_env, age_mixing = age_mixing)

    epienv = simplehabitatAE(298.0K, size(total_pop), area, Lockdown(20days))

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
    home_to_work = DataFrame(from=from, to=to, count=count)
    work = Commuting(home_to_work)
    movement = EpiMovement(home, work)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param, age_categories, movement_balance)
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
    #write_array(api, "simulation-outputs", "final-abundances", DataPipelineArray(abuns))

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
