using EcoSISTEM
using DataRegistryUtils
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using StatsBase
using Distributions
using AxisArrays
using HTTP
using Random
using DataFrames
using Plots
using SQLite

function run_model(db::SQLite.DB, times::Unitful.Time, interval::Unitful.Time, timestep::Unitful.Time; do_plot::Bool = false, do_download::Bool = true, save::Bool = false, savepath::String = pwd())

    # Download and read in population sizes for Scotland
    DataRegistryUtils.load_array!(db, "human/demographics/population/scotland", "/grid area/age/persons"; sql_alias="km_age_persons_arr")
    scotpop = get_3d_km_grid_axis_array(db, ["grid_area", "age_aggr"], "val", "scottish_population_view")
    scotpop = shrink_to_active(scotpop)

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

    # Prob of developing symptoms
    p_s = fill(read_estimate(db,
        "human/infection/SARS-CoV-2/%",
        "symptom-probability",
        key="value",
        data_type=Float64)[1], age_categories)

    param_tab = read_table(db, "prob_hosp_and_cfr/data_for_scotland", "cfr_byage")
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
        db,
        "human/infection/SARS-CoV-2/latent-period",
        "latent-period",
        key="value",
        data_type=Float64
    )[1] * Unitful.hr)

    # Time asymptomatic
    T_asym = days(read_estimate(
        db,
        "human/infection/SARS-CoV-2/asymptomatic-period",
        "asymptomatic-period",
        key="value",
        data_type=Float64
    )[1] * Unitful.hr)
    @show T_asym


    # Time pre-symptomatic
    T_presym = 1.5days
    # Time symptomatic
    T_sym = days(read_estimate(
        db,
        "human/infection/SARS-CoV-2/infectious-duration",
        "infectious-duration",
        key="value",
        data_type=Float64
    )[1] * Unitful.hr) - T_presym
    # Time in hospital
    T_hosp = read_estimate(
        db,
        "fixed-parameters/T_hos",
        "T_hos",
        key="value",
        data_type=Float64
    )[1] * days
    # Time to recovery if symptomatic
    T_rec = read_estimate(
        db,
        "fixed-parameters/T_rec",
        "T_rec",
        key="value",
        data_type=Float64
    )[1] * days

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

    cat_idx = reshape(1:(nrow(abun_h) * age_categories), age_categories, nrow(abun_h))
    transitiondat = DataFrame([
        (from="Susceptible", from_id=cat_idx[:, 1], to="Exposed", to_id=cat_idx[:, 2], prob=(env = beta_env, force = beta_force)),
        (from="Exposed", from_id=cat_idx[:, 2], to="Asymptomatic", to_id=cat_idx[:, 3], prob=mu_1),
        (from="Exposed", from_id=cat_idx[:, 2], to="Presymptomatic", to_id=cat_idx[:, 4], prob=mu_2),
        (from="Presymptomatic", from_id=cat_idx[:, 4], to="Symptomatic", to_id=cat_idx[:, 5], prob=mu_3),
        (from="Symptomatic", from_id=cat_idx[:, 5], to="Hospitalised", to_id=cat_idx[:, 6], prob=hospitalisation),
        (from="Asymptomatic", from_id=cat_idx[:, 3], to="Recovered", to_id=cat_idx[:, 7], prob=sigma_1),
        (from="Symptomatic", from_id=cat_idx[:, 5], to="Recovered", to_id=cat_idx[:, 7], prob=sigma_2),
        (from="Hospitalised", from_id=cat_idx[:, 6], to="Recovered", to_id=cat_idx[:, 7], prob=sigma_hospital),
        (from="Symptomatic", from_id=cat_idx[:, 5], to="Dead", to_id=cat_idx[:, 8], prob=death_home),
        (from="Hospitalised", from_id=cat_idx[:, 6], to="Dead", to_id=cat_idx[:, 8], prob=death_hospital)
    ])

    total_pop = dropdims(sum(Float64.(scotpop), dims=3), dims=3)
    total_pop[total_pop .≈ 0.0] .= NaN
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
    epilist = SpeciesList(traits, abun_v, abun_h, movement, transitiondat, param, age_categories, movement_balance)
    rel = Gauss{eltype(epienv.habitat)}()

    transitions = create_transition_list()
    addtransition!(transitions, UpdateEpiEnvironment(update_epi_environment!))
    addtransition!(transitions, SeedInfection(seedinfected!))

    initial_infecteds = 100
    # Create epi system with all information
    @time epi = Ecosystem(epilist, epienv, rel, scotpop, UInt32(1),
    initial_infected = initial_infecteds, transitions = transitions)

    for loc in epi.cache.ordered_active
        addtransition!(epi.transitions, ViralLoad(loc, param.virus_decay))
        for age in 1:age_categories
            # Add Force
            addtransition!(epi.transitions, ForceProduce(cat_idx[age, 3], loc, param.virus_growth[age, 1]))
            addtransition!(epi.transitions, ForceProduce(cat_idx[age, 4], loc, param.virus_growth[age, 2]))
            addtransition!(epi.transitions, ForceProduce(cat_idx[age, 5], loc, param.virus_growth[age, 3]))
            addtransition!(epi.transitions, ForceDisperse(cat_idx[age, 3], loc))
            addtransition!(epi.transitions, ForceDisperse(cat_idx[age, 4], loc))
            addtransition!(epi.transitions, ForceDisperse(cat_idx[age, 5], loc))
            # Exposure
            addtransition!(epi.transitions, Exposure(transitiondat[1, :from_id][age], loc,
                transitiondat[1, :to_id][age], transitiondat[1, :prob].force[age], transitiondat[1, :prob].env[age]))
            # Infected but asymptomatic
            addtransition!(epi.transitions, Infection(transitiondat[2, :from_id][age], loc,
                transitiondat[2, :to_id][age], transitiondat[2, :prob][age]))
            # Infected but presymptomatic
            addtransition!(epi.transitions, Infection(transitiondat[3, :from_id][age], loc,
                transitiondat[3, :to_id][age], transitiondat[3, :prob][age]))
            # Develop symptoms
            addtransition!(epi.transitions, DevelopSymptoms(transitiondat[4, :from_id][age], loc,
                transitiondat[4, :to_id][age], transitiondat[4, :prob][age]))
            # Hospitalise
            addtransition!(epi.transitions, Hospitalise(transitiondat[5, :from_id][age], loc,
                transitiondat[5, :to_id][age], transitiondat[5, :prob][age]))
            # Recover without symptoms
            addtransition!(epi.transitions, Recovery(transitiondat[6, :from_id][age], loc,
                transitiondat[6, :to_id][age], transitiondat[6, :prob][age]))
            # Recover with symptoms
            addtransition!(epi.transitions, Recovery(transitiondat[7, :from_id][age], loc,
                transitiondat[7, :to_id][age], transitiondat[7, :prob][age]))
            # Recover from hospital
            addtransition!(epi.transitions, Recovery(transitiondat[8, :from_id][age], loc,
                transitiondat[8, :to_id][age], transitiondat[8, :prob][age]))
            # Die from home
            addtransition!(epi.transitions, DeathFromInfection(transitiondat[9, :from_id][age], loc,
                transitiondat[9, :to_id][age], transitiondat[9, :prob][age]))
            # Die from hospital
            addtransition!(epi.transitions, DeathFromInfection(transitiondat[10, :from_id][age], loc,
                transitiondat[10, :to_id][age], transitiondat[10, :prob][age]))
        end
    end

    N_cells = size(epi.abundances.matrix, 2)

    # Turn off work moves for <20s and >70s
    epi.spplist.species.home_balance[cat_idx[1:2, :]] .= 1.0
    epi.spplist.species.home_balance[cat_idx[7:10, :]] .= 1.0
    epi.spplist.species.work_balance[cat_idx[1:2, :]] .= 0.0
    epi.spplist.species.work_balance[cat_idx[7:10, :]] .= 0.0

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

data_dir= "Epidemiology/data/"
config = "Epidemiology/data_config.yaml"
view_sql = "Epidemiology/Scotland_run_view.sql"
db = initialise_local_registry(data_dir, data_config = config, sql_file = view_sql)

times = 1month; interval = 1day; timestep = 1day
run_model(db, times, interval, timestep, do_plot = true)

#
# pollution = parse_pollution(api)
# pollution = pollution[5513m .. 470513m, 531500m .. 1221500m, "pm2-5"]
# ukclimateAE(pollution, size(total_pop), area, fill(true, total_pop), Lockdown(20days))
