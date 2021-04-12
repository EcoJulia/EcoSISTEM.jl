using RCall
using EcoSISTEM
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using Unitful
using HTTP
using AxisArrays
using Unitful.DefaultSymbols
using Plots
using DataFrames
import Dates.DateTime

# Download climate data and write to HDF5
startDate = DateTime("2020-01-01")
endDate = DateTime("2020-01-30")
uktemp = ClimatePref.MetOfficeDownload(:uk_daily, :temp_mean, "test/examples/temp", startDate, endDate)
outputfolder = "test/examples/climate.h5"
writeMet(uktemp, "temp", outputfolder)
uktemp = readMet(outputfolder, "temp")

function run_model(climatearray::AxisArray, times::Unitful.Time, interval::Unitful.Time, timestep::Unitful.Time; do_plot::Bool = false, do_download::Bool = true)
    # Download and read in population sizes for Scotland
    dir = EcoSISTEM.path("test", "TEMP")
    file = joinpath(dir, "demographics.h5")
    if do_download
        !isdir(EcoSISTEM.path("test", "TEMP")) && mkdir(EcoSISTEM.path("test", "TEMP"))
        io = open(EcoSISTEM.path("test", "TEMP", "demographics.h5"), "w")
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

    # Set initial population sizes for all pathogen categories
    abun_v = DataFrame([
        (name="Environment", initial=0),
        (name="Force", initial=fill(0, age_categories)),
    ])
    numvirus = nrow(abun_v)

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

    # Sum up age categories and turn into simple matrix
    total_pop = dropdims(sum(Float64.(scotpop), dims=3), dims=3)
    total_pop = AxisArray(total_pop, AxisArrays.axes(scotpop)[1], AxisArrays.axes(scotpop)[2])
    total_pop.data[total_pop .≈ 0.0] .= NaN

    # Prob of developing symptoms
    p_s = fill(0.96, age_categories)
    # Prob of hospitalisation
    p_h = fill(0.2, age_categories)
    # Case fatality ratio
    cfr_home = cfr_hospital = fill(0.1, age_categories)
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
    ageing = fill(0.0/day, age_categories - 1) # no ageing for now
    age_mixing = fill(1.0, age_categories, age_categories)

    param = (birth = birth_rates, death = death_rates, virus_growth = [virus_growth_asymp virus_growth_presymp virus_growth_symp], virus_decay = virus_decay, beta_force = beta_force, beta_env = beta_env, age_mixing = age_mixing)

    total_pop.data[isnan.(climatearray[:, :, 1])] .= NaN
    # ! Bug - this final argument is wrong, or a constructor is missing
    epienv = ukclimateAE(climatearray, area, NoControl(), total_pop)

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(20.0km, length(total_pop))
    cat_idx = reshape(1:(numclasses * age_categories), age_categories, numclasses)
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = AlwaysMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(279.0K, numvirus), fill(5.0K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h,
                      movement, transitions, param, age_categories)
    rel = Gauss{eltype(epienv.habitat)}()

    # Create epi system with all information
    epi = EpiSystem(epilist, epienv, rel)

    # Populate susceptibles according to actual population spread
    reshaped_pop =
        reshape(scotpop[1:size(epienv.active, 1), 1:size(epienv.active, 2), :],
                size(epienv.active, 1) * size(epienv.active, 2), size(scotpop, 3))'
    epi.abundances.matrix[cat_idx[:, 1], :] = reshaped_pop

    # Add in initial infections randomly (samples weighted by population size)
    # Define generator for all pair age x cell
    age_and_cells = Iterators.product(1:age_categories,
                                      1:size(epi.abundances.matrix, 2))
    # Take all susceptibles of each age per cell
    pop_weights = epi.abundances.matrix[vcat(cat_idx[:, 1]...), :]
    # It would be nice if it wasn't necessary to call collect here
    N_cells = size(epi.abundances.matrix, 2)
    samp = sample(collect(age_and_cells), weights(1.0 .* vec(pop_weights)), 100)
    age_ids = getfield.(samp, 1)
    cell_ids = getfield.(samp, 2)

    for i in eachindex(age_ids)
        if (epi.abundances.matrix[cat_idx[age_ids[i], 1], cell_ids[i]] > 0)
            # Add to exposed
            epi.abundances.matrix[cat_idx[age_ids[i], 2], cell_ids[i]] += 1
            # Remove from susceptible
            epi.abundances.matrix[cat_idx[age_ids[i], 1], cell_ids[i]] -= 1
        end
    end

    # Run simulation
    abuns = zeros(Int64, size(epi.abundances.matrix, 1), N_cells,
                  floor(Int, times/timestep) + 1)
    @time simulate_record!(abuns, epi, times, interval, timestep)

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
        display(plot_epiheatmaps(epi, abuns, steps = [21]))
    end
    return abuns
end

times = 2days; interval = 1day; timestep = 1day
abuns = run_model(uktemp, times, interval, timestep)
