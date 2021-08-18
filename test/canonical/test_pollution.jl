using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using StatsBase
using Test
using DataFrames
using EcoSISTEM: human, virus

import EcoSISTEM.getprob
@testset "Pollution" begin
    # Define how probabilities for transitions should be altered by pollution
    getprob(eco::Ecosystem, rule::Exposure) = (rule.force_prob * 2.0 * get_env(eco.abenv.habitat, rule.location), rule.virus_prob)
    getprob(eco::Ecosystem, rule::DevelopSymptoms) = rule.prob * 2.0 * get_env(eco.abenv.habitat, rule.location)
    getprob(eco::Ecosystem, rule::Hospitalise) = rule.prob * 2.0 * get_env(eco.abenv.habitat, rule.location)
    getprob(eco::Ecosystem, rule::DeathFromInfection) = rule.prob * 5.0 * get_env(eco.abenv.habitat, rule.location)


    # sort out settings to potentially save inputs/outputs of `simulate`
    do_save = (@isdefined do_save) ? do_save : false
    save_path = (@isdefined save_path) ? save_path : pwd()

    # Set up simple gridded environment
    grid = (4, 4)
    area = 525_000.0km^2
    epienv = simplehabitatAE(20.0, grid, area, NoControl())

    # Set initial population sizes for all pathogen categories
    abun_v = DataFrame([
        (name="Environment", initial=0),
        (name="Force", initial=0),
    ])
    numvirus = nrow(abun_v)

    # Set initial population sizes for all human categories
    susceptible = 5_000_000
    abun_h = DataFrame([
    (name="Susceptible", type=Susceptible, initial=susceptible),
    (name="Exposed", type=OtherDiseaseState, initial=0),
    (name="Asymptomatic", type=Infectious, initial=0),
    (name="Presymptomatic", type=Infectious, initial=0),
    (name="Symptomatic", type=Infectious, initial=0),
    (name="Hospitalised", type=OtherDiseaseState, initial=0),
    (name="Recovered", type=Removed, initial=0),
    (name="Dead", type=Removed, initial=0),
    ])
    numclasses = nrow(abun_h)

    # Prob of developing symptoms
    p_s = 1.0
    # Prob of hospitalisation
    p_h = 0.2
    # Case fatality ratio
    cfr_home = cfr_hosp = 1.0
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
    mu_1 = (1 - p_s) * 1/T_lat
    # Exposed -> Pre-symptomatic
    mu_2 = p_s * 1/T_lat
    # Pre-symptomatic -> symptomatic
    mu_3 = 1/T_presym
    # Symptomatic -> hospital
    hospitalisation = p_h * 1/T_sym
    # Asymptomatic -> recovered
    sigma_1 = 1/T_asym
    # Symptomatic -> recovered
    sigma_2 = (1 - p_h) * (1 - cfr_home) * 1/T_rec
    # Hospital -> recovered
    sigma_hospital = (1 - cfr_hosp) * 1/T_hosp
    # Symptomatic -> death
    death_home = cfr_home * 2/T_hosp
    # Hospital -> death
    death_hospital = cfr_hosp * 1/T_hosp

    ## High transmission & 100% case fatality
    birth = fill(0.0/day, numclasses)
    death = fill(0.0/day, numclasses)
    virus_growth_asymp = virus_growth_presymp = virus_growth_symp = 1.0/day
    virus_decay = 1/3days
    beta_force = 1e3/day
    beta_env = 1e3/day

    transitiondat = DataFrame([
        (from="Susceptible", from_id=1, to="Exposed", to_id=2, prob=(env = beta_env, force = beta_force)),
        (from="Exposed", from_id=2, to="Asymptomatic", to_id=3, prob=mu_1),
        (from="Exposed", from_id=2, to="Presymptomatic", to_id=4, prob=mu_2),
        (from="Presymptomatic", from_id=4, to="Symptomatic", to_id=5, prob=mu_3),
        (from="Symptomatic", from_id=5, to="Hospitalised", to_id=6, prob=hospitalisation),
        (from="Asymptomatic", from_id=3, to="Recovered", to_id=7, prob=sigma_1),
        (from="Symptomatic", from_id=5, to="Recovered", to_id=7, prob=sigma_2),
        (from="Hospitalised", from_id=6, to="Recovered", to_id=7, prob=sigma_hospital),
        (from="Symptomatic", from_id=5, to="Dead", to_id=8, prob=death_home),
        (from="Hospitalised", from_id=6, to="Dead", to_id=8, prob=death_hospital)
    ])

    param = (birth = birth, death = death, virus_growth = [virus_growth_asymp virus_growth_presymp virus_growth_symp], virus_decay = virus_decay, beta_force = beta_force, beta_env = beta_env)

    # Dispersal kernels for virus dispersal from different disease classes
    dispersal_dists = fill(500.0km, prod(grid))

    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(20.0μg/m^3, numvirus), fill(0.1μg/m^3, numvirus))
    epilist = SpeciesList(traits, abun_v, abun_h, movement, transitiondat, param)
    rel = Gauss{eltype(epienv.habitat)}()

    # Create epi system with all information
    new_exposed = 100
    new_virus = 1_000
    transitions = create_transition_list()
    addtransition!(transitions, UpdateEpiEnvironment(update_epi_environment!))
    addtransition!(transitions, SeedInfection(seedinfected!))

    epi = Ecosystem(epilist, epienv, rel, transitions = transitions)
    virus(epi.abundances)[1, 1] = new_virus
    human(epi.abundances)[2, 1] = new_exposed

    for loc in epi.cache.ordered_active
        addtransition!(epi.transitions, ViralLoad(loc, param.virus_decay))
        # Add Force
        addtransition!(epi.transitions, ForceProduce(3, loc, param.virus_growth[1]))
        addtransition!(epi.transitions, ForceProduce(4, loc, param.virus_growth[2]))
        addtransition!(epi.transitions, ForceProduce(5, loc, param.virus_growth[3]))
        addtransition!(epi.transitions, ForceDisperse(3, loc))
        addtransition!(epi.transitions, ForceDisperse(4, loc))
        addtransition!(epi.transitions, ForceDisperse(5, loc))
        # Exposure
        addtransition!(epi.transitions, Exposure(transitiondat[1, :from_id], loc, transitiondat[1, :to_id],
            transitiondat[1, :prob].force, transitiondat[1, :prob].env))
        # Infected but asymptomatic
        addtransition!(epi.transitions, Infection(transitiondat[2, :from_id], loc,
            transitiondat[2, :to_id], transitiondat[2, :prob]))
        # Infected but presymptomatic
        addtransition!(epi.transitions, Infection(transitiondat[3, :from_id], loc,
            transitiondat[3, :to_id], transitiondat[3, :prob]))
        # Develop symptoms
        addtransition!(epi.transitions, DevelopSymptoms(transitiondat[4, :from_id], loc,
            transitiondat[4, :to_id], transitiondat[4, :prob]))
        # Hospitalise
        addtransition!(epi.transitions, Hospitalise(transitiondat[5, :from_id], loc,
            transitiondat[5, :to_id], transitiondat[5, :prob]))
        # Recover without symptoms
        addtransition!(epi.transitions, Recovery(transitiondat[6, :from_id], loc,
            transitiondat[6, :to_id], transitiondat[6, :prob]))
        # Recover with symptoms
        addtransition!(epi.transitions, Recovery(transitiondat[7, :from_id], loc,
            transitiondat[7, :to_id], transitiondat[7, :prob]))
        # Recover from hospital
        addtransition!(epi.transitions, Recovery(transitiondat[8, :from_id], loc,
            transitiondat[8, :to_id], transitiondat[8, :prob]))
        # Die from home
        addtransition!(epi.transitions, DeathFromInfection(transitiondat[9, :from_id], loc,
            transitiondat[9, :to_id], transitiondat[9, :prob]))
        # Die from hospital
        addtransition!(epi.transitions, DeathFromInfection(transitiondat[10, :from_id], loc,
            transitiondat[10, :to_id], transitiondat[10, :prob]))
    end
    

    # Run simulation
    abuns = zeros(Int64, size(epi.abundances.matrix, 1), size(epi.abundances.matrix, 2), 366)
    times = 1year; interval = 1day; timestep = 1day
    @time simulate_record!(abuns, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "high_trans"))


    # Test everyone becomes infected and dies
    @test sum(abuns[1, :, 365]) == 0
    @test sum(abuns[end, :, 365]) == (susceptible + new_exposed)

end
