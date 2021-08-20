using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Test
using DataFrames
using LinearAlgebra
using StatsBase

@testset "mixing" begin

    # sort out settings to potentially save inputs/outputs of `simulate`
    do_save = (@isdefined do_save) ? do_save : false
    save_path = (@isdefined save_path) ? save_path : pwd()

    mixing_rate = [1.0, 0.5, 0.1]
    age_cats = 4
    abuns = Vector{Array{Int64, 3}}(undef, length(mixing_rate))
    sumabuns = Vector{Array{Int64, 2}}(undef, length(mixing_rate))
    for i in eachindex(mixing_rate)
        # Set up simple gridded environment
        grid = (4, 4)
        area = 525_000.0km^2
        epienv = simplehabitatAE(298.0K, grid, area, NoControl())

        # Set initial population sizes for all pathogen categories
        virus_env = 0
        virus_force = fill(0, age_cats)
        abun_v = DataFrame([
            (name="Environment", initial=virus_env),
            (name="Force", initial=virus_force),
        ])
        numvirus = sum(length.(abun_v.initial))

        # Set initial population sizes for all human categories
        susceptible = fill(Int64.(50_000_000/age_cats), age_cats)
        infected = fill(Int64.(10_000/age_cats), age_cats)
        dead = fill(0, age_cats)
        abun_h = DataFrame([
            (name="Susceptible", type=Susceptible, initial=susceptible),
            (name="Infected", type=Infectious, initial=infected),
            (name="Dead", type=Removed, initial=dead),
        ])
        numclasses = nrow(abun_h)
        numstates = sum(length.(abun_h.initial))

        # Set non-pathogen mediated transitions
        sigma = fill(0.02/day, age_cats)
        transitions = DataFrame([
            (from="Infected", to="Susceptible", prob=sigma),
        ])

        # Set simulation parameters
        birth = fill(0.0/day, numclasses, age_cats)
        death = fill(0.0/day, numclasses, age_cats)
        age_mixing = fill(mixing_rate[i], age_cats, age_cats)
        age_mixing[diagind(age_mixing)] .= 1.0 # 1 so we mix with ourselves
        beta_force = fill(1.0/day, age_cats)
        beta_env = fill(1.0/day, age_cats)
        virus_growth = fill(1e-2/day, age_cats)
        virus_decay = 1.0/2day
        param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force, age_mixing = age_mixing)

        # Dispersal kernels for virus and disease classes
        dispersal_dists = fill(1_000.0km, prod(grid))
        kernel = GaussianKernel.(dispersal_dists, 1e-10)
        movement = EpiMovement(kernel)

        # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
        traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
        epilist = SpeciesList(traits, abun_v, abun_h, movement, transitions, param, age_cats)

        # Create epi system with all information
        rel = Gauss{eltype(epienv.habitat)}()
        epi = Ecosystem(epilist, epienv, rel)

        # Run simulation
        times = 2years; interval = 1day; timestep = 1day
        abuns[i] = zeros(Int64, numstates, prod(grid), div(times, interval) + 1)
        thisabun = abuns[i]
        @time simulate_record!(thisabun, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "mixing_rate_$(mixing_rate[i])"))

        # Find correct indices in arrays
        row_sus = findfirst(==("Susceptible"), abun_h.name)
        row_inf = findfirst(==("Infected"), abun_h.name)
        row_dead = findfirst(==("Dead"), abun_h.name)

        true_indices = [0; cumsum(length.(abun_h.initial))]
        idx_sus = (true_indices[row_sus]+1):true_indices[row_sus+1]
        idx_inf = (true_indices[row_inf]+1):true_indices[row_inf+1]
        idx_dead = (true_indices[row_dead]+1):true_indices[row_dead+1]

        # Test no-one dies (death rate = 0)
        @test sum(thisabun[idx_dead, :, :]) == 0
        # Test overall population size stays constant (birth rate = death rate = 0)
        @test all(sum(thisabun[idx_sus, :, :] .+ thisabun[idx_inf, :, :],
                      dims = 2)[:,1,:] .==
                  abun_h.initial[row_sus] .+ abun_h.initial[row_inf])
        sumabuns[i] = sum(abuns[i], dims = 2)[:, 1, :]
    end

    # Check that more are infected when mixing rates are higher
    numclasses = Int(size(abuns[1], 1)/age_cats)
    cat_idx = reshape(1:(numclasses * age_cats), age_cats, numclasses)
    for j in 2:length(sumabuns)
        @test mean(sumabuns[j-1][cat_idx[:, 2], :] .>= sumabuns[j][cat_idx[:, 2], :]) >= 0.95
    end

end


@testset "mixing transitions" begin

    # sort out settings to potentially save inputs/outputs of `simulate`
    do_save = (@isdefined do_save) ? do_save : false
    save_path = (@isdefined save_path) ? save_path : pwd()

    mixing_rate = [1.0, 0.5, 0.1]
    age_cats = 4
    abuns = Vector{Array{Int64, 3}}(undef, length(mixing_rate))
    sumabuns = Vector{Array{Int64, 2}}(undef, length(mixing_rate))
    for i in eachindex(mixing_rate)
        # Set up simple gridded environment
        grid = (4, 4)
        area = 525_000.0km^2
        epienv = simplehabitatAE(298.0K, grid, area, NoControl())

        # Set initial population sizes for all pathogen categories
        virus_env = 0
        virus_force = fill(0, age_cats)
        abun_v = DataFrame([
            (name="Environment", initial=virus_env),
            (name="Force", initial=virus_force),
        ])
        numvirus = sum(length.(abun_v.initial))

        # Set initial population sizes for all human categories
        susceptible = fill(Int64.(50_000_000/age_cats), age_cats)
        infected = fill(Int64.(10_000/age_cats), age_cats)
        dead = fill(0, age_cats)
        abun_h = DataFrame([
            (name="Susceptible", type=Susceptible, initial=susceptible),
            (name="Infected", type=Infectious, initial=infected),
            (name="Dead", type=Removed, initial=dead),
        ])
        numclasses = nrow(abun_h)
        numstates = sum(length.(abun_h.initial))
        cat_idx = reshape(1:(nrow(abun_h) * age_cats), age_cats, nrow(abun_h))

        # Set simulation parameters
        birth = fill(0.0/day, numclasses, age_cats)
        death = fill(0.0/day, numclasses, age_cats)
        age_mixing = fill(mixing_rate[i], age_cats, age_cats)
        age_mixing[diagind(age_mixing)] .= 1.0 # 1 so we mix with ourselves
        beta_force = fill(1.0/day, age_cats)
        beta_env = fill(1.0/day, age_cats)
        virus_growth = fill(1e-2/day, age_cats)
        virus_decay = 1.0/2day
        param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force, age_mixing = age_mixing)

        # Set non-pathogen mediated transitions
        sigma = fill(0.02/day, age_cats)
        transitiondat = DataFrame([
            (from="Susceptible", from_id=cat_idx[:,1], to="Infected", to_id = cat_idx[:,2], prob=(env=beta_env, force = beta_force)),
            (from="Infected", from_id = cat_idx[:,2], to="Susceptible", to_id = cat_idx[:,1], prob=sigma),
        ])

        # Dispersal kernels for virus and disease classes
        dispersal_dists = fill(1_000.0km, prod(grid))
        kernel = GaussianKernel.(dispersal_dists, 1e-10)
        movement = EpiMovement(kernel)

        # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
        traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
        epilist = SpeciesList(traits, abun_v, abun_h, movement, transitiondat, param, age_cats)

        transitions = create_transition_list()
        addtransition!(transitions, SeedInfection(seedinfected!))
        addtransition!(transitions, UpdateEpiEnvironment(update_epi_environment!))
        
        for loc in eachindex(epienv.habitat.matrix)
            addtransition!(transitions, ViralLoad(loc, param.virus_decay))
            for age in 1:age_cats
                addtransition!(transitions, ForceProduce(cat_idx[age, 2], loc, param.virus_growth[age]))
                addtransition!(transitions, ForceDisperse(cat_idx[age, 2], loc))
                addtransition!(transitions, Exposure(transitiondat[1, :from_id][age], loc,
                    transitiondat[1, :to_id][age], transitiondat[1, :prob].force[age], transitiondat[1, :prob].env[age]))
                addtransition!(transitions, Recovery(transitiondat[2, :from_id][age], loc,
                    transitiondat[2, :to_id][age], transitiondat[2, :prob][age]))
            end
        end

        # Create epi system with all information
        rel = Gauss{eltype(epienv.habitat)}()
        epi = Ecosystem(epilist, epienv, rel, transitions = transitions)

        # Run simulation
        times = 2years; interval = 1day; timestep = 1day
        abuns[i] = zeros(Int64, numstates, prod(grid), div(times, interval) + 1)
        thisabun = abuns[i]
        @time simulate_record!(thisabun, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "mixing_rate_$(mixing_rate[i])"))

        # Find correct indices in arrays
        row_sus = findfirst(==("Susceptible"), abun_h.name)
        row_inf = findfirst(==("Infected"), abun_h.name)
        row_dead = findfirst(==("Dead"), abun_h.name)

        true_indices = [0; cumsum(length.(abun_h.initial))]
        idx_sus = (true_indices[row_sus]+1):true_indices[row_sus+1]
        idx_inf = (true_indices[row_inf]+1):true_indices[row_inf+1]
        idx_dead = (true_indices[row_dead]+1):true_indices[row_dead+1]

        # Test no-one dies (death rate = 0)
        @test sum(thisabun[idx_dead, :, :]) == 0
        # Test overall population size stays constant (birth rate = death rate = 0)
        @test all(sum(thisabun[idx_sus, :, :] .+ thisabun[idx_inf, :, :],
                      dims = 2)[:,1,:] .==
                  abun_h.initial[row_sus] .+ abun_h.initial[row_inf])
        sumabuns[i] = sum(abuns[i], dims = 2)[:, 1, :]
    end

    # Check that more are infected when mixing rates are higher
    numclasses = Int(size(abuns[1], 1)/age_cats)
    cat_idx = reshape(1:(numclasses * age_cats), age_cats, numclasses)
    for j in 2:length(sumabuns)
        @test mean(sumabuns[j-1][cat_idx[:, 2], :] .>= sumabuns[j][cat_idx[:, 2], :]) >= 0.95
    end

end
