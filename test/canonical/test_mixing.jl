using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Test
using DataFrames

@testset "mixing" begin
    # sort out settings to potentially save inputs/outputs of `simulate`
    do_save = (@isdefined do_save) ? do_save : false
    save_path = (@isdefined save_path) ? save_path : pwd()

    mixing_rate = [1.0, 0.5, 0.1]
    numclasses = 3
    age_cats = 4
    abuns = Vector{Array{Int64, 3}}(undef, length(mixing_rate))
    sumabuns = Vector{Array{Int64, 2}}(undef, length(mixing_rate))
    for i in eachindex(mixing_rate)
        numvirus = age_cats + 1
        # Set simulation parameters
        birth = fill(0.0/day, numclasses, age_cats)
        death = fill(0.0/day, numclasses, age_cats)
        age_mixing = fill(mixing_rate[i], age_cats, age_cats)
        beta_force = fill(1.0/day, age_cats)
        beta_env = fill(1.0/day, age_cats)
        sigma = fill(0.02/day, age_cats)
        virus_growth = fill(1e-2/day, age_cats)
        virus_decay = 1.0/2day
        param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force, age_mixing = age_mixing)
        paramDat = DataFrame([["Infected"], ["Susceptible"], [sigma]], [:from, :to, :prob])

        # Set up simple gridded environment
        grid = (4, 4)
        area = 525_000.0km^2
        epienv = simplehabitatAE(298.0K, grid, area, NoControl())

        # Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
        virus_env = 0
        susceptible = fill(Int64.(50_000_000/age_cats), age_cats)
        infected = fill(Int64.(10_000/age_cats), age_cats)
        dead = fill(0, age_cats)
        sus = ["Susceptible"]
        inf = ["Infected"]
        abun_h = (
          Susceptible = susceptible,
          Infected = infected,
          Dead = dead
        )
        disease_classes = (
            susceptible = ["Susceptible"],
            infectious = ["Infected"]
        )
        abun_v = (Environment = virus_env, Force = fill(0, age_cats))

        # Dispersal kernels for virus and disease classes
        dispersal_dists = fill(1_000.0km, prod(grid))
        kernel = GaussianKernel.(dispersal_dists, 1e-10)
        movement = EpiMovement(kernel)

        # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
        traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
        epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, paramDat, param, age_cats)

        # Create epi system with all information
        rel = Gauss{eltype(epienv.habitat)}()
        epi = EpiSystem(epilist, epienv, rel)

        # Run simulation
        times = 2years; interval = 1day; timestep = 1day
        abuns[i] = zeros(Int64, numclasses * age_cats, 16, convert(Int64, floor(times / interval)) + 1)
        thisabun = abuns[i]
        @time simulate_record!(thisabun, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "mixing_rate_$(mixing_rate[i])"))

        # Test no-one dies (death rate = 0)
        @test sum(thisabun[end, :, :]) == 0
        # Test overall population size stays constant (birth rate = death rate = 0)
        @test all(sum(thisabun[1:(2 * age_cats), :, :], dims = (1, 2)) .== sum(susceptible + infected))
        sumabuns[i] = sum(abuns[i], dims = 2)[:, 1, :]
    end

    # Check that more are infected when mixing rates are higher
    cat_idx = reshape(1:(numclasses * age_cats), age_cats, numclasses)
    for j in 2:length(sumabuns)
        @test sum(sum(sumabuns[j-1][cat_idx[:, 2], :], dims = 1) .>= sum(sumabuns[j][cat_idx[:, 2], :], dims = 1)) >= (0.95 * size(sumabuns[1], 2))
    end
end
