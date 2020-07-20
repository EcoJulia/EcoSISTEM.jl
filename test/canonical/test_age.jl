using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Test
using DataFrames

# sort out settings to potentially save inputs/outputs of `simulate`
do_save = (@isdefined do_save) ? do_save : false
save_path = (@isdefined save_path) ? save_path : pwd()

age_cats = [1, 2, 4]
numclasses = 3
abuns = Vector{Array{Int64, 3}}(undef, length(age_cats))
sumabuns = Vector{Array{Int64, 2}}(undef, length(age_cats))
for i in eachindex(age_cats)
    numvirus = age_cats[i] + 1
    # Set simulation parameters
    birth = fill(0.0/day, numclasses, age_cats[i])
    death = fill(0.0/day, numclasses, age_cats[i])
    age_mixing = fill(1.0, age_cats[i], age_cats[i])
    beta_force = fill(1.0/day, age_cats[i])
    beta_env = fill(1.0/day, age_cats[i])
    sigma = fill(0.02/day, age_cats[i])
    virus_growth = fill(1e-2/day, age_cats[i])
    virus_decay = 1.0/2day
    param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force, age_mixing = age_mixing)
    paramDat = DataFrame([["Infected"], ["Susceptible"], [sigma]], [:from, :to, :prob])

    # Set up simple gridded environment
    grid = (4, 4)
    area = 525_000.0km^2
    epienv = simplehabitatAE(298.0K, grid, area, NoControl())

    # Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
    virus_env = 0
    susceptible = fill(Int64.(50_000_000/age_cats[i]), age_cats[i])
    infected = fill(Int64.(10_000/age_cats[i]), age_cats[i])
    dead = fill(0, age_cats[i])
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
    abun_v = (Environment = virus_env, Force = fill(0, age_cats[i]))

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(100.0km, numclasses * age_cats[i])
    cat_idx = reshape(1:(numclasses * age_cats[i]), age_cats[i], numclasses)
    dispersal_dists[cat_idx[:, 2]] .= 1_000.0km
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = AlwaysMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, paramDat, param, age_cats[i])

    # Create epi system with all information
    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel)

    # Run simulation
    times = 2years; interval = 1day; timestep = 1day
    abuns[i] = zeros(Int64, numclasses * age_cats[i], 16, convert(Int64, floor(times / interval)) + 1)
    thisabun = abuns[i]
    @time simulate_record!(thisabun, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "age_cats_$(age_cats[i])"))

    # Test no-one dies (death rate = 0)
    @test sum(thisabun[end, :, :]) == 0
    # Test overall population size stays constant (birth rate = death rate = 0)
    @test all(sum(thisabun[1:(2 * age_cats[i]), :, :], dims = (1, 2)) .== sum(susceptible + infected))
    sumabuns[i] = sum(abuns[i], dims = 2)[:, 1, :]
end

# For each disease category, check trajectory is the same when we change grid size

for j in 2:length(sumabuns)
    for i in 1:numclasses
        cat_idx1 = reshape(1:(numclasses * age_cats[j - 1]), age_cats[j - 1], numclasses)
        cat_idx2 = reshape(1:(numclasses * age_cats[j]), age_cats[j], numclasses)
        @test isapprox(sum(sumabuns[j-1][cat_idx1[:, i], :], dims = 1), sum(sumabuns[j][cat_idx2[:, i], :], dims = 1), rtol = 5e-2)
    end
end
