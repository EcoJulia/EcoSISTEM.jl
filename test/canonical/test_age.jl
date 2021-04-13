using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Test
using DataFrames

@testset "age" begin

# sort out settings to potentially save inputs/outputs of `simulate`
do_save = (@isdefined do_save) ? do_save : false
save_path = (@isdefined save_path) ? save_path : pwd()

age_cats = [1, 2, 4]
abuns = Vector{Array{Int64, 3}}(undef, length(age_cats))
sumabuns = Vector{Array{Int64, 2}}(undef, length(age_cats))
for i in eachindex(age_cats)
    # Set up simple gridded environment
    grid = (4, 4)
    area = 525_000.0km^2
    epienv = simplehabitatAE(298.0K, grid, area, NoControl())

    # Set initial population sizes for all pathogen categories
    virus_env = 0
    virus_force = fill(0, age_cats[i])
    abun_v = DataFrame([
        (name="Environment", initial=virus_env),
        (name="Force", initial=virus_force),
    ])
    numvirus = sum(length.(abun_v.initial))

    # Set initial population sizes for all human categories
    susceptible = fill(Int64(50_000_000/age_cats[i]), age_cats[i])
    infected = fill(Int64(10_000/age_cats[i]), age_cats[i])
    dead = fill(0, age_cats[i])
    abun_h = DataFrame([
        (name="Susceptible", type=Susceptible, initial=susceptible),
        (name="Infected", type=Infectious, initial=infected),
        (name="Dead", type=Removed, initial=dead),
    ])
    numclasses = nrow(abun_h)
    numstates = sum(length.(abun_h.initial))

    # Set non-pathogen mediated transitions
    sigma = fill(0.02/day, age_cats[i])
    transitions = DataFrame([
        (from="Infected", to="Susceptible", prob=sigma),
    ])

    # Set simulation parameters
    birth = fill(0.0/day, numclasses, age_cats[i])
    death = fill(0.0/day, numclasses, age_cats[i])
    age_mixing = fill(1.0, age_cats[i], age_cats[i])
    beta_force = fill(10.0/day, age_cats[i])
    beta_env = fill(10.0/day, age_cats[i])
    virus_growth = fill(1e-2/day, age_cats[i])
    virus_decay = 1.0/2day
    param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force, age_mixing = age_mixing)

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(1_000.0km, prod(grid))
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param, age_cats[i])

    # Create epi system with all information
    rel = Gauss{eltype(epienv.habitat)}()
    epi = Ecosystem(epilist, epienv, rel)

    # Run simulation
    times = 2years; interval = 1day; timestep = 1day
    abuns[i] = zeros(Int64, numstates, prod(grid), div(times, interval) + 1)
    thisabun = abuns[i]
    @time epi_simulate_record!(thisabun, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "age_cats_$(age_cats[i])"))

    # Test no-one dies (death rate = 0)
    @test sum(thisabun[end, :, :]) == 0
    # Test overall population size stays constant (birth rate = death rate = 0)
    @test all(sum(thisabun[1:(2 * age_cats[i]), :, :], dims = (1, 2)) .== sum(susceptible + infected))
    sumabuns[i] = sum(abuns[i], dims = 2)[:, 1, :]
end

# For each disease category, check trajectory is the same when we change grid size

for j in 2:length(sumabuns)
    numclasses = 3
    for i in 1:numclasses
        cat_idx1 = reshape(1:(numclasses * age_cats[j - 1]),
                           age_cats[j - 1], numclasses)
        cat_idx2 = reshape(1:(numclasses * age_cats[j]),
                           age_cats[j], numclasses)
        @test isapprox(sum(sumabuns[j-1][cat_idx1[:, i], :], dims = 1),
                       sum(sumabuns[j][cat_idx2[:, i], :], dims = 1),
                       rtol = 5e-2)
    end
end

end
