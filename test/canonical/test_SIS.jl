using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Test

# sort out settings to potentially save inputs/outputs of `simulate`
do_save = (@isdefined do_save) ? do_save : false
save_path = (@isdefined save_path) ? save_path : pwd()

grid_sizes = [4, 8, 16]
numclasses = 3
numvirus = 2
abuns = Vector{Array{Int64, 3}}(undef, length(grid_sizes))
sumabuns = Vector{Array{Int64, 2}}(undef, length(grid_sizes))
for i in eachindex(grid_sizes)
    # Set simulation parameters
    birth = fill(0.0/day, numclasses)
    death = fill(0.0/day, numclasses)
    beta_force = 1.0/day
    beta_env = 1.0/day
    sigma = 0.02/day
    virus_growth = 1e-2/day
    virus_decay = 1.0/2day
    param = SISGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
    param = transition(param)

    # Set up simple gridded environment
    grid = (grid_sizes[i], grid_sizes[i])
    area = 525_000.0km^2
    epienv = simplehabitatAE(298.0K, grid, area, NoControl())

    # Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
    virus = 0
    susceptible = 1_000_000 * maximum(grid_sizes)^2
    infected = 250 * maximum(grid_sizes)^2
    dead = 0
    abun_h = (
      Susceptible = susceptible,
      Infected = infected,
      Dead = dead
    )
    disease_classes = (
        susceptible = ["Susceptible"],
        infectious = ["Infected"]
    )
    abun_v = (Environment = virus, Force = 0)

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(100.0km, numclasses)
    dispersal_dists[2] = 1_000.0km
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = AlwaysMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, param)

    # Create epi system with all information
    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel)

    # Run simulation
    times = 2years; interval = 1day; timestep = 1day
    abuns[i] = zeros(Int64, size(epi.abundances.matrix, 1), grid_sizes[i]^2,
                     convert(Int64, floor(times / interval)) + 1)
    thisabun = abuns[i]
    @time simulate_record!(thisabun, epi, times, interval, timestep; save=do_save, save_path=joinpath(save_path, "grid_size_$(grid_sizes[i])"))

    # Test no-one dies (death rate = 0)
    @test sum(thisabun[end, :, :]) == 0
    # Test overall population size stays constant (birth rate = death rate = 0)
    @test all(sum(thisabun[1:2, :, :], dims = (1, 2)) .== (susceptible + infected))
    sumabuns[i] = sum(abuns[i], dims = 2)[:, 1, :]
end

# For each disease category, check trajectory is the same when we change grid size
for j in 2:length(sumabuns)
    for i in 1:numclasses
        @test isapprox(sumabuns[j-1][i, :], sumabuns[j][i, :], rtol = 5e-2)
    end
end
