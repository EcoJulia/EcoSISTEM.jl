using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Test

grid_sizes = [1, 2, 4]
numclasses = 5
abuns = zeros(Int64, 3, numclasses, 16, 731)
for i in eachindex(grid_sizes)
    # Set simulation parameters
    birth = fill(0.0/day, numclasses)
    death = fill(0.0/day, numclasses)
    beta = 1e-4/day
    sigma = 0.02/day
    virus_growth = 1e-3/day
    virus_decay = 1.0/day
    param = SIRGrowth{typeof(unit(beta))}(birth, death, virus_growth, virus_decay, beta, sigma)
    param = transition(param)

    # Set up simple gridded environment
    grid = (grid_sizes[i], grid_sizes[i])
    area = 525_000.0km^2
    epienv = simplehabitatAE(298.0K, grid, area, NoControl())

    # Set initial population sizes for all categories: Virus, Susceptible, Infected, Recovered
    virus = 10_000 * grid_sizes[i]^2
    susceptible = 5_000_000
    infected = 1_000
    recovered = 0
    dead = 0
    abun = [virus, susceptible, infected, recovered, dead]

    # Dispersal kernels for virus and disease classes
    dispersal_dists = fill(100.0km, numclasses)
    dispersal_dists[3] = 700.0km
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = AlwaysMovement(kernel)

    # Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
    traits = GaussTrait(fill(298.0K, numclasses), fill(0.1K, numclasses))
    epilist = SIR(traits, abun, movement, param)

    # Create epi system with all information
    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel)

    # Run simulation
    times = 2years; interval = 1day; timestep = 1day
    thisabun = @view abuns[i, :, 1:grid_sizes[i]^2, :]
    @time simulate_record!(thisabun, epi, times, interval, timestep)
    # Test no-one dies (death rate = 0)
    @test sum(thisabun[end, :, :]) == 0
    # Test overall population size stays constant (birth rate = death rate = 0)
    @test all(sum(thisabun[2:4, :, :], dims = (1, 2)) .== (susceptible + infected + recovered))
end

abuns = sum(abuns, dims = 3)[:, :, 1, :]
# For each disease category, check trajectory is the same when we change grid size
for i in 2:numclasses
    @test isapprox(abuns[1, i, :], abuns[2, i, :], rtol = 1e-2)
    @test isapprox(abuns[2, i, :], abuns[3, i, :], rtol = 1e-2)
end
