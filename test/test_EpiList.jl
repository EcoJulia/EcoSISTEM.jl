using Simulation
using Test
using Unitful.DefaultSymbols
using Distributions
using Simulation.Units
@testset "EpiList" begin
    # Set initial population sizes for all pathogen categories
    abun_v = DataFrame([
        (name="Environment", initial=10),
        (name="Force", initial=10),
    ])
    numvirus = nrow(abun_v)

    # Set initial population sizes for all human categories
    susceptible = 1_000
    infected = 1
    abun_h = DataFrame([
        (name="Susceptible", type=Susceptible, initial=susceptible),
        (name="Infected", type=Infectious, initial=infected),
        (name="Recovered", type=Removed, initial=0),
        (name="Dead", type=Removed, initial=0),
    ])
    numclasses = nrow(abun_h)

    # Set non-pathogen mediated transitions
    sigma = 0.05/day
    transitions = DataFrame([
        (from="Infected", to="Recovered", prob=sigma),
    ])

    birth = [0.0/day; fill(1e-5/day, 3); 0.0/day]
    death = [0.0/day; fill(1e-5/day, 3); 0.0/day]
    beta_force = 5.0/day
    beta_env = 0.5/day
    virus_growth = 0.0001/day
    virus_decay = 0.07/day
    param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force)

    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    dispersal_dists = fill(2.0km, numclasses)
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    @test_nowarn EpiList(traits, abun_v, abun_h, movement, transitions, param)
    epilist = EpiList(traits, abun_v, abun_h, movement, transitions, param)
    @test epilist.virus.names[1] == "Environment"
    @test epilist.virus.names[2] == "Force"
    @test epilist.human.names[1] == "Susceptible"
    @test epilist.human.names[2] == "Infected"
    @test epilist.human.names[3] == "Recovered"
    @test epilist.human.names[4] == "Dead"

    @test length(epilist.human.names) == length(epilist.human.abun)
    @test length(epilist.virus.names) == length(epilist.virus.abun)
    @test length(epilist.virus.names) == length(epilist.virus.traits.mean)
    @test length(epilist.virus.names) == length(epilist.virus.traits.var)
end
