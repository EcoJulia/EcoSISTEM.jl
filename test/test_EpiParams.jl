using EcoSISTEM
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using DataFrames

@testset "Epi params" begin
    numclasses = 4
    age_categories = 1
    birth = fill(1e-5/day, numclasses)
    death = fill(1e-5/day, numclasses)
    beta_force = 5.0/day
    beta_env = 0.5/day
    sigma = 0.05/day
    virus_growth = 0.0001/day
    virus_decay = 0.07/day
    param = (birth = birth, death = death, virus_growth = virus_growth, virus_decay = virus_decay, beta_env = beta_env, beta_force = beta_force)
    paramDat = DataFrame([(from="Infected", to="Recovered", prob=sigma,
                           from_ind=2, to_ind=3)])
    transition_params = transition(param, paramDat, numclasses, [2], age_categories)
    @test size(transition_params.transition) == size(transition_params.transition_virus)
    @test transition_params.transition[end, :] == death

    @test transition_params.births == birth
    @test transition_params.transition_force[2, 1] == beta_force
    @test transition_params.transition_virus[2, 1] == beta_env
    @test transition_params.transition[3, 2] == sigma
end
