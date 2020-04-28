using Simulation
using Unitful.DefaultSymbols
using Test
using Simulation.Units
import Simulation.SIRGrowth

birth = [0.0/day; fill(1e-5/day, 3); 0.0/day]
death = [0.0/day; fill(1e-5/day, 3); 0.0/day]
beta = 0.0005/day
sigma = 0.05/day
virus_growth = 0.0001/day
virus_decay = 0.07/day
@test_nowarn SIRGrowth{typeof(unit(beta))}(birth, death, virus_growth, virus_decay, beta, sigma)
param = SIRGrowth{typeof(unit(beta))}(birth, death, virus_growth, virus_decay, beta, sigma)
@test length(param.birth) == length(param.death)
@test param.birth == birth
@test param.death == death
@test param.beta == beta
@test param.sigma == sigma

transition_params = transition(param)
@test size(transition_params.transition) == size(transition_params.transition_virus)
@test length(transition_params.virus_decay) == length(transition_params.virus_growth)
@test transition_params.transition[2, :] == birth
@test transition_params.transition[2, :] == death
