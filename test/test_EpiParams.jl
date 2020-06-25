using Simulation
using Unitful.DefaultSymbols
using Test
using Simulation.Units
import Simulation.SIRGrowth

birth = [fill(1e-5/day, 3); 0.0/day]
death = [fill(1e-5/day, 3); 0.0/day]
beta_force = 5.0/day
beta_env = 0.5/day
sigma = 0.05/day
virus_growth = 0.0001/day
virus_decay = 0.07/day
@test_nowarn SIRGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
param = SIRGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
@test length(param.birth) == length(param.death)
@test param.birth[:, 1] == birth
@test param.death[:, 1] == death
@test param.beta_force[1] == beta_force
@test param.beta_env[1] == beta_env
@test param.sigma[1] == sigma

transition_params = transition(param)
@test size(transition_params.transition) == size(transition_params.transition_virus)
@test transition_params.transition[end, :] == death
