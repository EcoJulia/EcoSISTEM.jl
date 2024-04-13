using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using DataFrames

# Set up simple gridded environment
grid = (10, 10)
area = 1_000.0km^2
epienv = simplehabitatAE(298.0K, grid, area, NoControl())

# Set initial population sizes for all pathogen categories
abun_v = DataFrame([
                       (name = "Environment", initial = 0),
                       (name = "Force", initial = 0)
                   ])
numvirus = nrow(abun_v)

# Set initial population sizes for all human categories
susceptible = 500_000 * prod(grid)
exposed = 0
infected = 100 * prod(grid)
abun_h = DataFrame([
                       (name = "Susceptible", type = Susceptible,
                        initial = susceptible),
                       (name = "Exposed", type = OtherDiseaseState,
                        initial = exposed),
                       (name = "Infected", type = Infectious,
                        initial = infected),
                       (name = "Recovered", type = Removed, initial = 0),
                       (name = "Dead", type = Removed, initial = 0)
                   ])
numclasses = nrow(abun_h)

# Set non-pathogen mediated transitions
mu = 1 / 7days
sigma = 1 / 7days

# Set simulation parameters
birth = fill(0.0 / day, numclasses)
death = fill(0.0 / day, numclasses)
beta_force = 10.0 / day
beta_env = 10.0 / day
virus_growth = 0.1 / day
virus_decay = 1.0 / 2days
param = (birth = birth, death = death, virus_growth = virus_growth,
         virus_decay = virus_decay, beta_env = beta_env,
         beta_force = beta_force)

transitiondat = DataFrame([
                              (from = "Susceptible", from_id = 1,
                               to = "Exposed", to_id = 2,
                               prob = (env = beta_env, force = beta_force)),
                              (from = "Exposed", from_id = 2, to = "Infected",
                               to_id = 3, prob = mu),
                              (from = "Infected", from_id = 3, to = "Recovered",
                               to_id = 4, prob = sigma)
                          ])

# Dispersal kernels for virus and disease classes
dispersal_dists = fill(10.0km, prod(grid))
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = EpiMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = SpeciesList(traits, abun_v, abun_h, movement, transitiondat, param)
rel = Gauss{eltype(epienv.habitat)}()

# Create list of transitions for the simulation
transitions = TransitionList()
addtransition!(transitions, UpdateEpiEnvironment(update_epi_environment!))
for loc in eachindex(epienv.habitat.matrix)
    addtransition!(transitions, ForceProduce(3, loc, param.virus_growth))
    addtransition!(transitions, ViralLoad(loc, param.virus_decay))
    addtransition!(transitions,
                   Exposure(transitiondat[1, :from_id], loc,
                            transitiondat[1, :to_id],
                            transitiondat[1, :prob].force,
                            transitiondat[1, :prob].env))
    addtransition!(transitions,
                   Infection(transitiondat[2, :from_id], loc,
                             transitiondat[2, :to_id], transitiondat[2, :prob]))
    addtransition!(transitions,
                   Recovery(transitiondat[3, :from_id], loc,
                            transitiondat[3, :to_id], transitiondat[3, :prob]))
    for spp in eachindex(epilist.species.names)
        addtransition!(transitions, ForceDisperse(spp, loc))
    end
end

# Create epi system with all information
epi = Ecosystem(epilist, epienv, rel, transitions = transitions)

# Run simulation
times = 1month;
interval = 1day;
timestep = 1day;
abuns = zeros(Int64, numclasses, prod(grid), floor(Int, times / timestep) + 1)
@test_nowarn simulate_record!(abuns, epi, times, interval, timestep);

for su in epi.transitions.setup
    @test_nowarn run_rule!(epi, su, timestep)
end
for st in epi.transitions.state
    @test_nowarn run_rule!(epi, st, timestep)
end
for pl in epi.transitions.place
    @test_nowarn run_rule!(epi, pl, timestep)
end
for wd in epi.transitions.winddown
    @test_nowarn run_rule!(epi, wd, timestep)
end
