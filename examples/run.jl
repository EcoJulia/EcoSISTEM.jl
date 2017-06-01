# Include simulation functions
using Simulation
using Distributions
using RCall

## Run simulation over a grid and plot
numSpecies=4
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = SimpleRequirement(repmat([2], numSpecies))

# Set probabilities
birth = 0.6
death = 0.6
l = 1.0
s = 0.0
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, timestep, l, s]

grid = (5,5)
gridSize = 1.0
totalK = 1000.0
individuals=100

# Create ecosystem
movement = GaussianMovement(0.2, numSpecies, 10e-4)
sppl = SpeciesList(numSpecies, numTraits, Multinomial(individuals, numSpecies),
                   energy_vec, movement)
abenv = simplenicheAE(numNiches, grid, totalK, gridSize)
eco = Ecosystem(sppl, abenv, false)
plot_move(eco, 2, 2, 1)
# Run simulation grid
abun = run_sim_spatial(eco, param, 100, 1, 1, 1, false)

plot_abun(abun, numSpecies, grid[1])
