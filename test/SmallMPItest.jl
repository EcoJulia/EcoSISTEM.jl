## Small scale test for MPI Ecosystems
## Set random number seeds for reproducible results

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using MPI
using Random
using Diversity
using JLD2
using Test

# Set up MPI and print threads
MPI.Init()
nthread = Threads.nthreads()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
println(Threads.nthreads())

# Set up initial parameters for ecosystem
numSpecies = 8; grid = (2, 2); req= 10.0kJ; individuals=1_000; area = 100.0*km^2; totalK = 10000.0kJ/km^2

# Set up how much energy each species consumes
# energy_vec = SolarRequirement(fill(req, numSpecies))
energy_vec = SolarRequirement(collect(1:numSpecies) .* 1.0kJ)

# Set probabilities
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.2
boost = 100.0

# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(3.0km, 10e-10), numSpecies)
movement = BirthOnlyMovement(kernel, NoBoundary())

# Create species list, including their temperature preferences, seed abundance and native status
opts = fill(274.0K, numSpecies)
vars = fill(0.5K, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)

# Create abiotic environment - even grid of one temperature
abenv = simplehabitatAE(274.0K, grid, totalK, area)
abenv.budget.matrix .= reshape(10_000.0kJ .* collect(1:prod(grid)), grid)

# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

# Create ecosystem
eco = MPIEcosystem(sppl, abenv, rel)

# Artifically fill ecosystem with individuals
eco.abundances.rows_matrix .= 10
sleep(rank)

# Set columns vector to zero and check synchronise from rows
eco.abundances.cols_vector .= 0
EcoSISTEM.synchronise_from_rows!(eco.abundances)
@test sum(eco.abundances.cols_vector) == sum(eco.abundances.rows_matrix)

# Set rows matrix to zero and check synchronise from cols
eco.abundances.rows_matrix .= 0
EcoSISTEM.synchronise_from_cols!(eco.abundances)
@test sum(eco.abundances.cols_vector) == sum(eco.abundances.rows_matrix)

# Set random seed to 0 for all rngs
for i in 1:length(eco.abundances.rngs)
    eco.abundances.rngs[i] = MersenneTwister(0)
end


# Simulation Parameters
burnin = 2years; times = 10years; timestep = 1month; record_interval = 3months; repeats = 1
lensim = length(0years:record_interval:times)

# Burnin
MPI.Barrier(comm)
@time simulate!(eco, burnin, timestep)

sleep(rank)

# Collect full abundance matrix together
true_abuns = gather_abundance(eco)
# On root node, print abundances and save out
if rank == 0
    print("$(rank):")
    println(true_abuns)
    JLD2.save("data/Test_abuns"*"$nthread.jld2", "abuns", true_abuns)
end

MPI.Finalize()
