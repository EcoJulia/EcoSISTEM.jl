using Pkg
Pkg.activate("examples")

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
println(Threads.nthreads())
numSpecies = 100; grid = (50, 50); req= 10.0kJ; individuals=1_000_000; area = 100000.0*km^2; totalK = 100.0kJ/km^2
# Set up initial parameters for ecosystem

# Set up how much energy each species consumes
energy_vec = SolarRequirement(fill(req, numSpecies))

# Set probabilities
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.2
boost = 1.0

# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(1.0km, 10e-10), numSpecies)
movement = BirthOnlyMovement(kernel, Torus())

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

# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

# Create ecosystem
eco = MPIEcosystem(sppl, abenv, rel)

sleep(rank)
print("$(rank):")
println(sum(eco.abundances.rows_matrix, dims = 2))

# Simulation Parameters
burnin = 20years; times = 50years; timestep = 1month; record_interval = 3months; repeats = 1
lensim = length(0years:record_interval:times)
# Burnin
MPI.Barrier(comm)
one = time_ns()
@time simulate!(eco, burnin, timestep)
two = time_ns()

sleep(rank + 1)
println("$(rank): $((two - one) / 10 ^ 9)")

print("$(rank):")
println(sum(eco.abundances.rows_matrix, dims = 2))
MPI.Finalize()
