using Pkg
Pkg.activate("examples")

using Simulation
using MyUnitful
using Unitful, Unitful.DefaultSymbols
using Distributions
using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
println(Threads.nthreads())
numSpecies = 100; grid = (50, 50); req= 10.0kJ; individuals=1000; area = 100000.0*km^2; totalK = 100.0kJ/km^2
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

# # Check synchronisation is working
# m1=deepcopy(eco.abundances.rows_matrix)
# Simulation.synchronise_from_rows!(eco.abundances);
# Simulation.synchronise_from_cols!(eco.abundances);
# m2=deepcopy(eco.abundances.rows_matrix)
# Simulation.synchronise_from_rows!(eco.abundances);
# eco.abundances.rows_matrix .= 0;
# Simulation.synchronise_from_cols!(eco.abundances);
# m3=deepcopy(eco.abundances.rows_matrix)
# eco.abundances.cols_vector .= 0;
# Simulation.synchronise_from_rows!(eco.abundances);
# Simulation.synchronise_from_cols!(eco.abundances);
# m4=deepcopy(eco.abundances.rows_matrix)
# Simulation.synchronise_from_rows!(eco.abundances);
# eco.abundances.reshaped_cols[1] .= 0;
# if MPI.Comm_size(MPI.COMM_WORLD) > 1
#     eco.abundances.reshaped_cols[2] .= 1;
# end
# Simulation.synchronise_from_cols!(eco.abundances);
# m5=deepcopy(eco.abundances.rows_matrix)
#
# # Check growth
# sleep(rank)
# println("$(rank):")
# println(all(m1 .== m2))
# println(all(m2 .== m3))
# println(all(m3 .== m4))
# if rank == 0
#     println(all(m5 .== 0))
# elseif rank == 1
#     println(all(m5 .== 1))
# else
#     println(all(m4 .== m5))
# end
#
# eco = MPIEcosystem(sppl, abenv, rel)
sleep(rank)
print("$(rank):")
println(eco.abundances.rows_matrix)

# Simulation Parameters
burnin = 5years; times = 50years; timestep = 1month; record_interval = 3months; repeats = 1
lensim = length(0years:record_interval:times)
# Burnin
MPI.Barrier(comm)
one = time_ns()
@time simulate!(eco, burnin, timestep)
two = time_ns()

sleep(rank + 1)
println("$(rank): $((two - one) / 10 ^ 9)")

print("$(rank):")
println(eco.abundances.rows_matrix)
MPI.Finalize()
