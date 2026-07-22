# SPDX-License-Identifier: LGPL-3.0-or-later

## Small scale test for MPI Ecosystems
using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using MPI
using Random
using Diversity
using JLD2
using Test

nt = Threads.nthreads();
@info "Total Memory: $(Sys.total_memory() / 2^30)GB, threads: $nt"

# Set up MPI and print threads
if !MPI.Initialized()
    MPI.Init()
end

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Set up initial parameters for ecosystem
numSpecies = 8;
grid = (4, 4);
demand = 10.0kJ;
individuals = 1_000;
area = 100.0 * km^2;
totalK = 10000.0kJ / km^2;

# Set up how much resource each species consumes
# resource_vec = SolarDemand(fill(demand, numSpecies))
resource_vec = SolarDemand(collect(1:numSpecies) .* 1.0kJ)

# Set probabilities
birth = 0.6 / year
death = 0.6 / year
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
traits = Bin(MeanTemperature, Normal, opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, resource_vec,
                   movement, param, native)

# Create abiotic environment - even grid of one temperature
habitat = simplehabitatAE(274.0K, grid, totalK, area)
habitat.supply.matrix .= reshape(10_000.0kJ .* collect(1:prod(grid)), grid)

# Set relationship between species and environment (gaussian)
rel = DistRel{typeof(1.0K)}()

# Create ecosystem
@test_nowarn MPIEcosystem(sppl, habitat, rel)
eco = MPIEcosystem(sppl, habitat, rel; seed = 0)

# Artifically fill ecosystem with individuals
eco.abundances.rows_matrix .= 10

# Set columns vector to zero and check synchronise from rows
eco.abundances.cols_vector .= 0
@test_nowarn EcoSISTEM.synchronise_from_rows!(eco.abundances)
@test sum(eco.abundances.cols_vector) == sum(eco.abundances.rows_matrix)

# Set rows matrix to zero and check synchronise from cols
eco.abundances.rows_matrix .= 0
@test_nowarn EcoSISTEM.synchronise_from_cols!(eco.abundances)
@test sum(eco.abundances.cols_vector) == sum(eco.abundances.rows_matrix)

## Reproducibility is provided by the per-species RNG streams seeded in the
## MPIEcosystem constructor (seed = 0 above), independent of the process/thread
## split; the global RNG is not used by the simulation.

# Simulation Parameters
burnin = 2years;
times = 10years;
timestep = 1month;
record_interval = 3months;
repeats = 1;
lensim = length((0years):record_interval:times)

# Burnin
MPI.Barrier(comm)
@test sum(getabundance(eco)) ≈ 1.0
@test sum(getmetaabundance(eco)) ≈ 1.0
@test_nowarn simulate!(eco, burnin, timestep)
@test sum(getabundance(eco)) ≈ 1.0
@test sum(getmetaabundance(eco)) ≈ 1.0

# Collect full abundance matrix together
true_abuns = gather_abundance(eco)
# On root node, print abundances and save out
if rank == 0
    @save joinpath(ARGS[1], "Test_abuns$nt.jld2") abuns=true_abuns
end

water_vec = WaterDemand(fill(2.0mm, numSpecies))
total_use = DemandCollection2(resource_vec, water_vec)

sppl = SpeciesList(numSpecies, traits, abun, total_use, movement, param, native)

# Create abiotic environment - even grid of one temperature
habitat1 = simplehabitatAE(274.0K, grid, totalK, area)

total_mm = 10.0mm / km^2
habitat2 = simplehabitatAE(274.0K, grid, total_mm, area)

supply = SupplyCollection2(habitat1.supply, habitat2.supply)
habitat = GridHabitat{typeof(habitat1.regime), typeof(supply)}(habitat1.regime,
                                                               habitat1.active,
                                                               supply,
                                                               habitat1.names)

# Set relationship between species and environment (gaussian)
rel = DistRel{typeof(1.0K)}()

# Create ecosystem
@test_nowarn MPIEcosystem(sppl, habitat, rel)
eco = MPIEcosystem(sppl, habitat, rel; seed = 0)

# Artifically fill ecosystem with individuals
eco.abundances.rows_matrix .= 10
sleep(rank)

# Set columns vector to zero and check synchronise from rows
eco.abundances.cols_vector .= 0
@test_nowarn EcoSISTEM.synchronise_from_rows!(eco.abundances)
@test sum(eco.abundances.cols_vector) == sum(eco.abundances.rows_matrix)

# Set rows matrix to zero and check synchronise from cols
eco.abundances.rows_matrix .= 0
@test_nowarn EcoSISTEM.synchronise_from_cols!(eco.abundances)
@test sum(eco.abundances.cols_vector) == sum(eco.abundances.rows_matrix)

## Reproducibility is provided by the per-species RNG streams seeded in the
## MPIEcosystem constructor (seed = 0 above), independent of the process/thread
## split; the global RNG is not used by the simulation.

# Simulation Parameters
burnin = 2years;
times = 10years;
timestep = 1month;
record_interval = 3months;
repeats = 1;
lensim = length((0years):record_interval:times)

# Burnin
MPI.Barrier(comm)
@test_nowarn simulate!(eco, burnin, timestep)

sleep(rank)

# Collect full abundance matrix together
true_abuns = gather_abundance(eco)
# On root node, print abundances and save out
if rank == 0
    @save joinpath(ARGS[1], "Test_abuns$nt.jld2") abuns=true_abuns
end

if !MPI.Finalized()
    MPI.Finalize()
end
