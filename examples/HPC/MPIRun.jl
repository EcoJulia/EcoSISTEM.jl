@info "Total Memory: $(Sys.total_memory() / 2^30)GB"
@info "Num threads: $(Threads.nthreads())"

start = time()
using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using MPI
using Statistics

const SAVEDIR = "/mnt/data/project0000/outputs/MPIRun"
mkpath(SAVEDIR)

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
rank == 0 && println("using: $((time() - start) * s)")
totMPI = MPI.Comm_size(comm)
io = open(joinpath(SAVEDIR, "output-cores$(totMPI*Threads.nthreads())-np$totMPI-$rank.txt"),
          write = true)
write(io,
      "rank $rank / $totMPI: $(Threads.nthreads()) threads @ $((time() - start) * s)\n")
close(io)

# Set up initial parameters for ecosystem
numSpecies = 100_000;
grid = (100, 100);
req = 1.0kJ;
individuals = 1_000_000_000;
area = 1_000_000.0 * km^2;
totalK = 1000.0kJ / km^2;

# Set up how much energy each species consumes
energy_vec = SolarRequirement(fill(req, numSpecies))

# Set probabilities
birth = 0.6 / year
death = 0.6 / year
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

rank == 0 && println("Startup: $((time() - start) * s)")

# Create ecosystem
eco = MPIEcosystem(sppl, abenv, rel)

io = open(joinpath(SAVEDIR, "output-cores$(totMPI*Threads.nthreads())-np$totMPI-$rank.txt"),
          append = true)
sppcounts = sum(eco.abundances.rows_matrix, dims = 2)
write(io,
      "numspp = $(size(eco.abundances.rows_matrix, 1)) " *
      "mean = $(mean(sppcounts)) std = $(std(sppcounts)) @ $((time() - start) * s)\n")
close(io)

# Simulation Parameters
burnin = 1years;
times = 50years;
timestep = 1month;
record_interval = 3months;
repeats = 1;
lensim = length((0years):record_interval:times)
# Burnin
rank == 0 && println("Start first burnin: $((time() - start) * s)")
MPI.Barrier(comm)
one = time_ns()
val = simulate!(eco, burnin, timestep)
two = time_ns()

io = open(joinpath(SAVEDIR, "output-cores$(totMPI*Threads.nthreads())-np$totMPI-$rank.txt"),
          append = true)
write(io,
      "time = $(convert(typeof(1.0s), (two - one) * ns)) @ $((time() - start) * s)\n")
sppcounts = sum(eco.abundances.rows_matrix, dims = 2)
write(io,
      "numspp = $(size(eco.abundances.rows_matrix, 1)) mean = $(mean(sppcounts)) std = $(std(sppcounts))\n")
close(io)
rank == 0 && println("Start second burnin: $((time() - start) * s)")
MPI.Barrier(comm)
one = time_ns()
val = @timed simulate!(eco, burnin, timestep)
two = time_ns()

io = open(joinpath(SAVEDIR, "output-cores$(totMPI*Threads.nthreads())-np$totMPI-$rank.txt"),
          append = true)
write(io, "$val @ $((time() - start) * s)\n")
write(io, "time = $(convert(typeof(1.0s), (two - one) * ns))\n")
sppcounts = sum(eco.abundances.rows_matrix, dims = 2)
write(io,
      "numspp = $(size(eco.abundances.rows_matrix, 1)) mean = $(mean(sppcounts)) std = $(std(sppcounts))\n")
close(io)
rank == 0 && println("End second burnin: $((time() - start) * s)")
MPI.Finalize()
