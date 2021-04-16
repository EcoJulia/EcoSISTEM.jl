start = time()
using UKclim
using JuliaDB
using JuliaDBMeta
using BritishNationalGrid
using Unitful
using EcoSISTEM.Units
using Unitful.DefaultSymbols
using AxisArrays
using Statistics
using EcoSISTEM
using Distributions
using Diversity
using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
rank == 0 && println("using: $((time() - start) * s)")
totMPI = MPI.Comm_size(comm)
io = open("output-cores$(totMPI*Threads.nthreads())-np$totMPI-$rank.txt", write = true)
write(io, "rank $rank / $totMPI: $(Threads.nthreads()) threads @ $((time() - start) * s)\n")
close(io)

traits = load("../data/BSBI_had_prefs_UK")
traits = filter(t-> !isnan(t.sun) & !isnan(t.rainfall) & !isnan(t.tas_st) & !isnan(t.rain_st), traits)
traits = filter(t -> (t.rain_st > 0) & (t.tas_st > 0), traits)
numSpecies = length(traits)
individuals = Int(1e10)

# Set up species requirements
sizes = abs.(rand(Normal(1.0, 0.1), numSpecies)) .* m^2
solarreq = collect(select(traits, :sun)) .* (kJ/km^2)
req1 = SolarRequirement(uconvert.(kJ, solarreq .* sizes))

waterreq = collect(select(traits, :rainfall)) .* (mm/km^2)
req2 = WaterRequirement(uconvert.(mm, waterreq .* sizes))

req = ReqCollection2(req1, req2)

tmean = collect(select(traits, :tas)) .* K
tsd = collect(select(traits, :tas_st)) .* K
tsd .+= 1e-3K
temp_traits = GaussTrait(tmean, tsd)

pmean = collect(select(traits, :rainfall)) .* mm
psd = collect(select(traits, :rain_st)) .* mm
psd .+= 1e-3mm
prec_traits = GaussTrait(pmean, psd)

av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)
movement = BirthOnlyMovement(kernel, NoBoundary())

abun = rand(Multinomial(individuals, numSpecies))
trts = TraitCollection2(temp_traits, prec_traits)

death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
birth_rates = death_rates
param = PopGrowth{typeof(unit(birth_rates[1]))}(birth_rates, death_rates, 1.0, 1e-3, 1.0)

native = fill(true, numSpecies)

sppl = SpeciesList(numSpecies, trts, abun, req, movement, param, native)


# Import HadUK grid data
dir = "../data/HadUK/tas/"
times = collect(2010year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "../data/HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "../data/HadUK/sun/"
sun = readHadUK(dir, "sun", times)

# Create abiotic environment
solar = uconvert.(kJ, 1km^2 .* sun.array .* 1000.0*(W/m^2))
sol = SolarTimeBudget(solar, 1)

rain = Array(rainfall.array)
water = WaterTimeBudget(rain, 1)

bud = BudgetCollection2(sol, water)

active = Array{Bool, 2}(.!isnan.(sun.array[:, :, 1]))
temp = hadAE(tas, sol, active)
rain = hadAE(rainfall, sol, active)
#lcae = lcAE(lc, 1000.0kJ/km^2, 242_495km^2)

hab = HabitatCollection2(temp.habitat, rain.habitat)
ae = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, temp.active, bud, temp.names)

rel1 = Gauss{eltype(ae.habitat.h1)}()
rel2 = Gauss{eltype(ae.habitat.h2)}()
rel = multiplicativeTR2(rel1, rel2)
eco = MPIEcosystem(sppl, ae, rel)

burnin = 2months; timestep = 1month; times = 1year
io = open("output-cores$(totMPI*Threads.nthreads())-np$totMPI-$rank.txt", append = true)
sppcounts = sum(eco.abundances.rows_matrix, dims = 2)
write(io, "numspp = $(size(eco.abundances.rows_matrix, 1)) " *
        "mean = $(mean(sppcounts)) std = $(std(sppcounts)) @ $((time() - start) * s)\n")
close(io)
rank == 0 && println("Start first burnin: $((time() - start) * s)")
MPI.Barrier(comm)
one = time_ns()
val = simulate!(eco, burnin, timestep)
two = time_ns()

io = open("output-cores$(totMPI*Threads.nthreads())-np$totMPI-$rank.txt", append = true)
write(io, "time = $(convert(typeof(1.0s), (two - one) * ns)) @ $((time() - start) * s)\n")
sppcounts = sum(eco.abundances.rows_matrix, dims = 2)
write(io, "numspp = $(size(eco.abundances.rows_matrix, 1)) mean = $(mean(sppcounts)) std = $(std(sppcounts))\n")
close(io)
rank == 0 && println("Start second burnin: $((time() - start) * s)")
MPI.Barrier(comm)
one = time_ns()
val = @timed simulate!(eco, times, timestep)
two = time_ns()

io = open("output-cores$(totMPI*Threads.nthreads())-np$totMPI-$rank.txt", append = true)
write(io, "$val @ $((time() - start) * s)\n")
write(io, "time = $(convert(typeof(1.0s), (two - one) * ns))\n")
sppcounts = sum(eco.abundances.rows_matrix, dims = 2)
write(io, "numspp = $(size(eco.abundances.rows_matrix, 1)) mean = $(mean(sppcounts)) std = $(std(sppcounts))\n")
close(io)
rank == 0 && println("End second burnin: $((time() - start) * s)")
MPI.Finalize()
