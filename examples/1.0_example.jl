## STATIC VERSION
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using RCall
using JLD
using StatsBase
using GLM
using Diversity
using Simulation
using DataFrames
using Distributions
using DataStructures

# Set up Ecosystem as discrete habitat with species having a trait preference
# for one of the two niche types
numSpecies = 150
numTraits = 2
numNiches = 2
numInvasive = 1

sus_mean = 1.0
sus_var = 0.5

size_mean = 5.0
size_var = 15.0

# Set probabilities
birth = fill(0.0/month, numSpecies+numInvasive)
death = fill(0.0/month, numSpecies+numInvasive)
long = 1.0
surv = 0.0
boost = 1000.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = NoGrowth{typeof(unit(0.0/month))}(birth, death, long, surv, boost)

# Set up a large grid of 4km^2 grid squares and enough energy to sustain
# 1 million individuals
grid = (50, 50)
area = 10000.0km^2
totalK = 1000000.0 * (numSpecies + numInvasive)
individuals=20000 * numSpecies

# Create movement type - all individuals are allowed to move and have a wide range
kernel = GaussianKernel(0.0km, numSpecies+numInvasive, 10e-04)
movement = NoMovement(kernel)

# Set up scenario of total habitat loss at certain rate
habloss = 1.0 /10year
declines = 1.0 /10year
scenario = [SimpleScenario(UniformDecline, declines),
    SimpleScenario(ProportionalDecline, declines),
    SimpleScenario(RareDecline, declines),
    SimpleScenario(CommonDecline, declines),
    SimpleScenario(LargeDecline, declines),
    SimpleScenario(Invasive, declines),
    SimpleScenario(Invasive, declines),
    SimpleScenario(SusceptibleDecline, declines),
    SimpleScenario(RandHabitatLoss!, habloss),
    SimpleScenario(ClustHabitatLoss!, habloss)]
divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson,
    mean_abun, geom_mean_abun, pd]
q = 1.0

native = fill(true, numSpecies + numInvasive)
native[numSpecies+numInvasive] = false
pop_mass = rand(Normal(-0.75, 0.1))
sppl = SpeciesList(numSpecies + numInvasive, 2, pop_mass, size_mean,
size_var, area, movement, param, native, [0.5, 0.5])
Simulation.resettraits!(sppl.types.tree)
sppl.susceptible = ContinuousEvolve(sus_mean, sus_var, sppl.types.tree).mean
abenv = simplenicheAE(numNiches, grid, totalK, area)
rel = Match{eltype(abenv.habitat)}()
eco = Ecosystem(trait_populate!, sppl, abenv, rel)

function runsim(times::Unitful.Time)
    interval = 1month; reps = 1
    lensim = length(0month:interval:times)
    abun = zeros(1, length(divfuns), lensim, length(scenario), reps)
    for j in 1:reps
        for i in 1:length(scenario)
            native = fill(true, numSpecies + numInvasive)
            native[numSpecies+numInvasive] = false
            pop_mass = rand(Normal(-0.75, 0.1))
            sppl = SpeciesList(numSpecies + numInvasive, 2, pop_mass, size_mean,
            size_var, area, movement, param, native, [0.5, 0.5])
            Simulation.resettraits!(sppl.types.tree)
            sppl.susceptible = ContinuousEvolve(sus_mean, sus_var, sppl.types.tree).mean
            if i == 7
                reroot!(sppl.types.tree, "151")
            end
            abenv = simplenicheAE(numNiches, grid, totalK, area)
            rel = Match{eltype(abenv.habitat)}()
            eco = Ecosystem(trait_populate!, sppl, abenv, rel)
            thisabun = view(abun, :, :, :, i, j);
            simulate_record_diversity!(thisabun, eco, times, interval, timestep,
            scenario[i], divfuns, q)
        end
    end
    abun
end


times = 1year
div = runsim(times)


using Simulation
using JLD
using Diversity
using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using AxisArrays
using Simulation.ClimatePref

# Load in temperature profiles
Temp = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/Temperature.jld",
 "Temperature")

 ## Run simulation over a grid and plot
 numSpecies=3

 # Set up how much energy each species consumes
 energy_vec = SolarRequirement(fill(2.0*kJ, numSpecies))

 # Set probabilities
 birth = 0.6/month
 death = 0.6/month
 long = 1.0
 surv = 0.5
 boost = 1000.0
 timestep = 1.0month

 # Collect model parameters together (in this order!!)
 param = EqualPop(birth, death, long, surv, boost)

 grid = (94, 60)
 area = 10000.0km^2
 totalK = 1000000.0
 individuals=1000

 # Load data for land cover
 file = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/World.tif"
 world = extractfile(file)
 europe = world[-10° .. 60°, 35° .. 80°]
 eu = ustrip.(europe)
 @rput eu
 R"image(eu)"

 dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_1980"
 tempax1 = extractERA(dir1, "t2m", collect(1.0month:1month:10year))
 dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_1990"
 tempax2 = extractERA(dir2, "t2m", collect(121month:1month:20year))
 dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_2000"
 tempax3 = extractERA(dir3, "t2m", collect(241month:1month:30year))
 dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_2010"
 tempax4 = extractERA(dir4, "t2m", collect(361month:1month:38year))


 dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
 srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
 srad = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3}, srad.array[-10° .. 60°, 35° .. 80°,:])
 srad = SolarTimeBudget(srad, 1)
 active = Array{Bool, 2}(.!isnan.(eu))


 testtemp = tempax1
 testtemp.array = tempax1.array[-10° .. 60°, 35° .. 80°, :]
 # Create ecosystem
 kernel = GaussianKernel(10.0km, numSpecies, 10e-4)
 movement = BirthOnlyMovement(kernel)


 native = fill(true, numSpecies)
 traits = Array(transpose(Temp[1:3,:]))
 traits = TempBin(traits)
 abun = rand(Multinomial(individuals, numSpecies))
 sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
 movement, param, native)
 abenv = eraAE(testtemp, srad, active)
 #abenv = tempgradAE(-10.0°C, 30.0°C,  grid, totalK, area, 0.0°C/month, active)
 #abenv.habitat.matrix = transpose(tempax1.array[-10° .. 60°, 35° .. 80°, 1])
 rel = Trapeze{eltype(abenv.habitat)}()
 eco = Ecosystem(sppl, abenv, rel)


 function runsim(eco::Ecosystem, times::Unitful.Time)
     burnin = 1year; interval = 1month
     lensim = length(0month:interval:times)
     abun = generate_storage(eco, lensim, 1)
     simulate!(eco, burnin, timestep)
     eco.abenv.habitat.time = 1
     #resetrate!(eco, 0.1°C/month)
     simulate_record!(abun, eco, times, interval, timestep)
 end
 times = 1year
 abun = runsim(eco, times)
