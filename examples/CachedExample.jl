using Simulation
using JLD
using Diversity
using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using AxisArrays
#using ClimatePref

# Load in temperature profiles
Temp = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/Temperature.jld",
 "Temperature")
 Rain = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/Rainfall.jld",
  "Rainfall")

## Run simulation over a grid and plot
numSpecies=3

# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(repmat([2*day^-1*kJ*m^-2], numSpecies))
energy_vec2 = WaterRequirement(repmat([0.5*mm], numSpecies))
energy_vec = ReqCollection2(energy_vec1, energy_vec2)

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

dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_1980"
tempax1 = extractERA(dir1, "t2m", collect(1.0month:1month:10year))

dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
prec = extractworldclim(joinpath(dir, "wc2.0_5m_prec"))
prec.array = prec.array[-10° .. 60°, 35° .. 80°,:]
x = prec.array.axes[1]
y = prec.array.axes[2]
t = prec.array.axes[3]
prec.array = AxisArray(1.0.*(prec.array),
    Axis{:longitude}(x), Axis{:latitude}(y), Axis{:time}(t))
water = WaterBudget(Array{typeof(1.0*mm), 3}(prec.array), 1)


dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3},
            srad.array[-10° .. 60°, 35° .. 80°,:])
srad = SolarBudget(srad, 1)
active = Array{Bool, 2}(.!isnan.(eu))
bud = BudgetCollection2(srad, water)

testtemp = tempax1
testtemp.array = tempax1.array[-10° .. 60°, 35° .. 80°, :]

# Create ecosystem
kernel = GaussianKernel(10.0km, numSpecies, 10e-4)
movement = BirthOnlyMovement(kernel)

common_species = ["Trifolium repens", "Urtica dioica", "Achillea millefolium"]
native = repmat([true], numSpecies)
traits1 = Array(transpose(Temp[common_species,:]))
traits1 = TempBin(traits1)
traits2 = Array(transpose(Rain[common_species,:]))
traits2 = RainBin(traits2)
abun = Multinomial(individuals, numSpecies)
sppl = SpeciesList(numSpecies, TraitCollection2(traits1, traits2), abun, energy_vec,
movement, param, native)
abenv1 = eraAE(testtemp, srad, active)
abenv2 = worldclimAE(prec, srad, active)
hab = HabitatCollection2(abenv1.habitat, abenv2.habitat)
abenv = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, abenv1.active,
    bud, abenv1.names)
rel1 = Trapeze{eltype(abenv.habitat.h1)}()
rel2 = Unif{eltype(abenv.habitat.h2)}()
rel = multiplicativeTR2(rel1, rel2)
eco = Ecosystem(sppl, abenv, rel)

times = 1month:1month:10years
cache = CachedEcosystem(eco, "examples/Cache", times)
div = norm_sub_alpha(cache, 1.0)
for i in times[2:end]
    abundances(cache, i)
    newdiv = norm_sub_alpha(cache, 1.0)
    append!(div, newdiv)
end
div[:time] = vcat(map(x -> fill(x, size(eco.abundances.matrix, 2)), times) ...)
