using Traceur
using Simulation
using Distributions
using Unitful
using Unitful.DefaultSymbols
using ClimatePref
using MyUnitful
using AxisArrays


# Set up initial parameters for ecosystem
numSpecies = 1

# Set up how much energy each species consumes
energy_vec = SolarRequirement([1000 .* kJ])

# Set probabilities
birth = 0.6/year
death = 0.6/year
l = 1.0
s = 0.0
boost = 1.0
timestep = 1month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, l, s, boost)

grid = (20, 20)
area = 4.0km^2
dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad.array = srad.array[-10° .. 60°, 35° .. 80°]
meansrad = mean(srad.array[.!isnan.(srad.array)])
totalK = uconvert(kJ, meansrad * month * (area/(grid[1]*grid[2])))
individuals = 100000

# Create ecosystem
kernel = GaussianKernel(0.5km, numSpecies, 10e-10)
movement = BirthOnlyMovement(kernel, Torus())

opts = rand(Normal(274.0, 10.0), numSpecies) * K
vars = rand(Uniform(0, 25/9), numSpecies) * K
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
abun = rand(Multinomial(individuals, numSpecies))
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
abenv = simplehabitatAE(274.0K, grid, totalK, area)
rel = Gauss{typeof(1.0K)}()
eco = Ecosystem(sppl,abenv,rel)


@trace Simulation.calc_lookup_moves(eco.spplist.movement.boundary, 1, 1, 1, eco, 100)
