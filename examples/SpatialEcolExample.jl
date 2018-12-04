# Spatial Ecology plotting examples

using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using Diversity
using Distributions
using Simulation
using EcoBase
using SpatialEcology
using Plots
using DataFrames

# Set up Ecosystem as discrete habitat with species having a trait preference
# for one of the two niche types
numSpecies = 150
numNiches = 2

# Demographics
birth = 0.0/month
death = 0.0/month
long = 1.0
surv = 0.0
boost = 1000.0
timestep = 1.0month
param = EqualPop(birth, death, long, surv, boost)

# Abiotic Environment
grid = (50, 50)
area = 10000.0km^2
totalK = 1000000.0 * numSpecies
abenv = simplenicheAE(numNiches, grid, totalK, area)

# Draw abundances from Log normal
individuals = 20000 * numSpecies
probs = rand(LogNormal(10.0, 1.5), numSpecies)
probs /= sum(probs)
abun = rand(Multinomial(individuals, probs))

# Species characteristics
kernel = GaussianKernel(1.0km, numSpecies, 10e-04)
movement = BirthOnlyMovement(kernel)
native = fill(true, numSpecies)
energy = SimpleRequirement(fill(2.0, numSpecies))
sppl = SpeciesList(numSpecies, numNiches, abun, energy, movement, param, native)

# Set up relationship between habitat and species
rel = Match{eltype(abenv.habitat)}()
eco = Ecosystem(trait_populate!, sppl, abenv, rel)

# Plotting ecosystem results in error: 'type Ecosystem has no field site'
plot(eco)

## Can plot separately by converting dataframe of abundances into an assemblage

# Normalised alpha at q = 1
alpha_div = norm_sub_alpha(eco, 1.0)

# Convert ecosystem abundances into dataframe
abun_data = DataFrame(transpose(eco.abundances.matrix))
# Convert grid cell positions to x,y coordinates (making sure they are Float values
# because otherwise they are interpreted as point locations)
coords = convert_coords.(parse.(alpha_div[:partition_name]), grid[1])
abun_data[:lon] = [ x[1] for x in coords ] * 1.0
abun_data[:lat] = [ x[2] for x in coords ] * 1.0
abun_data[:coords] = map((x,y) ->"$x" * "_$y", abun_data[:lon], abun_data[:lat])

# Now plot as an assemblage
assem = Assemblage(abun_data[1:150], abun_data[151:153], sitecolumns = false)
plot(assem)

# Plot alpha diversity on top of this
plot(alpha_div[:diversity], assem, color = :fire)
