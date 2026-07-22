# SPDX-License-Identifier: LGPL-3.0-or-later

#### SINGLE SPECIES ####
# Code to run single species across Africa with WorldClim data.
using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using RasterDataSources
using AxisArrays
using Unitful
using Unitful.DefaultSymbols
using StatsBase
using Plots
using DataPipeline

# Initialise datapipeline
handle = DataPipeline.initialise()

# Download temperature and precipitation data
path = link_read!(handle, "AfricaModel/WorldClim")
newpath = EcoSISTEM.unziptemp(path)
world = read(WorldClim{BioClim}, [1, 12];
             cut = EcoSISTEM.ClimatePref.boundingbox("Africa"; round = 5°))
# `world.array` is a (latitude, longitude, layer) stack cut to Africa: index 1 = bio1 (annual mean
# temperature), index 2 = bio12 (annual precipitation).
africa_temp = world.array[:, :, 1]
bio_africa = uconvert.(K, africa_temp .* °C)
bio_africa = ClimateRaster(WorldClim{BioClim},
                           AxisArray(bio_africa,
                                     AxisArrays.axes(africa_temp)))
africa_water = world.array[:, :, 2] .* mm
africa_water = ClimateRaster(WorldClim{BioClim},
                             AxisArray(africa_water,
                                       AxisArrays.axes(africa_temp)))
bio_africa_water = WaterSupply(africa_water)

# Find which grid cells are land
active = Matrix{Bool}(.!isnan.(bio_africa.array))

heatmap(africa_temp)

# Set up initial parameters for ecosystem
numSpecies = 1;
grid = size(active);
dem = 0.1mm;
individuals = 0;
area = 64e6km^2;
totalK = 1000.0kJ / km^2;

# Set up how much water each species consumes
resource_vec = WaterDemand(fill(dem, numSpecies))

# Set rates for birth and death
birth = 0.6 / year
death = 0.6 / year
longevity = 1.0
survival = 0.2
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
movement = AlwaysMovement(kernel, Torus())

# Create species list, including their temperature preferences, seed abundance and native status
opts = fill(280.0K, numSpecies)
vars = fill(10.0K, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, resource_vec,
                   movement, param, native)

# Create abiotic environment - with temperature and water resource
abenv = bioclimAE(bio_africa, bio_africa_water, active)

# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

# Create ecosystem and fill every active grid square with an individual
eco = Ecosystem(sppl, abenv, rel)
rand_start = findall(active)
for i in rand_start
    eco.abundances.grid[1, i[1], i[2]] += 1
end

# Run simulation
times = 10years;
timestep = 1month;
record_interval = 1month;
repeats = 1;
lensim = length((0years):record_interval:times)
abuns = zeros(Int64, numSpecies, prod(grid), lensim)
@time simulate_record!(abuns, eco, times, record_interval, timestep);

# Reshape abundances for plotting
abuns = reshape(abuns[1, :, :, 1], grid[1], grid[2], lensim)

# Plot start and end abundances, next to temperature and rainfall
africa_startabun = Float64.(abuns[:, :, 1])
africa_startabun[.!(active)] .= NaN
africa_endabun = Float64.(abuns[:, :, end])
africa_endabun[.!(active)] .= NaN
heatmap(africa_startabun,
        clim = (0, maximum(abuns)),
        background_color = :lightblue,
        background_color_outside = :white,
        grid = false,
        color = cgrad(:algae, scale = :exp),
        layout = (@layout [a b; c d]))
heatmap!(africa_endabun,
         clim = (0, maximum(abuns)),
         background_color = :lightblue,
         background_color_outside = :white,
         grid = false,
         color = cgrad(:algae, scale = :exp),
         subplot = 2)

africa_temp = world.array[:, :, 1]
africa_water = world.array[:, :, 2]
heatmap!(africa_temp, grid = false, subplot = 3)
heatmap!(africa_water, grid = false, subplot = 4)

path = link_write!(handle, "Africa-plot")
Plots.pdf(path)

DataPipeline.finalise(handle)
