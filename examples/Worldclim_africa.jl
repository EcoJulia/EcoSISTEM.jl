#### SINGLE SPECIES ####
# Code to run single species across Africa with WorldClim data.
using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using SimpleSDMLayers
using Unitful
using Unitful.DefaultSymbols
using StatsBase
using Plots

# Download temperature and precipitation data
africa_temp = SimpleSDMPredictor(WorldClim, BioClim, 1, left = -25, right = 50, bottom = -35, top = 40)
bio_africa = Worldclim_bioclim(africa_temp, °C)
africa_water = SimpleSDMPredictor(WorldClim, BioClim, 12, left = -25, right = 50, bottom = -35, top = 40)
bio_africa_water = WaterBudget(Worldclim_bioclim(africa_water, mm))

# Find which grid cells are land
active =  Array{Bool, 2}(.!isnan.(bio_africa.array))

plot(africa_temp)

# Set up initial parameters for ecosystem
numSpecies = 1; grid = size(active); req= 10.0mm; individuals=0; area = 64e6km^2; totalK = 1000.0kJ/km^2

# Set up how much water each species consumes
energy_vec = WaterRequirement(fill(req, numSpecies))

# Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.2
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
movement = AlwaysMovement(kernel, Torus())

# Create species list, including their temperature preferences, seed abundance and native status
opts = fill(284.0K, numSpecies)
vars = fill(4.0K, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
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
times = 10years; timestep = 1month; record_interval = 1month; repeats = 1
lensim = length(0years:record_interval:times)
abuns = zeros(Int64, numSpecies, prod(grid), lensim)
@time simulate_record!(abuns, eco, times, record_interval, timestep);

# Reshape abundances for plotting
abuns = reshape(abuns[1, :, :, 1], grid[1], grid[2], lensim)

# Create a gif (warning, slow!)
anim = @animate for i in 1:lensim
    africa_abun = Float64.(abuns[:, :, i])
    africa_abun[.!(active)] .= NaN
    heatmap(africa_abun, clim = (0, maximum(abuns)), 
    background_color = :lightblue, background_color_outside=:white, 
    grid = false, color = cgrad(:algae, scale = :exp))
end
gif(anim, "examples/Africa.gif", fps = 15)

# Plot start and end abundances, next to temperature and rainfall
africa_startabun = Float64.(abuns[:, :, 1])
africa_startabun[.!(active)] .= NaN
africa_endabun = Float64.(abuns[:, :, end])
africa_endabun[.!(active)] .= NaN
heatmap(africa_startabun, clim = (0, maximum(abuns)), 
    background_color = :lightblue, background_color_outside=:white, 
    grid = false, color = cgrad(:algae, scale = :exp), 
    layout = (@layout [a b; c d]))
heatmap!(africa_endabun, clim = (0, maximum(abuns)), 
    background_color = :lightblue, background_color_outside=:white, 
    grid = false, color = cgrad(:algae, scale = :exp), 
    subplot = 2)
plot!(africa_temp, grid = false, subplot = 3)
plot!(africa_water, grid = false, subplot = 4)
