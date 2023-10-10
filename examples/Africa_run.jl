Sys.total_memory()/1e9 >= 100 || error("You do not have enough memory to run these examples!")

#### SINGLE SPECIES ####

using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Distances
using StatsBase
using Plots
using LinearAlgebra
file = "Africa.tif"
africa = readfile(file, -25.0°, 50.0°, -35.0°, 40.0°)
active = Array{Bool, 2}(fill(0, size(africa)))

xs = 1:size(africa, 1)
ys = 1:size(africa, 2)
radius = 50
for x in xs
    for y in ys
        if norm((x-radius,y-radius)) < radius
            active[x, y] = 1
        end
    end
end

heatmap(active)

# Set up initial parameters for ecosystem
numSpecies = 1; grid = size(africa); req= 10.0kJ; individuals=0; area = 64e6km^2; totalK = 1000.0kJ/km^2

# Set up how much energy each species consumes
energy_vec = SolarRequirement(fill(req, numSpecies))


# Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.0
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
movement = AlwaysMovement(kernel, Torus())


# Create species list, including their temperature preferences, seed abundance and native status
opts = fill(274.0K, numSpecies)
vars = fill(0.5K, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
sppl.params.birth

# Create abiotic environment - even grid of one temperature
abenv = simplehabitatAE(274.0K, grid, totalK, area, active)


# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

#Create ecosystem
eco = Ecosystem(sppl, abenv, rel)
rand_start = rand(findall(active), 1)[1]
eco.abundances.grid[1, rand_start[1], rand_start[2]] = 100

# Simulation Parameters
times = 100years; timestep = 1month; record_interval = 1month; repeats = 1
lensim = length(0years:record_interval:times)
abuns = zeros(Int64, numSpecies, prod(grid), lensim)
@time simulate_record!(abuns, eco, times, record_interval, timestep);

abuns = reshape(abuns[1, :, :, 1], grid[1], grid[2], lensim)

anim = @animate for i in 1:lensim
    africa_abun = Float64.(abuns[:, :, i])
    africa_abun[.!(active)] .= NaN
    heatmap(africa_abun, clim = (0, 700_000), background_color = :lightblue, background_color_outside=:white, grid = false, color = cgrad(:algae, scale = :exp), aspect_ratio = 1)
end
gif(anim, "Africa.gif", fps = 30)

#### SPECIALIST VERSUS GENERALIST ####

specialist_vars = [0.5K, 1.0K, 5.0K, 10.0K, 25.0K, 50.0K]
velocity = zeros(typeof(1.0km/year), length(specialist_vars))
rand_start = rand(findall(active), 1)[1]
for i in eachindex(specialist_vars)
    # Set up initial parameters for ecosystem
    numSpecies = 2; grid = size(africa); req= 10.0kJ; individuals=0; area = 64e6km^2; totalK = 1000.0kJ/km^2

    # Set up how much energy each species consumes
    energy_vec = SolarRequirement(fill(req, numSpecies))


    # Set rates for birth and death
    birth = 0.6/year
    death = 0.6/year
    longevity = 1.0
    survival = 0.1
    boost = 1.0
    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
    movement = AlwaysMovement(kernel, Torus())


    # Create species list, including their temperature preferences, seed abundance and native status
    opts = fill(274.0K, numSpecies)
    vars = [50.0K, specialist_vars[i]]
    traits = GaussTrait(opts, vars)
    native = fill(true, numSpecies)
    # abun = rand(Multinomial(individuals, numSpecies))
    abun = fill(div(individuals, numSpecies), numSpecies)
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    sppl.params.birth

    # Create abiotic environment - even grid of one temperature
    abenv = simplehabitatAE(274.0K, grid, totalK, area, active)

    # Set relationship between species and environment (gaussian)
    rel = Gauss{typeof(1.0K)}()

    #Create ecosystem
    eco = Ecosystem(sppl, abenv, rel)
    eco.abundances.grid[1, rand_start[1], rand_start[2]] = 100

    # Simulation Parameters
    burnin = 100years; times = 100years; timestep = 1month; record_interval = 1month; repeats = 1
    lensim = length(0years:record_interval:times)
    simulate!(eco, burnin,timestep)
    eco.abundances.grid[2, rand_start[1], rand_start[2]] = 100
    abuns = zeros(Int64, numSpecies, prod(grid), lensim)
    @time simulate_record!(abuns, eco, times, record_interval, timestep);

    abuns = reshape(abuns[:, :, :, 1], numSpecies, grid[1], grid[2], lensim)
    origin = [rand_start[1], rand_start[2]]
    dest = findall(abuns[2, :, :, end] .> 0)
    dists = [euclidean(origin, [dest[i][1], dest[i][2]]) for i in length(dest)] .* getgridsize(eco)
    velocity[i] = mean(dists) / 100years
    # inst_velocity = map(1:lensim) do t
    #     dest = findall(abuns[2, :, :, t] .> 0)
    #     dists = [euclidean(origin, [dest[i][1], dest[i][2]]) for i in length(dest)] .* getgridsize(eco)
    #     return maximum(dists)/month
    # end
    #velocity[i] = mean(inst_velocity)
    # anim = @animate for i in 1:lensim
    #     africa_abun1 = Float64.(abuns[1, :, :, i])
    #     africa_abun1[.!(active)] .= NaN
    #     africa_abun2 = Float64.(abuns[2, :, :, i])
    #     africa_abun2[.!(active)] .= NaN
    #     heatmap(africa_abun1, clim = (0, 700_000),
    #     background_color = :lightblue,
    #     background_color_outside=:white,
    #     grid = false, color = cgrad(:algae, scale = :exp),
    #     aspect_ratio = 1, layout = (@layout [a b]), subplot = 1)
    #     heatmap!(africa_abun2, clim = (0, 700_000),
    #     background_color = :lightblue,
    #     background_color_outside=:white,
    #     grid = false, color = cgrad(:algae, scale = :exp),
    #     aspect_ratio = 1, subplot = 2)
    # end
    # gif(anim, "examples/Biodiversity/Africa_$i.gif", fps = 30)
end

plot(ustrip.(abs.(specialist_vars .- 50.0K)), ustrip.(velocity),
    xlab = "Selective advantage", ylab = "Average invasion speed (km/year)",
    label = "", grid = false)
Plots.pdf("InvasionCircle.pdf")


#### SPECIALIST VERSUS MANY GENERALISTS ####
using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using JLD2
using Printf
file = "Africa.tif"
africa = readfile(file, -25.0°, 50.0°, -35.0°, 40.0°)
active = Array{Bool, 2}(fill(0, size(africa)))

xs = 1:size(africa, 1)
ys = 1:size(africa, 2)
radius = 50
for x in xs
    for y in ys
        if norm((x-radius,y-radius)) < radius
            active[x, y] = 1
        end
    end
end

heatmap(active)

# Set up initial parameters for ecosystem
numSpecies = 50_000; grid = size(africa); req= 10.0kJ; individuals=3*10^8; area = 64e6km^2; totalK = 1000.0kJ/km^2

# Set up how much energy each species consumes
energy_vec = SolarRequirement(fill(req, numSpecies))


# Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.1
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
movement = AlwaysMovement(kernel, Torus())


# Create species list, including their temperature preferences, seed abundance and native status
opts = fill(274.0K, numSpecies)
vars = fill(50.0K, numSpecies)
vars[50_000] = 0.5K
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
sppl.params.birth

# Create abiotic environment - even grid of one temperature
abenv = simplehabitatAE(274.0K, grid, totalK, area, active)


# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

#Create ecosystem
eco = Ecosystem(sppl, abenv, rel)
eco.abundances.matrix[50_000, :] .= 0


# Simulation Parameters
burnin = 100years; times = 100years; timestep = 1month; record_interval = 12months;
lensim = length(0years:record_interval:times)
@time simulate!(eco, burnin, timestep)
rand_start = rand(findall(active), 1)[1]
eco.abundances.grid[50_000, rand_start[1], rand_start[2]] = 100
@time simulate!(eco, times, timestep, record_interval, "/home/claireh/sdc/Africa/specialist2", "Africa_run");

#### 50,000 SPECIES COEXISTING #####

using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using JLD2
using Printf

file = "Africa.tif"
africa = readfile(file, -25.0°, 50.0°, -35.0°, 40.0°)
active = Array{Bool, 2}(fill(0, size(africa)))

xs = 1:size(africa, 1)
ys = 1:size(africa, 2)
radius = 50
for x in xs
    for y in ys
        if norm((x-radius,y-radius)) < radius
            active[x, y] = 1
        end
    end
end

heatmap(active)

# Set up initial parameters for ecosystem
numSpecies = 50_000; grid = size(africa); req= 10.0kJ; individuals=3*10^8; area = 64e6km^2; totalK = 1000.0kJ/km^2

# Set up how much energy each species consumes
energy_vec = SolarRequirement(fill(req, numSpecies))


# Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.1
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
movement = AlwaysMovement(kernel, Torus())


# Create species list, including their temperature preferences, seed abundance and native status
opts = fill(274.0K, numSpecies)
vars = fill(50.0K, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
sppl.params.birth

# Create abiotic environment - even grid of one temperature
abenv = simplehabitatAE(274.0K, grid, totalK, area, active)


# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

#Create ecosystem
eco = Ecosystem(sppl, abenv, rel)

# Simulation Parameters
burnin = 10years; times = 100years; timestep = 1month; record_interval = 12months;
lensim = length(0years:record_interval:times)
@time simulate!(eco, burnin, timestep)
@time simulate!(eco, times, timestep, record_interval, "/home/claireh/sdc/Africa/coexist2", "Africa_run_coexist");

using JLD2
using Plots
using Diversity
@load "examples/Africa_run_coexist100.jld2" abun
meta = Metacommunity(abun)
div = norm_sub_alpha(meta, 0)
sumabuns = reshape(div[!, :diversity], 100, 100)

heatmap(sumabuns,
    grid = false, color = :algae,
    aspect_ratio = 1, layout = (@layout [a b; c d]),
    clim = (0, 50_000), margin = 0.5 * Plots.mm,
    title = "A", titleloc = :left)

@load "examples/Africa_run50.jld2" abun
meta = Metacommunity(abun)
div = norm_sub_alpha(meta, 0)
sumabuns = reshape(div[!, :diversity], 100, 100)
heatmap!(sumabuns,
    background_color = :lightblue,
    background_color_outside=:white,
    grid = false, color = :algae,
    aspect_ratio = 1, subplot = 2,
    clim = (0, 50_000), right_margin = 2.0 * Plots.mm,
    title = "B", titleloc = :left)

@load "examples/Africa_run100.jld2" abun
meta = Metacommunity(abun)
div = norm_sub_alpha(meta, 0)
sumabuns = reshape(div[!, :diversity], 100, 100)
heatmap!(sumabuns,
    background_color = :lightblue,
    background_color_outside=:white,
    grid = false, color = :algae,
    aspect_ratio = 1, subplot = 3,
    clim = (0, 50_000), right_margin = 2.0 * Plots.mm,
    title = "C", titleloc = :left)

using Diversity.Ecology
@load "examples/Africa_run50.jld2" abun
meta = Metacommunity(abun)
diver = shannon(meta)
sumabuns = reshape(diver[!, :diversity], 100, 100)
heatmap!(sumabuns,
    background_color = :lightblue,
    background_color_outside=:white,
    grid = false, color = :algae,
    aspect_ratio = 1, subplot = 4,
        right_margin = 2.0 * Plots.mm,
    title = "D", titleloc = :left, clim = (0, 10))
Plots.pdf("examples/Africa.pdf")
