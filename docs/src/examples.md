# Examples

Here we describe several properties that the model ecosystem should be able to recreate, namely:

- Species are more abundant when more resources are available to them.
- Species' abundances scale with area and are invariant to grid size.
- Species with larger average dispersal distances can move further and faster across the landscape.
- Species have a competitive advantage when their niche preference is close to that of the climate.
- Specialist species with a narrow niche width have a competitive advantage over generalists with a broad niche width, for the same niche preference.
- Large numbers of species are sustained over large areas.

What follows are several examples of this operating in practice.

```julia
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Diversity
using JLD2
using OnlineStats
using Plots
using Distributions
using Diversity
plotlyjs()
```

## 1. Different niche preferences and widths

Here we compare the abundances of species with different temperature preferences and tolerances. On a small patch,we explored the abundance of species given different niche preferences, with all other parameters kept equal. We found that species with niche preference nearer to the 25°C optimum were more abundant, when all species were given the same niche widths. Additionally, when all species had a preference for the 25°C climate and a range of niche widths, those with broader niche widths (generalists) were less abundant than species with narrow (specialists). If the temperature in the ecosystem was then increased by 1°C, those with the narrowest niches went extinct, and the generalists became more abundant, with a preference for those with a niche width of around 1°C as we would expect.

``` julia
numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2);
individuals = 100_000_000; area = 100.0*km^2;
totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
bud = BudgetCollection2(abenv1.budget, abenv2.budget)
abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

vars = fill(2.0, numSpecies) .* K
opts = 298.0K .+ vars .* range(-3, stop = 3, length = numSpecies)

av_dist = fill(2.4, numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death = 0.15/ year
birth = death
l = 1.0
s = 0.1
boost = 1.0

size_mean = 1.0m^2
# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

energy_vec = ReqCollection2(energy_vec1, energy_vec2)
param = EqualPop(birth, death, l, s , boost)

# Create ecosystem

movement = BirthOnlyMovement(kernel, Torus())

traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
abun = rand(Multinomial(individuals, numSpecies))
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
rel = Gauss{typeof(first(opts))}()
eco = Ecosystem(sppl, abenv, rel)

times = 10years; timestep = 1month
lensim = length(0month:timestep:times)

simulate!(eco, times, timestep)
endabun = eco.abundances.matrix
temps = map(eachindex(opts)) do i
    repeat([opts[i]], endabun[i])
end

temps = vcat(temps...)
edges = collect(292.0:1:304) .* K
h = Hist(edges)
fit!(h, temps)
bar(ustrip.(uconvert.(°C, edges)), h.counts,
grid = false, xlab = "Temperature preference (°C)",
ylab = "Abundance",
guidefontsize = 16, tickfontsize= 16, size = (1200, 1000),
titlefontsize = 16, title = "A", titleloc = :left,
margin = 10.0*Plots.mm, label = "", layout = (@layout[a; b c]))

## DIFFERENT VARS ##

numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2);
individuals = 100_000_000; area = 100.0*km^2;
totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
bud = BudgetCollection2(abenv1.budget, abenv2.budget)
abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

vars = range(0.0001, stop = 5, length = numSpecies) .* K
opts = fill(298.0K, numSpecies)

av_dist = fill(2.4, numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death = 0.15/ year
birth = death
l = 1.0
s = 0.1
boost = 1.0

size_mean = 1.0m^2
# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

energy_vec = ReqCollection2(energy_vec1, energy_vec2)
param = EqualPop(birth, death, l, s , boost)

# Create ecosystem

movement = BirthOnlyMovement(kernel, Torus())

traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
abun = rand(Multinomial(individuals, numSpecies))
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
rel = Gauss{typeof(first(opts))}()
eco = Ecosystem(sppl, abenv, rel)

times = 10years; timestep = 1month
lensim = length(0month:timestep:times)

simulate!(eco, times, timestep)
endabun = eco.abundances.matrix
widths = map(eachindex(vars)) do i
    repeat([vars[i]], endabun[i])
end
widths = vcat(widths...)
edges = collect(0.1:0.2:5) .* K
h = Hist(edges)
fit!(h, widths)
bar!(edges./K, h.counts, grid = false,
xlab = "Niche width (°C)", ylab = "Abundance",
guidefontsize = 16, tickfontsize= 16, titlefontsize = 16,
margin = 10.0*Plots.mm, left_margin = 20.0 * Plots.mm, label = "",
subplot = 2,
title = "B", titleloc = :left, ylim = (0, 32_000))

## DIFFERENT VARS MISMATCH ##
numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2);
individuals = 100_000_000; area = 100.0*km^2;
totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

abenv1 = simplehabitatAE(299.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(299.0K, grd, totalK[2], area)
bud = BudgetCollection2(abenv1.budget, abenv2.budget)
abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

vars = range(0.0001, stop = 5, length = numSpecies) .* K
opts = fill(298.0K, numSpecies)

av_dist = fill(2.4, numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death = 0.15/ year
birth = death
l = 1.0
s = 0.1
boost = 1.0

size_mean = 1.0m^2
# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

energy_vec = ReqCollection2(energy_vec1, energy_vec2)
param = EqualPop(birth, death, l, s , boost)

# Create ecosystem

movement = BirthOnlyMovement(kernel, Torus())

traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
abun = rand(Multinomial(individuals, numSpecies))
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
rel = Gauss{typeof(first(opts))}()
eco = Ecosystem(sppl, abenv, rel)

times = 10years; timestep = 1month
lensim = length(0month:timestep:times)

simulate!(eco, times, timestep)
endabun = eco.abundances.matrix
widths = map(eachindex(vars)) do i
    repeat([vars[i]], endabun[i])
end
widths = vcat(widths...)
edges = collect(0.1:0.2:5) .* K
h = Hist(edges)
fit!(h, widths)
bar!(edges./K, h.counts, grid = false,
xlab = "Niche width (°C)", ylab = "",
guidefontsize = 16, tickfontsize= 16, titlefontsize = 16,
margin = 10.0*Plots.mm, label = "", subplot = 3,
title = "C", titleloc = :left, ylim = (0, 32_000))
```

![](Opt_var_panel.svg)
*Abundance of species across 100km² patch ecosystem with 100 species, (A) with a different temperature preferences and a homogeneous climate of 25°C, (B)  with different niche widths and a temperature preference for 25°C and (C) different niche widths with a shifted homogeneous climate of 26°C.*

## 2. Varying resources, grid sizes, areas and number of species

Firstly, we confirmed that abundance depended upon the amount of available resource. Here, we simulated an island ecosystem with two resources, water and sunlight, each on a gradient West to East and South to North, respectively. All species were seeded with the same resource requirements and vital rates. Abundance increased in squares with greater amounts of water and sunlight, with some edge effects. Next, we investigated the relationship between abundance and area size. As expected, ecosystems with greater areas could support more individuals, and these abundances were invariant to the resolution of the grid. We also tested in an ecosystem in which species demographic and dispersal rates and resource requirements varied. Under these circumstances, some species not favoured for the conditions go extinct, but most species survive to the end of the simulation.

``` julia
## MORE ENERGY MORE ABUNDANCE ##
numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2); individuals = 100_000_000; area = 100.0*km^2; totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
bud = BudgetCollection2(abenv1.budget, abenv2.budget)
abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)
gsize = size(abenv.budget.b1.matrix, 1)
sol_range = collect(range(0.0kJ, stop = 4.5e11kJ, length = gsize))
 map(1:gsize) do seq
   abenv.budget.b1.matrix[seq, :] .= sol_range[seq]
 end
abenv.budget.b1.matrix

gsize = size(abenv.budget.b2.matrix, 1)
water_range = collect(range(0.0mm, stop = 192mm, length = gsize))
map(1:gsize) do seq
    abenv.budget.b2.matrix[:, seq] .= water_range[seq]
end
abenv.budget.b2.matrix

vars = fill(2.0, numSpecies) .* K
opts = fill(298.0, numSpecies) .* K
av_dist = fill(2.4, numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death = 0.15/ year
birth = death
l = 1.0
s = 0.1
boost = 1.0

size_mean = 1.0m^2
# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

energy_vec = ReqCollection2(energy_vec1, energy_vec2)
param = EqualPop(birth, death, l, s , boost)

# Create ecosystem

movement = BirthOnlyMovement(kernel, NoBoundary())

traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
abun = rand(Multinomial(individuals, numSpecies))
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
rel = Gauss{typeof(first(opts))}()
eco = Ecosystem(sppl, abenv, rel)

times = 10years; timestep = 1month
lensim = length(0month:timestep:times)

simulate!(eco, times, timestep)
endabun = sum(eco.abundances.matrix, dims = 1)
endabun = reshape(endabun, 10, 10)

heatmap(sol_range./kJ, water_range./mm, endabun, grid = false,
 xlab = "Solar energy", ylab = "Water", size = (1600, 1200),
guidefontsize = 16, tickfontsize= 16, titlefontsize=24,
margin = 10.0*Plots.mm, legendfontsize = 16, label = "", left_margin = 20.0 * Plots.mm,
layout = (@layout [a b; c d]), title = "A", titleloc = :left)



## INVARIANT TO GRID SIZE ##
times = 10years; timestep = 1month
lensim = length(0month:timestep:times)
endabuns = zeros(Int64, 4)
grids = [1,2,5,10]
for i in eachindex(grids)
    numSpecies = 100; grd = (grids[i],grids[i]); req=(450000.0kJ/m^2, 192.0nm/m^2);
    individuals = 100_000_000; area = 100.0*km^2;
    totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

    abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
    abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
    bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    vars = fill(2.0, numSpecies) .* K
    opts = fill(298.0, numSpecies) .* K
    av_dist = fill(2.4, numSpecies) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)

    death = 0.15/ year
    birth = death
    l = 1.0
    s = 0.1
    boost = 1.0

    size_mean = 1.0m^2
    # Set up how much energy each species consumes
    energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
    energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

    energy_vec = ReqCollection2(energy_vec1, energy_vec2)
    param = EqualPop(birth, death, l, s , boost)

    # Create ecosystem

    movement = BirthOnlyMovement(kernel, NoBoundary())

    traits = GaussTrait(opts, vars)
    native = fill(true, numSpecies)
    abun = rand(Multinomial(individuals, numSpecies))
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    rel = Gauss{typeof(first(opts))}()
    eco = Ecosystem(sppl, abenv, rel)
    simulate!(eco, times, timestep)
    endabuns[i] = sum(eco.abundances.matrix)
end

bar!(string.(grids), endabuns,
grid = false, xlab = "Number of grid squares",
ylab = "Total abundance",
guidefontsize = 16,tickfontsize= 16, titlefontsize=24,
margin = 10.0*Plots.mm, label = "", left_margin = 20.0 * Plots.mm,
subplot = 2, title = "B", titleloc = :left)

## ABUNDANCE SCALES WITH AREA ##
times = 10years; timestep = 1month
lensim = length(0month:timestep:times)
endabuns = zeros(Int64, 4)
areas = [10.0,20.0,50.0,100.0]
for i in eachindex(areas)
    numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2);
    individuals = 100_000_000; area = areas[i].*km^2;
    totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

    abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
    abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
    bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    vars = fill(2.0, numSpecies) .* K
    opts = fill(298.0, numSpecies) .* K
    av_dist = fill(2.4, numSpecies) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)

    death = 0.15/ year
    birth = death
    l = 1.0
    s = 0.1
    boost = 1.0

    size_mean = 1.0m^2
    # Set up how much energy each species consumes
    energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
    energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

    energy_vec = ReqCollection2(energy_vec1, energy_vec2)
    param = EqualPop(birth, death, l, s , boost)

    # Create ecosystem

    movement = BirthOnlyMovement(kernel, NoBoundary())

    traits = GaussTrait(opts, vars)
    native = fill(true, numSpecies)
    abun = rand(Multinomial(individuals, numSpecies))
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    rel = Gauss{typeof(first(opts))}()
    eco = Ecosystem(sppl, abenv, rel)
    simulate!(eco, times, timestep)
    endabuns[i] = sum(eco.abundances.matrix)
end

bar!(string.(areas), endabuns, grid = false, xlab = "Area (km²)",
ylab = "Total abundance", guidefontsize = 16,
tickfontsize= 16, titlefontsize=24, margin = 10.0*Plots.mm,
label = "", subplot = 3, title = "C", titleloc = :left,
left_margin = 20.0 *Plots.mm)

## Sustain large number of species ##
times = 10years; timestep = 1month
lensim = length(0month:timestep:times)
reps = 10
species = [100, 500, 1_000, 5_000]
SR = zeros(Float64, length(species), reps)
for r in 1:reps
    for i in eachindex(species)
        numSpecies = species[i]; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2);
        individuals = 100_000_000; area = 100.0km^2;
        totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

        abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
        abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
        bud = BudgetCollection2(abenv1.budget, abenv2.budget)
        abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

        vars = rand(Uniform(1.0, 5.0), numSpecies) .* K
        opts = 298.0K .+ vars .* range(-3, stop = 3, length = numSpecies)
        av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
        kernel = GaussianKernel.(av_dist, 10e-10)

        death = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
        birth = death
        l = 1.0
        s = 0.1
        boost = 1.0

        size_mean = rand(Normal(1.0, 0.5), numSpecies) .* m^2
        # Set up how much energy each species consumes
        energy_vec1 = SolarRequirement(abs.(req[1] .* size_mean))
        energy_vec2 = WaterRequirement(abs.(req[2] .* size_mean))

        energy_vec = ReqCollection2(energy_vec1, energy_vec2)
        param = PopGrowth{typeof(unit(birth[1]))}(birth, death, l, s , boost)

        # Create ecosystem

        movement = BirthOnlyMovement(kernel, NoBoundary())

        traits = GaussTrait(opts, vars)
        native = fill(true, numSpecies)
        abun = rand(Multinomial(individuals, numSpecies))
        sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
            movement, param, native)
        rel = Gauss{typeof(first(opts))}()
        eco = Ecosystem(sppl, abenv, rel)
        simulate!(eco, times, timestep)
        SR[i, r] = sum(sum(eco.abundances.matrix, dims = 2) .> 0)
        print(".")
    end
end

meanSR = dropdims(mean(SR, dims = 2), dims = 2)
sdSR = dropdims(std(SR, dims = 2), dims = 2)

bar!(string.(species), meanSR ./species, yerr= sdSR ./ species, grid = false, xlab = "Number of species introduced",
ylab = "% Species survived", guidefontsize = 16,
tickfontsize= 16, titlefontsize=24, margin = 10.0*Plots.mm,
label = "",  title = "D", subplot = 4, titleloc = :left,
left_margin = 20.0 *Plots.mm, ylim = (0, 1))
```

![](Abundance.svg)
*Model testing on island ecosystems. (A) Total abundance of 100 species, with varying resources of water and solar energy across the grid. (B) Total abundance of 100 species, with increasing area size. (C) Total abundance of 100 species, with increasing grid square resolution. (D) Percentage of species survived after 10 years of simulation.*

## 3. Dispersal

Finally, we verify that dispersal is functioning as expected. This figure shows the overall abundance of an island populated with two species at opposite extremes of the ecosystem after ten years of simulation. The species moved faster and further into the unpopulated island centre when they had higher average dispersal distances, though again with some edge effects.

```julia
times = 50years; timestep = 1month
lensim = length(0month:timestep:times)
distances = [0.5, 1.0, 2.0, 4.0]
endabuns = zeros(Int64, 10, 10, length(distances))
for i in eachindex(distances)
    numSpecies = 2; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2);
    individuals = 0; area = 100.0km^2;
    totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

    abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
    abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
    bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    vars = fill(2.0, numSpecies) .* K
    opts = fill(298.0, numSpecies) .* K
    av_dist = fill(distances[i], numSpecies) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)

    death = 0.15/ year
    birth = death
    l = 1.0
    s = 0.1
    boost = 1.0

    size_mean = 1.0m^2
    # Set up how much energy each species consumes
    energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
    energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

    energy_vec = ReqCollection2(energy_vec1, energy_vec2)
    param = EqualPop(birth, death, l, s , boost)

    # Create ecosystem

    movement = BirthOnlyMovement(kernel, NoBoundary())

    traits = GaussTrait(opts, vars)
    native = fill(true, numSpecies)
    abun = rand(Multinomial(individuals, numSpecies))
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    rel = Gauss{typeof(first(opts))}()
    eco = Ecosystem(sppl, abenv, rel)
    eco.abundances.grid[1, :, 1] .= 100.0
    eco.abundances.grid[2, :, 10] .= 100.0
    simulate!(eco, times, timestep)
    endabuns[:, :, i] = sum(eco.abundances.matrix, dims = 1)
end

heatmap(grid = false, xlab = "Distance (km)",
ylab = "Distance (km)", size = (1200, 800),
guidefontsize = 12,tickfontsize= 12, titlefontsize=18,
margin = 10.0*Plots.mm, legendfontsize = 12, label = "",
layout = (@layout [a b; c d]), link = :both)
titles = ["A", "B", "C", "D"]
for i in 1:4
    m = distances[i]
    display(heatmap!(1:10,1:10, endabuns[:, :, i],
    grid = false, xlab = "Distance (km)", ylab = "Distance (km)",
    guidefontsize = 16, tickfontsize= 16,
    titlefontsize=24, title = titles[i], margin = 10.0*Plots.mm,
    label = "", subplot = i, titleloc = :left,
    clim = (0, 1.5e4), link = :both))
end
```

![](DispersalSD.svg)
*Total abundance of two species in island ecosystems after 10 years of simulation, with species populated at opposite sides of the island. Those with higher dispersal distances moved further away from their starting populations at a faster rate. (A) Mean dispersal distance of 0.5km, (B) mean dispersal distance of 1km, (C) Mean dispersal distance of 2km, Mean dispersal distance of 4km.*
