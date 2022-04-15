using Test
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using OnlineStats
using Distributions

## DIFFERENT OPTS ##

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
longevity = 1.0
survival = 0.1
boost = 1.0

size_mean = 1.0m^2
# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

energy_vec = ReqCollection2(energy_vec1, energy_vec2)
param = EqualPop(birth, death, longevity, survival , boost)

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
# Check abundances are greater at temperature optimum than edges
@test h.counts[edges[1:end-1] .== 298.0K][1] > h.counts[edges[1:end-1] .== 292.0K][1]
@test h.counts[edges[1:end-1] .== 298.0K][1] > h.counts[edges[1:end-1] .== 303.0K][1]


## DIFFERENT VARS ##

numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2);
individuals = 100_000_000; area = 100.0*km^2;
totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
bud = BudgetCollection2(abenv1.budget, abenv2.budget)
abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

collect(range(0.0001K, stop = 5.0K, length = numSpecies))
opts = fill(298.0K, numSpecies)

av_dist = fill(2.4, numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death = 0.15/ year
birth = death
longevity = 1.0
survival = 0.1
boost = 1.0

size_mean = 1.0m^2
# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

energy_vec = ReqCollection2(energy_vec1, energy_vec2)
param = EqualPop(birth, death, longevity, survival , boost)

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

# Check abundances are greater at smaller niche widths
@test h.counts[edges[1:end-1] .== 0.1K][1] > h.counts[edges[1:end-1] .== 4.7K][1]

## DIFFERENT VARS MISMATCH ##
numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2);
individuals = 100_000_000; area = 100.0*km^2;
totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

abenv1 = simplehabitatAE(299.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(299.0K, grd, totalK[2], area)
bud = BudgetCollection2(abenv1.budget, abenv2.budget)
abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

vars = collect(range(0.0001, stop = 5, length = numSpecies) .* K)
opts = fill(298.0K, numSpecies)

av_dist = fill(2.4, numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death = 0.15/ year
birth = death
longevity = 1.0
survival = 0.1
boost = 1.0

size_mean = 1.0m^2
# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

energy_vec = ReqCollection2(energy_vec1, energy_vec2)
param = EqualPop(birth, death, longevity, survival, boost)

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

# Check abundances are greater at smaller niche widths (with temperature shifted by 1.0K)
@test h.counts[edges[1:end-1] .== 1.1K][1] > h.counts[edges[1:end-1] .== 4.7K][1]
@test h.counts[edges[1:end-1] .== 1.1K][1] > h.counts[edges[1:end-1] .== 0.1K][1]

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
longevity = 1.0
survival = 0.1
boost = 1.0

size_mean = 1.0m^2
# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

energy_vec = ReqCollection2(energy_vec1, energy_vec2)
param = EqualPop(birth, death, longevity, survival, boost)

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

# Check abundances greater in areas of greater resource
@test all(endabun[1, :] .< endabun[end, :])
@test all(endabun[:, 1] .< endabun[:, end])

## INVARIANT TO GRID SIZE ##
times = 10years; timestep = 1month
lensim = length(0month:timestep:times)
endabuns = zeros(Int64, 4)
grids = [1,2,5,10]
for i in eachindex(grids)
    local numSpecies = 100; local grd = (grids[i],grids[i]); local req=(450000.0kJ/m^2, 192.0nm/m^2);
    local individuals = 100_000_000; local area = 100.0*km^2;
    local totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

    local abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
    local abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
    local bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    local abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    local vars = fill(2.0, numSpecies) .* K
    local opts = fill(298.0, numSpecies) .* K
    local av_dist = fill(2.4, numSpecies) .* km
    local kernel = GaussianKernel.(av_dist, 10e-10)

    local death = 0.15/ year
    local birth = death
    local longevity = 1.0
    local survival = 0.1
    local boost = 1.0

    local size_mean = 1.0m^2
    # Set up how much energy each species consumes
    local energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
    local energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

    local energy_vec = ReqCollection2(energy_vec1, energy_vec2)
    local param = EqualPop(birth, death, longevity, survival, boost)

    # Create ecosystem

    local movement = BirthOnlyMovement(kernel, NoBoundary())

    local traits = GaussTrait(opts, vars)
    local native = fill(true, numSpecies)
    local abun = rand(Multinomial(individuals, numSpecies))
    local sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    local rel = Gauss{typeof(first(opts))}()
    local eco = Ecosystem(sppl, abenv, rel)
    simulate!(eco, times, timestep)
    endabuns[i] = sum(eco.abundances.matrix)
end

# Check answers are similar despite different grid sizes
@test all(isapprox.(endabuns[1], endabuns, rtol = 0.1))

## ABUNDANCE SCALES WITH AREA ##
times = 10years; timestep = 1month
lensim = length(0month:timestep:times)
endabuns = zeros(Int64, 4)
areas = [10.0,20.0,50.0,100.0]
for i in eachindex(areas)
    local numSpecies = 100; local grd = (10,10); local req=(450000.0kJ/m^2, 192.0nm/m^2);
    local individuals = 100_000_000; local area = areas[i].*km^2;
    local totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

    local abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
    local abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
    local bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    local abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    local vars = fill(2.0, numSpecies) .* K
    local opts = fill(298.0, numSpecies) .* K
    local av_dist = fill(2.4, numSpecies) .* km
    local kernel = GaussianKernel.(av_dist, 10e-10)

    local death = 0.15/ year
    local birth = death
    local longevity = 1.0
    local survival = 0.1
    local boost = 1.0

    local size_mean = 1.0m^2
    # Set up how much energy each species consumes
    local energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
    local energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

    local energy_vec = ReqCollection2(energy_vec1, energy_vec2)
    local param = EqualPop(birth, death, longevity, survival , boost)

    # Create ecosystem

    local movement = BirthOnlyMovement(kernel, NoBoundary())

    local traits = GaussTrait(opts, vars)
    local native = fill(true, numSpecies)
    local abun = rand(Multinomial(individuals, numSpecies))
    local sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    local rel = Gauss{typeof(first(opts))}()
    local eco = Ecosystem(sppl, abenv, rel)
    simulate!(eco, times, timestep)
    endabuns[i] = sum(eco.abundances.matrix)
end
# Test larger areas can support more individuals
@test endabuns[1] < endabuns[2] < endabuns[3] < endabuns[4]

## Sustain large number of species ##
times = 10years; timestep = 1month
lensim = length(0month:timestep:times)
species = [50, 100, 200, 500]
SR = zeros(Float64, length(species))
for i in eachindex(species)
    local numSpecies = species[i]; local grd = (10,10); local req=(450000.0kJ/m^2, 192.0nm/m^2);
    local individuals = 100_000_000; local area = 100.0km^2;
    local totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

    local abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
    local abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
    local bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    local abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    local vars = rand(Uniform(1.0, 5.0), numSpecies) .* K
    local opts = 298.0K .+ vars .* range(-3, stop = 3, length = numSpecies)
    local av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
    local kernel = GaussianKernel.(av_dist, 10e-10)

    local death = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
    local birth = death
    local longevity = 1.0
    local survival = 0.1
    local boost = 1.0

    local size_mean = rand(Normal(1.0, 0.05), numSpecies) .* m^2
    # Set up how much energy each species consumes
    local energy_vec1 = SolarRequirement(abs.(req[1] .* size_mean))
    local energy_vec2 = WaterRequirement(abs.(req[2] .* size_mean))

    local energy_vec = ReqCollection2(energy_vec1, energy_vec2)
    local param = PopGrowth{typeof(unit(birth[1]))}(birth, death, longevity, survival, boost)

    # Create ecosystem

    local movement = BirthOnlyMovement(kernel, NoBoundary())

    local traits = GaussTrait(opts, vars)
    local native = fill(true, numSpecies)
    local abun = rand(Multinomial(individuals, numSpecies))
    local sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    local rel = Gauss{typeof(first(opts))}()
    local eco = Ecosystem(sppl, abenv, rel)
    simulate!(eco, times, timestep)
    SR[i] = sum(sum(eco.abundances.matrix, dims = 2) .> 0)
end

# Test most species can survive
@test all(SR./species .>= 0.9)

## DISPERSAL ##
times = 50years; timestep = 1month
lensim = length(0month:timestep:times)
distances = [0.5, 1.0, 2.0, 4.0]
endabuns = zeros(Int64, 10, 10, length(distances))
for i in eachindex(distances)
    local numSpecies = 2; local grd = (10,10); local req=(450000.0kJ/m^2, 192.0nm/m^2);
    local individuals = 0; local area = 100.0km^2;
    local totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

    local abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
    local abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
    local bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    local abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    local vars = fill(2.0, numSpecies) .* K
    local opts = fill(298.0, numSpecies) .* K
    local av_dist = fill(distances[i], numSpecies) .* km
    local kernel = GaussianKernel.(av_dist, 10e-10)

    local death = 0.15/ year
    local birth = death
    local longevity = 1.0
    local survival = 0.1
    local boost = 1.0

    local size_mean = 1.0m^2
    # Set up how much energy each species consumes
    local energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
    local energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))

    local energy_vec = ReqCollection2(energy_vec1, energy_vec2)
    local param = EqualPop(birth, death, longevity, survival, boost)

    # Create ecosystem

    local movement = BirthOnlyMovement(kernel, NoBoundary())

    local traits = GaussTrait(opts, vars)
    local native = fill(true, numSpecies)
    local abun = rand(Multinomial(individuals, numSpecies))
    local sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
        movement, param, native)
    local rel = Gauss{typeof(first(opts))}()
    local eco = Ecosystem(sppl, abenv, rel)
    eco.abundances.grid[1, :, 1] .= 100.0
    eco.abundances.grid[2, :, 10] .= 100.0
    simulate!(eco, times, timestep)
    endabuns[:, :, i] = sum(eco.abundances.matrix, dims = 1)
end

# Test that species disperse further into the middle of the grid with larger dispersal distances
@test sum(endabuns[:, 1, 1] .> endabuns[:, 1, 2] .> endabuns[:, 1, 3] .> endabuns[:, 1, 4]) >= 9
@test sum(endabuns[:, end, 1] .> endabuns[:, end, 2] .> endabuns[:, end, 3] .> endabuns[:, end, 4]) >= 9
@test sum(endabuns[:, 5, 1] .< endabuns[:, 5, 2] .< endabuns[:, 5, 3] .< endabuns[:, 5, 4]) >= 9
