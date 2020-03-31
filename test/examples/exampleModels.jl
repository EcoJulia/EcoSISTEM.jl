include("CreateEco.jl")
include("RunSim.jl")
include("Scenarios.jl")
using Diversity
using JLD
using OnlineStats
using Plots
plotlyjs()

cachefolder = joinpath(dirname(pathof(Simulation)), "../test/examples/")

## DIFFERENT TEMPERATURE OPTIMUMS ##
# 100 species with a range of different niche preferences for temperature. Check those closer to the environmental temperature have a higher abundance.

numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2); individuals = 100_000_000; area = 100.0*km^2; totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

scenario = [SimpleScenario(TempIncrease!, 0.0K/10year)]

divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson, mean_abun, geom_mean_abun, norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta, norm_meta_rho, raw_meta_rho, meta_gamma]
q = 1.0
scenario_names = ["DiffOpts"]

simDict = Dict("times" => 10years, "burnin" => 0year, "interval" => 1month, "timestep" => 1month, "scenarios" => scenario, "divfuns" => divfuns, "q" => q, "reps" => 1, "scenario_names" => scenario_names, "cacheInterval" => 1years)
lensim = length(0month:simDict["interval"]:simDict["times"])

abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
bud = BudgetCollection2(abenv1.budget, abenv2.budget)
abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

vars = fill(2.0, numSpecies) .* K
opts = 298.0K .+ vars .* range(-3, stop = 3, length = numSpecies)

av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death_rates = fill(0.15, numSpecies) ./ year
birth_rates = death_rates

paramDict = Dict("numSpecies" => numSpecies, "numInvasive" => 0, "numIndiv" => individuals, "reqs" => req, "opts" => opts, "vars" => vars, "birth" => birth_rates, "death" => death_rates, "s" => 5e-2, "boost" => 1.0, "kernel" => kernel, "totalK" => totalK, "bound" => Torus())

diver = zeros(length(divfuns), lensim, length(scenario))
runsim!(diver, abenv, paramDict, simDict, 1, cachefolder)

endabun = mapslices(sum, JLD.load(joinpath(cachefolder, "cache/DiffOpts10001.jld"), "abun"), dims = 2)[:,1]
temps = map(eachindex(opts)) do i
    repeat([opts[i]], endabun[i])
end

temps = vcat(temps...)
edges = collect(292.0:1:304) .* K
h = Hist(edges)
fit!(h, temps)
bar(ustrip.(uconvert.(°C, edges)), h.counts, grid = false, xlab = "Temperature preference (°C)", ylab = "Abundance", size = (1200, 800), guidefontsize = 22,tickfontsize= 20, titlefontsize=12, margin = 2.0*Plots.mm, legendfontsize = 12, label = "")

## DIFFERENT NICHE WIDTHS ##
# 100 species with a range of different niche widths (variances) for temperature. Check those with narrower niche widths have a higher abundance.

numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2); individuals = 100_000_000; area = 100.0*km^2; totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

scenario = [SimpleScenario(TempIncrease!, 0.0K/10year)]

divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson, mean_abun, geom_mean_abun, norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta, norm_meta_rho, raw_meta_rho, meta_gamma]
q = 1.0
scenario_names = ["DiffVars"]

simDict = Dict("times" => 10years, "burnin" => 0year, "interval" => 1month, "timestep" => 1month, "scenarios" => scenario, "divfuns" => divfuns, "q" => q, "reps" => 1, "scenario_names" => scenario_names, "cacheInterval" => 1years)
lensim = length(0month:simDict["interval"]:simDict["times"])

abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
bud = BudgetCollection2(abenv1.budget, abenv2.budget)
abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

vars = range(0.0001, stop = 5, length = numSpecies) .* K
opts = fill(298.0K, numSpecies)

av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
birth_rates = death_rates

paramDict = Dict("numSpecies" => numSpecies, "numInvasive" => 0, "numIndiv" => individuals, "reqs" => req, "opts" => opts, "vars" => vars, "birth" => birth_rates, "death" => death_rates, "s" => 5e-2, "boost" => 1.0, "kernel" => kernel, "totalK" => totalK, "bound" => Torus())

diver = zeros(length(divfuns), lensim, length(scenario))
runsim!(diver, abenv, paramDict, simDict, 1, cachefolder)

endabun = mapslices(sum, JLD.load(joinpath(cachefolder, "cache/DiffVars10001.jld"), "abun"), dims = 2)[:,1]
widths = map(eachindex(vars)) do i
    repeat([vars[i]], endabun[i])
end
widths = vcat(widths...)
edges = collect(0.1:0.2:5) .* K
h = Hist(edges)
fit!(h, widths)
bar(edges./K, h.counts, grid = false, xlab = "Niche width (°C)", ylab = "Abundance", size = (1200, 800), guidefontsize = 22,tickfontsize= 20, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, label = "")

## DIFFERENT NICHE WIDTHS MISMATCH ##
# Increase temperature of environment by 1 degree C but keep all else the same. Check that there is a corresponding shift in abundance by 1 degree.

numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2); individuals = 100_000_000; area = 100.0*km^2; totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

scenario = [SimpleScenario(TempIncrease!, 0.0K/10year)]

divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson, mean_abun, geom_mean_abun, norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta, norm_meta_rho, raw_meta_rho, meta_gamma]
q = 1.0
scenario_names = ["DiffVarsMismatch"]

simDict = Dict("times" => 10years, "burnin" => 0year, "interval" => 1month, "timestep" => 1month, "scenarios" => scenario, "divfuns" => divfuns, "q" => q, "reps" => 1, "scenario_names" => scenario_names, "cacheInterval" => 1years)
lensim = length(0month:simDict["interval"]:simDict["times"])

abenv1 = simplehabitatAE(299.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(299.0K, grd, totalK[2], area)
bud = BudgetCollection2(abenv1.budget, abenv2.budget)
abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

vars = range(0.0001, stop = 5, length = numSpecies) .* K
opts = fill(298.0K, numSpecies)

av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
birth_rates = death_rates

paramDict = Dict("numSpecies" => numSpecies, "numInvasive" => 0, "numIndiv" => individuals, "reqs" => req, "opts" => opts, "vars" => vars, "birth" => birth_rates, "death" => death_rates, "s" => 5e-2, "boost" => 1.0, "kernel" => kernel, "totalK" => totalK, "bound" => Torus())

diver = zeros(length(divfuns), lensim, length(scenario))
runsim!(diver, abenv, paramDict, simDict, 1, cachefolder)

endabun = mapslices(sum, JLD.load(joinpath(cachefolder, "cache/DiffVarsMismatch10001.jld"), "abun"), dims = 2)[:,1]
widths = map(eachindex(vars)) do i
    repeat([vars[i]], endabun[i])
end
widths = vcat(widths...)
edges = collect(0.1:0.2:5) .* K
h = Hist(edges)
fit!(h, widths)
bar(edges./K, h.counts, grid = false, xlab = "Niche width (°C)", ylab = "Abundance", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, label = "", ylim = (0, 1e7))

## ABUNDANCE SCALES WITH ENERGY ##
# Check that ecosystems with greater available sunlight and water can support a greater number of individuals.

numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2); individuals = 100_000_000; area = 100.0*km^2; totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

scenario = [SimpleScenario(TempIncrease!, 0.0K/10year)]

divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson, mean_abun, geom_mean_abun, norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta, norm_meta_rho, raw_meta_rho, meta_gamma]
q = 1.0
scenario_names = ["Energy"]

simDict = Dict("times" => 10years, "burnin" => 0year, "interval" => 1month, "timestep" => 1month, "scenarios" => scenario, "divfuns" => divfuns, "q" => q, "reps" => 1, "scenario_names" => scenario_names, "cacheInterval" => 1years)
lensim = length(0month:simDict["interval"]:simDict["times"])

abenv1 = simplehabitatAE(299.0K, grd, totalK[1], area)
abenv2 = simplehabitatAE(299.0K, grd, totalK[2], area)
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

av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
kernel = GaussianKernel.(av_dist, 10e-10)

death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
birth_rates = death_rates

paramDict = Dict("numSpecies" => numSpecies, "numInvasive" => 0, "numIndiv" => individuals, "reqs" => req, "opts" => opts, "vars" => vars, "birth" => birth_rates, "death" => death_rates, "s" => 5e-2, "boost" => 1.0, "kernel" => kernel, "totalK" => totalK, "bound" => NoBoundary())

diver = zeros(length(divfuns), lensim, length(scenario))
runsim!(diver, abenv, paramDict, simDict, 1, cachefolder, false)

endabun = mapslices(sum, JLD.load(joinpath(cachefolder, "cache/Energy10001.jld"), "abun"), dims = 1)[1,:]
endabun = reshape(endabun, 10, 10)

heatmap(sol_range./kJ, water_range./mm, endabun, grid = false, xlab = "Solar energy", ylab = "Water", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, label = "")

## INVARIANT TO GRID SIZE ##
# Check that same ecosystems with same total area support the same number of species at different grid sizes.

for i in [1,2,5,10]
    numSpecies = 100; grd = (i,i); req=(450000.0kJ/m^2, 192.0nm/m^2); individuals = 100_000_000; area = 100.0*km^2; totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

    scenario = [SimpleScenario(TempIncrease!, 0.0K/10year)]

    divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson, mean_abun, geom_mean_abun, norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta, norm_meta_rho, raw_meta_rho, meta_gamma]
    q = 1.0
    scenario_names = ["GridSize$i"]

    simDict = Dict("times" => 10years, "burnin" => 0year, "interval" => 1month, "timestep" => 1month, "scenarios" => scenario, "divfuns" => divfuns, "q" => q, "reps" => 1, "scenario_names" => scenario_names, "cacheInterval" => 1years)
    lensim = length(0month:simDict["interval"]:simDict["times"])

    abenv1 = simplehabitatAE(299.0K, grd, totalK[1], area)
    abenv2 = simplehabitatAE(299.0K, grd, totalK[2], area)
    bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    vars = fill(2.0, numSpecies) .* K
    opts = fill(298.0, numSpecies) .* K

    av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)

    death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
    birth_rates = death_rates

    paramDict = Dict("numSpecies" => numSpecies, "numInvasive" => 0, "numIndiv" => individuals, "reqs" => req, "opts" => opts, "vars" => vars, "birth" => birth_rates, "death" => death_rates, "s" => 5e-2, "boost" => 1.0, "kernel" => kernel, "totalK" => totalK, "bound" => Torus())

    diver = zeros(length(divfuns), lensim, length(scenario))
    runsim!(diver, abenv, paramDict, simDict, 1, cachefolder)
end

endabun1 = sum(JLD.load(joinpath(cachefolder, "cache/GridSize110001.jld"), "abun"))
endabun2 = sum(JLD.load(joinpath(cachefolder, "cache/GridSize210001.jld"), "abun"))
endabun3 = sum(JLD.load(joinpath(cachefolder, "cache/GridSize510001.jld"), "abun"))
endabun4 = sum(JLD.load(joinpath(cachefolder, "cache/GridSize1010001.jld"), "abun"))

bar(["1","4","25","100"], [endabun1, endabun2, endabun3, endabun4], grid = false, xlab = "Number of grid squares", ylab = "Total abundance", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, label = "")

## ABUNDANCE SCALES WITH AREA ##
# Check that larger areas support larger numbers of individuals.

for i in [10.0,20.0,50.0,100.0]
    numSpecies = 100; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2); individuals = 100_000_000; area = i.*km^2; totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

    scenario = [SimpleScenario(TempIncrease!, 0.0K/10year)]

    divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson, mean_abun, geom_mean_abun, norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta, norm_meta_rho, raw_meta_rho, meta_gamma]
    q = 1.0
    a = Int64(i)
    scenario_names = ["Area$a"]

    simDict = Dict("times" => 10years, "burnin" => 0year, "interval" => 1month, "timestep" => 1month, "scenarios" => scenario, "divfuns" => divfuns, "q" => q, "reps" => 1, "scenario_names" => scenario_names, "cacheInterval" => 1years)
    lensim = length(0month:simDict["interval"]:simDict["times"])

    abenv1 = simplehabitatAE(299.0K, grd, totalK[1], area)
    abenv2 = simplehabitatAE(299.0K, grd, totalK[2], area)
    bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    vars = fill(2.0, numSpecies) .* K
    opts = fill(298.0, numSpecies) .* K

    av_dist = rand(Uniform(0.6, 2.4), numSpecies) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)

    death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies)) ./year
    birth_rates = death_rates

    paramDict = Dict("numSpecies" => numSpecies, "numInvasive" => 0, "numIndiv" => individuals, "reqs" => req, "opts" => opts, "vars" => vars, "birth" => birth_rates, "death" => death_rates, "s" => 5e-2, "boost" => 1.0, "kernel" => kernel, "totalK" => totalK, "bound" => Torus())

    diver = zeros(length(divfuns), lensim, length(scenario))
    runsim!(diver, abenv, paramDict, simDict, 1, cachefolder)
end

endabun1 = sum(JLD.load(joinpath(cachefolder, "cache/Area1010001.jld"), "abun"))
endabun2 = sum(JLD.load(joinpath(cachefolder, "cache/Area2010001.jld"), "abun"))
endabun3 = sum(JLD.load(joinpath(cachefolder, "cache/Area5010001.jld"), "abun"))
endabun4 = sum(JLD.load(joinpath(cachefolder, "cache/Area10010001.jld"), "abun"))

bar(["10","20","50","100"], [endabun1, endabun2, endabun3, endabun4], grid = false, xlab = "Area (km²)", ylab = "Total abundance", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, label = "")

## DISPERSAL ##
# 2 species with different dispersal distances seeded at opposite ends of the ecosystem. How long does it take them to mix?

distances = [0.5, 1.0, 2.0, 4.0]
for i in eachindex(distances)
    numSpecies = 2; grd = (10,10); req=(450000.0kJ/m^2, 192.0nm/m^2); individuals = 100_000_000; area = 100.0*km^2; totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

    scenario = [SimpleScenario(TempIncrease!, 0.0K/10year)]

    divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson, mean_abun, geom_mean_abun, norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta, norm_meta_rho, raw_meta_rho, meta_gamma]
    q = 1.0
    a = Int64(i)
    scenario_names = ["Dispersal$a"]

    simDict = Dict("times" => 10years, "burnin" => 0year, "interval" => 1month, "timestep" => 1month, "scenarios" => scenario, "divfuns" => divfuns, "q" => q, "reps" => 1, "scenario_names" => scenario_names, "cacheInterval" => 1years)
    lensim = length(0month:simDict["interval"]:simDict["times"])

    abenv1 = simplehabitatAE(299.0K, grd, totalK[1], area)
    abenv2 = simplehabitatAE(299.0K, grd, totalK[2], area)
    bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    vars = fill(2.0, numSpecies) .* K
    opts = fill(298.0, numSpecies) .* K

    av_dist = fill(distances[i], numSpecies) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)

    death_rates = fill(0.15, numSpecies) ./year
    birth_rates = death_rates

    paramDict = Dict("numSpecies" => numSpecies, "numInvasive" => 0, "numIndiv" => individuals, "reqs" => req, "opts" => opts, "vars" => vars, "birth" => birth_rates, "death" => death_rates, "s" => 5e-2, "boost" => 1.0, "kernel" => kernel, "totalK" => totalK, "bound" => NoBoundary())

    diver = zeros(length(divfuns), lensim, length(scenario))
    dispersalrun!(diver, abenv, paramDict, simDict, 1, cachefolder)
end

display(heatmap(grid = false, xlab = "Distance (km)", ylab = "Distance (km)", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=10, margin = 2.0*Plots.mm, legendfontsize = 12, label = "", layout = (@layout [a b; c d]), clim = (0, 1e6), link = :both))
for i in 1:4
    endabun = mapslices(sum, JLD.load(joinpath(cachefolder, "cache/Dispersal$i"*"10001.jld"), "abun"), dims = 1)[1,:]
    endabun = reshape(endabun, 10, 10)
    m = distances[i]
    display(heatmap!(1:10,1:10, endabun, grid = false, xlab = "Distance (km)", ylab = "Distance (km)", size = (1200, 800), guidefontsize = 12,tickfontsize= 12, titlefontsize=18, margin = 2.0*Plots.mm, legendfontsize = 12, label = "", layout = (@layout [a b; c d]), subplot = i, clim = (0, 1e6), link = :both, title ="Average dispersal of $m km"))
end
