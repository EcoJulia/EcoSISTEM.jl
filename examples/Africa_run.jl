@info "Total Memory: $(Sys.total_memory() / 2^30)GB"
@info "Num threads: $(Threads.nthreads())"

@assert (Sys.total_memory() / 2^30)>=100 "You do not have enough memory to run these examples!"

# Load packages
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Distances
using StatsBase
using Plots
using LinearAlgebra
using JLD2
using Printf
using Diversity
using Diversity.Ecology

const AFRICA_FILE = joinpath(@__DIR__, "Africa.tif")
const SAVEDIR = "/mnt/data/project0000/outputs/Africa_run"
mkpath(SAVEDIR)

function get_active_circle(africa, radius)
    active = Matrix{Bool}(fill(false, size(africa)))
    xs = Base.axes(africa, 1)
    ys = Base.axes(africa, 2)
    for x in xs
        for y in ys
            if norm((x - radius, y - radius)) < radius
                active[x, y] = true
            end
        end
    end
    return active
end

"""
    function run_single(africa, active)

Run a single species model on a circle the size of Africa and save as
a 100-year animated gif.
"""
function run_single(africa, active;
                    savedir = SAVEDIR)

    # Set up initial parameters for ecosystem
    num_species = 1
    grid = size(africa)
    req = 10.0kJ
    individuals = 0
    area = 64e6km^2
    totalK = 1000.0kJ / km^2

    # Set up how much energy each species consumes
    energy_vec = SolarRequirement(fill(req, num_species))

    # Set rates for birth and death
    birth = 0.6 / year
    death = 0.6 / year
    longevity = 1.0
    survival = 0.0
    boost = 1.0
    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(15.0km, 10e-10), num_species)
    movement = AlwaysMovement(kernel, Torus())

    # Create species list, including their temperature preferences, seed abundance and native status
    opts = fill(274.0K, num_species)
    vars = fill(0.5K, num_species)
    traits = GaussTrait(opts, vars)
    native = fill(true, num_species)
    # abun = rand(Multinomial(individuals, num_species))
    abun = fill(div(individuals, num_species), num_species)
    sppl = SpeciesList(num_species, traits, abun, energy_vec,
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
    times = 100years
    timestep = 1month
    record_interval = 1month

    lensim = length((0years):record_interval:times)
    abuns = zeros(Int64, num_species, prod(grid), lensim)
    @time simulate_record!(abuns, eco, times, record_interval, timestep)

    abuns = reshape(abuns[1, :, :, 1], grid[1], grid[2], lensim)

    anim = @animate for i in 1:lensim
        africa_abun = Float64.(abuns[:, :, i])
        africa_abun[.!(active)] .= NaN
        heatmap(africa_abun, clim = (0, 700_000), background_color = :lightblue,
                background_color_outside = :white, grid = false,
                color = cgrad(:algae, scale = :exp), aspect_ratio = 1)
    end
    gif(anim, joinpath(savedir, "Africa.gif"), fps = 30)
    return nothing
end

"""
    function specialist_vs_generalist(africa, active)

Run a model of a specialist invading a generalist on a circle the size of
Africa and save as a pdf.
"""
function specialist_vs_generalist(africa, active;
                                  savedir = SAVEDIR)
    specialist_vars = [0.5K, 1.0K, 5.0K, 10.0K, 25.0K, 50.0K]
    velocity = zeros(typeof(1.0km / year), length(specialist_vars))
    rand_start = rand(findall(active), 1)[1]
    for i in eachindex(specialist_vars)
        # Set up initial parameters for ecosystem
        num_species = 2
        grid = size(africa)
        req = 10.0kJ
        individuals = 0
        area = 64e6km^2
        totalK = 1000.0kJ / km^2

        # Set up how much energy each species consumes
        energy_vec = SolarRequirement(fill(req, num_species))

        # Set rates for birth and death
        birth = 0.6 / year
        death = 0.6 / year
        longevity = 1.0
        survival = 0.1
        boost = 1.0
        # Collect model parameters together
        param = EqualPop(birth, death, longevity, survival, boost)

        # Create kernel for movement
        kernel = fill(GaussianKernel(15.0km, 10e-10), num_species)
        movement = AlwaysMovement(kernel, Torus())

        # Create species list, including their temperature preferences, seed abundance and native status
        opts = fill(274.0K, num_species)
        vars = [50.0K, specialist_vars[i]]
        @debug "Generalist $(vars[1]), specialist $(vard[2])"
        traits = GaussTrait(opts, vars)
        native = fill(true, num_species)
        # abun = rand(Multinomial(individuals, num_species))
        abun = fill(div(individuals, num_species), num_species)
        sppl = SpeciesList(num_species, traits, abun, energy_vec,
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
        burnin = 100years
        times = 100years
        timestep = 1month
        record_interval = 1month

        lensim = length((0years):record_interval:times)
        @time simulate!(eco, burnin, timestep)
        eco.abundances.grid[2, rand_start[1], rand_start[2]] = 100
        abuns = zeros(Int64, num_species, prod(grid), lensim)
        @time simulate_record!(abuns, eco, times, record_interval, timestep)

        abuns = reshape(abuns[:, :, :, 1], num_species, grid[1], grid[2],
                        lensim)
        origin = [rand_start[1], rand_start[2]]
        dest = findall(abuns[2, :, :, end] .> 0)
        dists = [euclidean(origin, [dest[i][1], dest[i][2]])
                 for i in eachindex(dest)] .* getgridsize(eco)
        velocity[i] = mean(dists) / 100.0years
        # inst_velocity = map(1:lensim) do t
        #     dest = findall(abuns[2, :, :, t] .> 0)
        #     dists = [euclidean(origin, [dest[i][1], dest[i][2]]) for i in eachindex(dest)] .* getgridsize(eco)
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
        plot(ustrip.(abs.(specialist_vars .- 50.0K)), ustrip.(velocity),
             xlab = "Selective advantage",
             ylab = "Average invasion speed (km/year)", label = "",
             grid = false)
        Plots.pdf(joinpath(savedir, "InvasionCircle.pdf"))
    end
    return nothing
end

"""
    function specialist_vs_many(africa, active)

Run a model of a specialist invading many (default 50,000) generalisations
on a circle the size of Africa and save as JLD2 files.
"""
function specialist_vs_many(africa, active, num_species = 50_000;
                            savedir = SAVEDIR)
    # Set up initial parameters for ecosystem
    grid = size(africa)
    req = 10.0kJ
    individuals = 3 * 10^8
    area = 64e6km^2
    totalK = 1000.0kJ / km^2

    # Set up how much energy each species consumes
    energy_vec = SolarRequirement(fill(req, num_species))

    # Set rates for birth and death
    birth = 0.6 / year
    death = 0.6 / year
    longevity = 1.0
    survival = 0.1
    boost = 1.0
    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(15.0km, 10e-10), num_species)
    movement = AlwaysMovement(kernel, Torus())

    # Create species list, including their temperature preferences, seed abundance and native status
    opts = fill(274.0K, num_species)
    vars = fill(50.0K, num_species)
    vars[50_000] = 0.5K
    traits = GaussTrait(opts, vars)
    native = fill(true, num_species)
    # abun = rand(Multinomial(individuals, num_species))
    abun = fill(div(individuals, num_species), num_species)
    sppl = SpeciesList(num_species, traits, abun, energy_vec,
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
    burnin = 100years
    times = 100years
    timestep = 1month
    record_interval = 12months
    lensim = length((0years):record_interval:times)
    @time simulate!(eco, burnin, timestep)
    rand_start = rand(findall(active), 1)[1]
    eco.abundances.grid[50_000, rand_start[1], rand_start[2]] = 100
    mkpath(joinpath(savedir, "specialist"))
    @time simulate!(eco, times, timestep, record_interval,
                    joinpath(savedir, "specialist"), "Africa_run")
    return nothing
end

"""
    function specialist_vs_many(africa, active)

Run a model of many (default 50,000) generalists on a circle the size of
Africa and save as a pdf.
"""
function run_many(africa, active, num_species = 50_000; savedir = SAVEDIR)
    # Set up initial parameters for ecosystem
    grid = size(africa)
    req = 10.0kJ
    individuals = 3 * 10^8
    area = 64e6km^2
    totalK = 1000.0kJ / km^2

    # Set up how much energy each species consumes
    energy_vec = SolarRequirement(fill(req, num_species))

    # Set rates for birth and death
    birth = 0.6 / year
    death = 0.6 / year
    longevity = 1.0
    survival = 0.1
    boost = 1.0
    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(15.0km, 10e-10), num_species)
    movement = AlwaysMovement(kernel, Torus())

    # Create species list, including their temperature preferences, seed abundance and native status
    opts = fill(274.0K, num_species)
    vars = fill(50.0K, num_species)
    traits = GaussTrait(opts, vars)
    native = fill(true, num_species)
    # abun = rand(Multinomial(individuals, num_species))
    abun = fill(div(individuals, num_species), num_species)
    sppl = SpeciesList(num_species, traits, abun, energy_vec,
                       movement, param, native)
    sppl.params.birth

    # Create abiotic environment - even grid of one temperature
    abenv = simplehabitatAE(274.0K, grid, totalK, area, active)

    # Set relationship between species and environment (gaussian)
    rel = Gauss{typeof(1.0K)}()

    #Create ecosystem
    eco = Ecosystem(sppl, abenv, rel)

    # Simulation Parameters
    burnin = 10years
    times = 100years
    timestep = 1month
    record_interval = 12months
    lensim = length((0years):record_interval:times)
    @time simulate!(eco, burnin, timestep)
    mkpath(joinpath(savedir, "coexist"))
    @time simulate!(eco, times, timestep, record_interval,
                    joinpath(savedir, "coexist"), "Africa_run_coexist")

    @load joinpath(savedir, "coexist", "Africa_run_coexist100.jld2") abun
    meta = Metacommunity(abun)
    diver = norm_sub_alpha(meta, 0)
    sumabuns = reshape(diver[!, :diversity], 100, 100)

    heatmap(sumabuns,
            grid = false, color = :algae,
            aspect_ratio = 1, layout = (@layout [a b; c d]),
            clim = (0, 50_000), margin = 0.5 * Plots.mm,
            title = "A", titleloc = :left)

    @load joinpath(savedir, "specialist", "Africa_run50.jld2") abun
    meta = Metacommunity(abun)
    diver = norm_sub_alpha(meta, 0)
    sumabuns = reshape(diver[!, :diversity], 100, 100)
    heatmap!(sumabuns,
             background_color = :lightblue,
             background_color_outside = :white,
             grid = false, color = :algae,
             aspect_ratio = 1, subplot = 2,
             clim = (0, 50_000), right_margin = 2.0 * Plots.mm,
             title = "B", titleloc = :left)

    @load joinpath(savedir, "specialist", "Africa_run100.jld2") abun
    meta = Metacommunity(abun)
    diver = norm_sub_alpha(meta, 0)
    sumabuns = reshape(diver[!, :diversity], 100, 100)
    heatmap!(sumabuns,
             background_color = :lightblue,
             background_color_outside = :white,
             grid = false, color = :algae,
             aspect_ratio = 1, subplot = 3,
             clim = (0, 50_000), right_margin = 2.0 * Plots.mm,
             title = "C", titleloc = :left)

    @load joinpath(savedir, "specialist", "Africa_run50.jld2") abun
    meta = Metacommunity(abun)
    diver = shannon(meta)
    sumabuns = reshape(diver[!, :diversity], 100, 100)
    heatmap!(sumabuns,
             background_color = :lightblue,
             background_color_outside = :white,
             grid = false, color = :algae,
             aspect_ratio = 1, subplot = 4,
             right_margin = 2.0 * Plots.mm,
             title = "D", titleloc = :left, clim = (0, 10))
    Plots.pdf(joinpath(savedir, "Africa.pdf"))
    return nothing
end

const AFRICA_TIF = readfile(AFRICA_FILE, -25.0°, 50.0°, -35.0°, 40.0°)
const RADIUS = 50
const ACTIVE = get_active_circle(AFRICA_TIF, RADIUS)

heatmap(ACTIVE)

#### SINGLE SPECIES ####
@info "Run a single species"
flush(stdout)
run_single(AFRICA_TIF, ACTIVE)

#### SPECIALIST VERSUS GENERALIST ####
@info "Run a single specialist against a generalist"
flush(stdout)
specialist_vs_generalist(AFRICA_TIF, ACTIVE)

#### SPECIALIST VERSUS MANY GENERALISTS ####
@info "Run a single specialist against many generalists"
flush(stdout)
specialist_vs_many(AFRICA_TIF, ACTIVE)

#### 50,000 SPECIES COEXISTING #####
@info "Run many generalists"
flush(stdout)
run_many(AFRICA_TIF, ACTIVE)
