### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 108951ec-3ecf-4f0c-b5e4-e79d00b1bfac
begin
    using EcoSISTEM
    using EcoSISTEM.Units
    using Unitful
    using Unitful.DefaultSymbols
    using Diversity
    using OnlineStats
    using Plots
    using Distributions
    using Diversity
    using SpatialEcology
    using RasterDataSources
    using Rasters: ProjString
end

# ╔═╡ 39b75180-f384-11eb-3449-4f7c2ad25d99
md"# _Interactive Africa!_

Welcome to EcoSISTEM.jl! In the examples folder, you'll notice an example of running simulations of 50,000 plant species across Africa. This is very computationally intensive! Here we have an easy-to-run, interactive example for you to run.

The first time round things may be a little slow - this is because Julia is compiling the code. After that, the time will whizz by!"

# ╔═╡ e9e9065b-cf7e-4d37-8162-88f0076ad1eb
begin
    spp_slider = @bind numSpecies html"<input type='range' min='2' max='10' step='1' value='5'>"
    invasive_slider = @bind niche_width html"<input type='range' min='0.02' max='50' step='10' value='10'>"

    md"""**Please set model parameters:**

    Number of species: $(spp_slider)

    Invasive ability: $invasive_slider
    """
end

# ╔═╡ 220e4af6-f228-4b8d-a77e-0ddbf5fc6705
begin
    file = pkgdir(EcoSISTEM, "data", "Africa.tif")
    africa = readfile(file)

    # The shipped Africa raster is a *landmask*: real values on land, `NaN` at sea. We hand it to
    # the study area as a layer so that its own gaps decide which cells are active - no hand-built
    # `.!isnan.(...)` matrix, and no chance of the mask and the grid disagreeing.
    #
    # `in_memory_raster` is how a raster you already hold becomes a layer spec: a raster carries
    # values but no niche axis, so nothing about it says what its numbers mean. Here they mean
    # nothing in particular - the raster is a shape - so it stays `Unclassified`.
    #
    # It is a *geographic* (WGS 84) raster, and a simulation needs a projected grid: dispersal
    # assumes one uniform cell size, whereas a degree cell shrinks towards the poles. Giving a
    # projected `crs` reprojects it onto one.
    # `readfile` hands back a plain array-with-coordinates; `ClimateRaster` is what pairs it with a
    # statement of where it came from, and `SyntheticData` is the honest answer for a shipped
    # landmask that belongs to no catalogued dataset.
    africa_shape = EcoSISTEM.in_memory_raster(ClimateRaster(EcoSISTEM.SyntheticData,
                                                            africa),
                                              axis = EcoSISTEM.NicheAxis)
    studyarea = StudyArea(regime = africa_shape,
                          crs = ProjString("+proj=aea +lat_1=20 +lat_2=-23 " *
                                           "+lat_0=0 +lon_0=25 +datum=WGS84 " *
                                           "+units=m +no_defs"),
                          cellsize = 100km, verbosity = :silent)
    active = studyarea.report.active

    # Set up initial parameters for ecosystem
    grd = size(active)
    demand = 10.0kJ / day
    individuals = 3 * 10^8
    area = 64e6km^2
    totalK = 1000.0kJ / km^2 / day

    # Set up how much resource each species consumes
    resource_vec = Demand{SolarRadiation}(fill(demand, numSpecies))

    # Set rates for birth and death
    birth = 0.6 / year
    death = 0.6 / year
    longevity = 1.0
    survival = 0.1
    boost = 1.0
    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
    movement = AlwaysMovement(kernel)

    # Create species list, including their temperature preferences, seed abundance and native status
    opts = fill(274.0K, numSpecies)
    vars = fill(50.0K, numSpecies)
    vars[end] = (1 / niche_width) * K
    tolerance = NicheTolerance(Temperature, Normal, opts, vars)
    native = fill(true, numSpecies)
    abun = fill(div(individuals, numSpecies), numSpecies)
    sppl = SpeciesList(numSpecies, tolerance, abun, resource_vec,
                       movement, param, native)
    sppl.params.birth

    # Create abiotic environment - an even grid of one temperature, on the Africa-shaped grid
    # decided above. `GridHabitat` chooses nothing: it samples the layers named here onto
    # the grid the study area already settled, mask and all.
    habitat = GridHabitat(regime = UniformSpec(274.0K,
                                               axis = Temperature),
                          supply = UniformSpec(totalK,
                                               axis = SolarRadiation),
                          area = studyarea)

    # Set nichefit between species and environment (gaussian)
    nichefit = NicheSuitability{Temperature, typeof(1.0K)}()

    # Create ecosystem
    eco = Ecosystem(sppl, habitat, nichefit)
    eco.abundances.matrix[end, :] .= 0

    # Simulation Parameters
    burnin = 1year
    times = 10year
    timestep = 1month_mean_duration
    record_interval = 1month_mean_duration
    lensim = length((0year):record_interval:times)
    abuns = zeros(numSpecies, prod(grd), lensim)

    # Run simulation for burnin and then add invasive species
    @time simulate!(eco, burnin, timestep)
    eco.abundances.grid[end, 50, 50] = 100
    @time simulate_record!(abuns, eco, times, record_interval,
                           timestep)

    # **Not `plot(eco)`.** That recipe throws a `BoundsError` on any grid that is not square -
    # reproduced on a 4×6 grid, where it reaches for `[5, 1]` - so it worked here only while this
    # notebook used the shipped 100×100 Africa raster unprojected. Reprojecting onto a real grid
    # made it 76×67 and the latent transposition surfaced immediately.
    #
    # Mean abundance per cell, masked to the active area - the same shape the section below uses,
    # and one that reads the `(y, x)` grid the right way round.
    mean_abuns = reshape(mean(eco.abundances.matrix, dims = 1)[1, :], grd)
    mean_abuns = Float64.(mean_abuns)
    mean_abuns[.!active] .= NaN
    heatmap(mean_abuns, title = "Mean abundance per cell", grid = false,
            background_color = :lightblue)
end

# ╔═╡ 7e16f197-874b-482d-80b6-13a62ddda1f7
begin
    africa_bio = read(WorldClim{BioClim}, 1,
                      cut = EcoSISTEM.boundingbox("Africa",
                                                  round = 5°))
    africa_temp = africa_bio.array
    plot(africa_temp)
end

# ╔═╡ f9d43c58-4888-402e-873a-81f3c4ffd367
begin
    temp_slider = @bind meantemp html"<input type='range' min='270' max='300' step='10' value='290'>"

    md"""**Please set model parameters:**

    Temperature preference: $(temp_slider)
    """
end

# ╔═╡ ee925e21-b0b6-478e-a3a0-573e8497b9f6
begin
    # The grid is decided before anything is built on it. Real climate data has to be simulated on a
    # *projected* grid: dispersal assumes one uniform cell size, whereas a degree grid's cells shrink
    # towards the poles, so `Ecosystem` refuses a geographic one. Africa Albers Equal Area Conic is
    # the natural choice for the whole continent - equal-area, so a cell really does mean the same
    # amount of ground everywhere, which is what the ecology depends on.
    albers_new = ProjString("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 " *
                            "+lon_0=25 +datum=WGS84 +units=m +no_defs")
    regime_new = SourceSpec(WorldClim{BioClim}, 1)
    totalK_new = 1000.0kJ / km^2 / day
    supply_new = UniformSpec(totalK_new, axis = SolarRadiation)
    studyarea_new = StudyArea(regime = regime_new,
                              within = EcoSISTEM.boundingbox("Africa",
                                                             round = 5°),
                              crs = albers_new, cellsize = 50km,
                              verbosity = :silent)

    # Set up initial parameters for ecosystem
    grd_new = size(studyarea_new.report.active)
    req_new = 10.0kJ / day
    individuals_new = 3 * 10^8

    # Set up how much resource each species consumes
    energy_vec_new = Demand{SolarRadiation}(fill(demand, numSpecies))

    # Set rates for birth and death
    birth_new = 0.6 / year
    death_new = 0.6 / year
    longevity_new = 1.0
    survival_new = 0.1
    boost_new = 1.0
    # Collect model parameters together
    param_new = EqualPop(birth, death, longevity, survival, boost)

    # Create kernel for movement
    kernel_new = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
    movement_new = AlwaysMovement(kernel)

    # Create species list, including their temperature preferences, seed abundance and native status
    opts_new = fill(meantemp * 1.0K, numSpecies)
    vars_new = fill(5.0K, numSpecies)
    trts_new = NicheTolerance(Temperature, Normal, opts_new, vars_new)
    native_new = fill(true, numSpecies)
    abun_new = fill(div(individuals_new, numSpecies), numSpecies)
    sppl_new = SpeciesList(numSpecies, trts_new, abun_new, energy_vec_new,
                           movement_new, param_new, native_new)
    sppl.params.birth

    # Build the environment on that decided grid: the regime is the real bio1 temperature layer,
    # reprojected onto it, and the resource a flat solar supply.
    abenv_new = GridHabitat(regime = regime_new, supply = supply_new,
                            area = studyarea_new,
                            topology = Torus())

    # Set nichefit between species and environment (gaussian)
    rel_new = NicheSuitability{Temperature, typeof(1.0K)}()

    # Create ecosystem
    eco_new = Ecosystem(sppl_new, abenv_new, rel_new)

    # Simulation Parameters
    burnin_new = 1year
    times_new = 10year
    timestep_new = 1year
    record_interval_new = 1year
    lensim_new = length((0year):record_interval_new:times_new)
    abuns_new = zeros(numSpecies, prod(grd_new), lensim_new)

    # Run simulation for burnin and then add invasive species
    @time simulate!(eco_new, burnin_new, timestep_new)
    @time simulate_record!(abuns_new, eco_new, times_new, timestep_new,
                           record_interval_new)

    mean_abuns = reshape(mean(eco_new.abundances.matrix, dims = 1)[1, :],
                         grd_new)
    mean_abuns[.!eco_new.habitat.active] .= NaN
    heatmap(mean_abuns)
end

# ╔═╡ Cell order:
# ╟─39b75180-f384-11eb-3449-4f7c2ad25d99
# ╠═108951ec-3ecf-4f0c-b5e4-e79d00b1bfac
# ╠═e9e9065b-cf7e-4d37-8162-88f0076ad1eb
# ╠═220e4af6-f228-4b8d-a77e-0ddbf5fc6705
# ╠═7e16f197-874b-482d-80b6-13a62ddda1f7
# ╠═f9d43c58-4888-402e-873a-81f3c4ffd367
# ╠═ee925e21-b0b6-478e-a3a0-573e8497b9f6
