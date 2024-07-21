# SPDX-License-Identifier: LGPL-3.0-or-later

### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end
    
# ╔═╡ 108951ec-3ecf-4f0c-b5e4-e79d00b1bfac
begin
    using EcoSISTEM
    using EcoSISTEM.Units
    using EcoSISTEM.ClimatePref
    using Unitful
    using Unitful.DefaultSymbols
    using Diversity
    using OnlineStats
    using Plots
    using Distributions
    using Diversity
    using SpatialEcology
    using RasterDataSources
    using AxisArrays
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
    file = "../examples/Africa.tif"
    africa = readfile(file, -25.0°, 50.0°, -35.0°, 40.0°)
    active = Matrix{Bool}(.!isnan.(africa))
    # Set up initial parameters for ecosystem
    grd = size(africa)
    req = 10.0kJ
    individuals = 3 * 10^8
    area = 64e6km^2
    totalK = 1000.0kJ / km^2

    # Set up how much energy each species consumes
    energy_vec = SolarRequirement(fill(req, numSpecies))

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
    movement = AlwaysMovement(kernel, Torus())

    # Create species list, including their temperature preferences, seed abundance and native status
    opts = fill(274.0K, numSpecies)
    vars = fill(50.0K, numSpecies)
    vars[end] = (1 / niche_width) * K
    trts = GaussTrait(opts, vars)
    native = fill(true, numSpecies)
    abun = fill(div(individuals, numSpecies), numSpecies)
    sppl = SpeciesList(numSpecies, trts, abun, energy_vec,
                       movement, param, native)
    sppl.params.birth

    # Create abiotic environment - even grid of one temperature
    abenv = simplehabitatAE(274.0K, grd, totalK, area, active)

    # Set relationship between species and environment (gaussian)
    rel = Gauss{typeof(1.0K)}()

    # Create ecosystem
    eco = Ecosystem(sppl, abenv, rel)
    eco.abundances.matrix[end, :] .= 0

    # Simulation Parameters
    burnin = 1years
    times = 10years
    timestep = 1month
    record_interval = 1month
    lensim = length((0years):record_interval:times)
    abuns = zeros(numSpecies, prod(grd), lensim)

    # Run simulation for burnin and then add invasive species
    @time simulate!(eco, burnin, timestep)
    eco.abundances.grid[end, 50, 50] = 100
    @time simulate_record!(abuns, eco, times, timestep, record_interval)

    plot(eco, clim = (0, numSpecies))
end

# ╔═╡ 7e16f197-874b-482d-80b6-13a62ddda1f7
begin
    if !isdir("assets")
        mkdir("assets")
    end
    ENV["RASTERDATASOURCES_PATH"] = "assets"
    getraster(WorldClim{BioClim})
    world = readbioclim("assets/WorldClim/BioClim/")
    africa_temp = world.array[-25° .. 50°, -35° .. 40°, 1]
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
    temp = uconvert.(K, africa_temp .* °C)
    africa_new = Worldclim_bioclim(AxisArray(temp,
                                             AxisArrays.axes(africa_temp)))
    active_new = Matrix{Bool}(.!isnan.(africa))

    # Set up initial parameters for ecosystem
    grd_new = size(africa_new.array)
    req_new = 10.0kJ
    individuals_new = 3 * 10^8
    area_new = 64e6km^2
    totalK_new = 1000.0kJ / km^2

    # Set up how much energy each species consumes
    energy_vec_new = SolarRequirement(fill(req, numSpecies))

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
    movement_new = AlwaysMovement(kernel, Torus())

    # Create species list, including their temperature preferences, seed abundance and native status
    opts_new = fill(meantemp * 1.0K, numSpecies)
    vars_new = fill(5.0K, numSpecies)
    trts_new = GaussTrait(opts_new, vars_new)
    native_new = fill(true, numSpecies)
    abun_new = fill(div(individuals_new, numSpecies), numSpecies)
    sppl_new = SpeciesList(numSpecies, trts_new, abun_new, energy_vec_new,
                           movement_new, param_new, native_new)
    sppl.params.birth

    # Create abiotic environment - even grid of one temperature
    abenv_new = bioclimAE(africa_new, totalK_new, area_new)

    # Set relationship between species and environment (gaussian)
    rel_new = Gauss{typeof(1.0K)}()

    # Create ecosystem
    eco_new = Ecosystem(sppl_new, abenv_new, rel_new)

    # Simulation Parameters
    burnin_new = 1years
    times_new = 10years
    timestep_new = 1year
    record_interval_new = 1year
    lensim_new = length((0years):record_interval_new:times_new)
    abuns_new = zeros(numSpecies, prod(grd_new), lensim_new)

    # Run simulation for burnin and then add invasive species
    @time simulate!(eco_new, burnin_new, timestep_new)
    @time simulate_record!(abuns_new, eco_new, times_new, timestep_new,
                           record_interval_new)

    mean_abuns = reshape(mean(eco_new.abundances.matrix, dims = 1)[1, :],
                         grd_new)
    mean_abuns[.!eco_new.abenv.active] .= NaN
    heatmap(mean_abuns')
end

# ╔═╡ Cell order:
# ╟─39b75180-f384-11eb-3449-4f7c2ad25d99
# ╠═108951ec-3ecf-4f0c-b5e4-e79d00b1bfac
# ╠═e9e9065b-cf7e-4d37-8162-88f0076ad1eb
# ╠═220e4af6-f228-4b8d-a77e-0ddbf5fc6705
# ╠═7e16f197-874b-482d-80b6-13a62ddda1f7
# ╠═f9d43c58-4888-402e-873a-81f3c4ffd367
# ╠═ee925e21-b0b6-478e-a3a0-573e8497b9f6
