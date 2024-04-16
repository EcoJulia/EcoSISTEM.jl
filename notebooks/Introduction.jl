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

# ╔═╡ 496434ae-1bc8-4472-8981-2d9fc08b688a
# Load in necessary packages
begin
    using EcoSISTEM
    using EcoSISTEM.Units
    using Unitful
    using Unitful.DefaultSymbols
    using Diversity
    using OnlineStats
    using Plots
    using Distributions
    using SpatialEcology
end

# ╔═╡ 6b280ae0-6af9-4d3c-91c8-f02a20d28d18
md"# _A quick introduction to EcoSISTEM.jl!_

Welcome to the EcoSISTEM package. You can use us to simulate communities of plants! In any simulation, we run multiple species competing, dispersing and dying over time and space. You can set up any type of Ecosystem you like, but there are a few in-built functions to make things easier.

Have a go below and remember you can write `?` to check the documentation for any Julia function or variable."

# ╔═╡ 39be3104-92c0-40fe-ad10-905d6c889917
begin
    spp_slider = @bind numSpecies html"<input type='range' min='2' max='10' step='1' value='5'>"
    grid_slider = @bind numGrid html"<input type='range' min='1' max='10' step='1' value='5'>"

    md"""**Please set model parameters:**

    Number of species: $(spp_slider)

    Number of grid cells: $grid_slider
    """
end

# ╔═╡ fa91699e-926f-452b-bae0-219fe489709f
md"You will see: (A) the abundances of each species over time, and (B) the total abundance of all species over space at the end of the simulation. These simulations are stochastic, so the exact numbers will differ between runs.
"

# ╔═╡ 368a2a64-2aee-4282-b182-e5182b949392
begin

    # Set up abiotic environment
    grd = (numGrid, numGrid)
    area = 100.0 * km^2
    totalK = (4.5e11kJ / km^2, 192.0mm / km^2)
    abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
    abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
    bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat,
                                                                abenv1.active,
                                                                bud,
                                                                abenv1.names)

    # Species characteristics
    individuals = 100_000
    abun = rand(Multinomial(individuals, numSpecies))
    death = 0.15 / year
    birth = death
    long = 1.0
    surv = 0.1
    boost = 1.0
    param = EqualPop(birth, death, long, surv, boost)
    native = fill(true, numSpecies)

    # Dispersal
    av_dist = fill(2.4, numSpecies) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)
    movement = BirthOnlyMovement(kernel, Torus())

    # Resource requirements
    req = (450000.0kJ / m^2, 192.0nm / m^2)
    size_mean = 1.0m^2
    energy_vec1 = SolarRequirement(fill(req[1] * size_mean, numSpecies))
    energy_vec2 = WaterRequirement(fill(req[2] * size_mean, numSpecies))
    energy_vec = ReqCollection2(energy_vec1, energy_vec2)

    # Habitat preferences
    vars = fill(2.0, numSpecies) .* K
    opts = 298.0K .+ vars .* range(-3, stop = 3, length = numSpecies)
    trts = GaussTrait(opts, vars)

    # Species list
    sppl = SpeciesList(numSpecies, trts, abun, energy_vec,
                       movement, param, native)

    # Trait relationship
    rel = Gauss{typeof(first(opts))}()

    # Build ecosystem
    eco = Ecosystem(sppl, abenv, rel)

    # Run simulation
    times = 10years
    timestep = 1month
    interval = 1month
    lensim = length((0month):timestep:times)
    abuns = zeros(Int64, numSpecies, prod(grd), lensim)
    simulate_record!(abuns, eco, times, timestep, interval)

    # Plot abundances
    sumabun = sum(abuns, dims = 2)[:, 1, :]
    floor_plot = zeros(lensim)
    p1 = plot(sumabun[1, :], ribbon = (sumabun[1, :], floor_plot),
              guidefontsize = 16, tickfontsize = 16, size = (1000, 1000),
              titlefontsize = 16,
              grid = false, layout = (2, 1), label = "", titleloc = :left)
    for i in 2:(numSpecies)
        before_spp = sum(sumabun[1:(i - 1), :], dims = 1)[1, :]
        p1 = plot!(sumabun[i, :] .+ before_spp,
                   ribbon = (sumabun[i, :], floor_plot),
                   xlabel = "Time", ylabel = "Population", label = "",
                   title = "A", subplot = 1)
    end
    sumabun_space = sum(eco.abundances.matrix, dims = 1)[1, :]
    p1 = heatmap!(reshape(sumabun_space, grd[1], grd[2]), subplot = 2,
                  aspect_ratio = 1, axis = false, xlabel = "Grid location",
                  ylabel = "Grid location", guidefontsize = 16,
                  tickfontsize = 16,
                  titlefontsize = 16, grid = false,
                  right_margin = 20.0 * Plots.mm, label = "", title = "B",
                  titleloc = :left)
    p1
end

# ╔═╡ ff15c0e9-2fee-4c42-a255-924cd29c5ecb
md"
_P.S. you'll notice that the number of grid cells doesn't really change the results of plot A - this is because the overall area remains the same and the simulation is built to be invariant to grid size. Neat!_
"

# ╔═╡ 288fb7fb-60d5-4e4c-914c-39eaa54577d9
md"## The nuts and bolts

If you have taken a peak at the code to plot the results above (by clicking the eye icon), you'll notice it is fairly dense! Let's go through it bit by bit and see how we built an `Ecosystem`, the base type which holds all the information about species and their environment.

### _Let's start with the environment_

Something simple perhaps! We'll build a 10 by 10 grid, because in our world everything is gridded. This abiotic environment will house information on the habitat (like climate) and what resources are available (e.g. sunlight and water for the plants). It is its own Julia type, a `GridAbioticEnv`. For a simple example, we'll build a grid of temperature as the habitat and water as the resource:
"

# ╔═╡ 16a70493-683a-45bc-baa0-83e917fa774b
begin
    # Here is our grid
    grid = (10, 10)

    # The total area we would like it to be
    area_size = 100.0 * km^2

    # The total water it will have available
    totalW = 200.0mm / km^2

    # Overall temperature it will be
    totalT = 298.0K

    # Perfect, now we can build a simple habitat!
    temp_env = simplehabitatAE(totalT, grid, totalW, area_size)

    # Let's plot it to see what it looks like
    heatmap(temp_env.habitat.matrix ./ K, clim = (278, 308), title = "Habitat",
            layout = 2)
    heatmap!(temp_env.budget.matrix ./ mm, title = "Resource", subplot = 2,
             clim = (100, 200))
end

# ╔═╡ d1602be3-2724-4441-8dc3-8fa7713ed249
md"Very boring! But we have succeeded in making a simple environment with a temperature of 298K (or 25 degrees C) with 200 mm of water in each square.

Let's try something slightly more adventurous:"

# ╔═╡ e1ca4d60-383f-4c2a-945b-c7665b51bff9
begin
    # A temperature gradient spanning 10 degrees either side of total temperature from above. We can also give it a rate over which to change over time
    temp_change_rate = 0.2K / month
    temp_grad_env = tempgradAE(totalT - 10.0K, totalT + 10.0K, grid, totalW,
                               area_size, temp_change_rate)

    # Let's plot it to see what it looks like now
    heatmap(temp_grad_env.habitat.matrix' ./ K, clim = (278, 308),
            title = "Habitat")
end

# ╔═╡ 8f7b85a1-1e64-4d1a-90f6-55c6307a03cf
md"That's better! How about something even fancier?"

# ╔═╡ 01a11afc-7c9b-41f5-b308-003303dfa72a
begin
    # A temperature peak spanning 10 degrees either side of total temperature from above. We can also give it a rate over which to change over time
    temp_peak_env = peakedgradAE(totalT - 10.0K, totalT + 10.0K, grid, totalW,
                                 area_size, temp_change_rate)

    # Let's plot it to see what it looks like now
    heatmap(temp_peak_env.habitat.matrix' ./ K, clim = (278, 308),
            title = "Habitat")
end

# ╔═╡ 5928ca0d-479e-4714-acfe-ff2d9c43533e
md"
### _Now to add the species!_

Similarly to the abiotic environment, all information about species is housed in its own Julia type, called a `SpeciesList`. We'll set one up now:
"

# ╔═╡ 797cc60f-ce47-4d64-b7d9-1fa16a74db7a
begin
    # First things first - how many species?
    numSpp = 10

    # And how many individuals overall
    numInd = 100_000

    # Great, now we'll randomly allocate those individuals evenly (but randomly!) 
    # across the species
    start_abuns = rand(Multinomial(numInd, numSpp))

    # Are they all native?
    is_native = fill(true, numSpp)

    # Next up, we need to decide on some base rates for the simulation to run on - 
    # like birth and death
    death_rates = 0.15 / year
    birth_rates = death_rates

    # There are also a couple of parameters that should be set about how the species 
    # lifespan relates to their resource consumption (longevity) and how well they 
    # survive in their current environment based on their traits (survival)
    longevity = 1.0
    survival = 0.1

    # Finally, how many times boost can they get to their reproduction from available 
    # resources?
    resource_boost = 1.0

    # Let's assume that all species are equal for these parameters
    parameters = EqualPop(birth_rates, death_rates, longevity, survival,
                          resource_boost)

    # Now we have to consider movement - let's say the species move an average of 
    # 2.4km according to a Gaussian kernel
    average_dist = fill(2.4, numSpp) .* km
    gauss_kernel = GaussianKernel.(average_dist, 10e-10)

    # Because we are considering plants, let's assume the seed production is combined 
    # with dispersal, so only new seeds move
    move = BirthOnlyMovement(gauss_kernel, Torus())

    # We must also decide how much water each species needs per timestep
    water_req = (100.0nm / m^2)
    size = 1.0m^2
    water_vec = WaterRequirement(fill(water_req * size, numSpp))

    # Plus, their niche width - the range of habitats they find suitable
    niche_width = fill(2.0, numSpp) .* K

    # And what is their niche optimum - here we have a range around 25 degrees
    optima = 298.0K .+ niche_width .* range(-3, stop = 3, length = numSpp)

    # Let's combine these pieces of information together so that each species has a 
    # Gaussian curve representing their match to the current temperature
    traits = GaussTrait(optima, niche_width)

    # Okay, now we are ready for our species list!
    species_list = SpeciesList(numSpp, traits, start_abuns, water_vec,
                               move, parameters, is_native)
end

# ╔═╡ 8e33c75b-d128-4b25-9b2b-5cff452c3a4e
md"Phew, that was a lot, but now we can look at some of those pieces of information we put in. For example, the match of the different species to their environment. "

# ╔═╡ b9091796-d377-49c7-8601-a84ec0006f65
begin
    p3 = plot(288:0.1:308,
              pdf.(Normal(optima[1] / K, niche_width[1] / K),
                   collect(288:0.1:308)),
              label = "", xlabel = "Temperature (K)",
              ylabel = "Match to environment")
    for j in 2:numSpp
        plot!(288:0.1:308,
              pdf.(Normal(optima[j] / K, niche_width[j] / K),
                   collect(288:0.1:308)),
              subplot = 1, label = "")
    end
    p3 = vline!([298], colour = :black, label = "Current temperature")
    p3
end

# ╔═╡ ccbf7d75-1c1d-43ad-998e-aa033bcd1ebc
md"### _The final touches ..._

The last thing we need to specify is the relationship between the species and their environment - for this example, let's assume that we want this match to be a Gaussian curve like above.
"

# ╔═╡ 41206ea4-77ed-4b87-8acf-8d2a9ee170db
trait_relationship = Gauss{typeof(first(optima))}()

# ╔═╡ 1d59bcaa-3670-4765-8981-9e2dd6d99436
md"Now we can build our ecosystem! EcoSISTEM is integrated with SpatialEcology.jl, so we can take advantage of their plotting system, by simply calling plot on our `Ecosystem` object. This will show us the species richness over space."

# ╔═╡ ffd5ca3f-00b3-4444-9f89-77062358c439
example_eco = Ecosystem(species_list, temp_grad_env, trait_relationship)

# ╔═╡ 7492f556-957c-48a3-8597-9c7c7ba20547
begin
    using Diversity.Ecology
    shannon(example_eco.abundances.matrix)
end

# ╔═╡ ca77e96e-dda3-4c5a-b063-e782535045e5
plot(example_eco, clim = (0, 10), title = "Species richness")

# ╔═╡ 43b32242-7577-4593-b8ea-9f3fcd0cfda5
md"We can see that all species are present in the grid as we have set it up. Let's simulate it over a 10 year period and monthly timesteps to see what will happen. Remember that we added in a change in temperature over time!"

# ╔═╡ fa1ea836-6750-43d0-b574-d1490ecd6ebf
simulate!(example_eco, 10year, 1month)

# ╔═╡ 6b0fb54e-6d52-40ed-87d5-96aff5f3d93d
plot(example_eco, clim = (0, 10))

# ╔═╡ 46bd023d-d3d7-4b8f-b8da-ef1800a2b939
md"Our first simulation run! You'll notice that in the hotter locations, we see species richness drop more quickly.

We can also play around with other measures of diversity through the Diversity.jl package. Here we have a measure of representativeness, called rho diversity, for every grid cell. You'll notice that as we go down the list, grid cells get less representative, as they contain fewer and more rare species."

# ╔═╡ 858fd79e-d41e-4302-8197-c6c628919844
norm_sub_rho(example_eco, 1.0)

# ╔═╡ 078406b4-7546-4da7-a13a-57acbd1f6925
md"We can also access more traditional diversity measures like Shannon and Simpson through the `Diversity.Ecology` module."

# ╔═╡ 81cd8c4b-bf8a-424d-840d-3834beae575b
simpson(example_eco.abundances.matrix)

# ╔═╡ f1abf877-98ce-4d0f-8408-ad6e203c6fbd
md"The final thing you need to know for this tutorial is that you can store the outputs of the model on a certain time interval, like so:"

# ╔═╡ 97a3ef6d-e7ff-4851-95b5-c8482e953b70
begin
    simulation_time = 10years
    time_step = 1month
    record_interval = 1month
    len_sim = length((0month):time_step:simulation_time)
    record_abuns = zeros(Int64, numSpp, prod(grid), len_sim)
    simulate_record!(record_abuns, example_eco, simulation_time, time_step,
                     record_interval)
end

# ╔═╡ 0d0e0197-d826-4f62-943c-79783e5fa701
md"# _That's all folks!_

There are more examples that can be found in the examples folder. Be aware to check first how big the simulation is before running it on your laptop!
"

# ╔═╡ Cell order:
# ╟─6b280ae0-6af9-4d3c-91c8-f02a20d28d18
# ╠═496434ae-1bc8-4472-8981-2d9fc08b688a
# ╠═39be3104-92c0-40fe-ad10-905d6c889917
# ╟─fa91699e-926f-452b-bae0-219fe489709f
# ╟─368a2a64-2aee-4282-b182-e5182b949392
# ╟─ff15c0e9-2fee-4c42-a255-924cd29c5ecb
# ╟─288fb7fb-60d5-4e4c-914c-39eaa54577d9
# ╠═16a70493-683a-45bc-baa0-83e917fa774b
# ╟─d1602be3-2724-4441-8dc3-8fa7713ed249
# ╠═e1ca4d60-383f-4c2a-945b-c7665b51bff9
# ╟─8f7b85a1-1e64-4d1a-90f6-55c6307a03cf
# ╠═01a11afc-7c9b-41f5-b308-003303dfa72a
# ╟─5928ca0d-479e-4714-acfe-ff2d9c43533e
# ╠═797cc60f-ce47-4d64-b7d9-1fa16a74db7a
# ╟─8e33c75b-d128-4b25-9b2b-5cff452c3a4e
# ╠═b9091796-d377-49c7-8601-a84ec0006f65
# ╟─ccbf7d75-1c1d-43ad-998e-aa033bcd1ebc
# ╠═41206ea4-77ed-4b87-8acf-8d2a9ee170db
# ╟─1d59bcaa-3670-4765-8981-9e2dd6d99436
# ╠═ffd5ca3f-00b3-4444-9f89-77062358c439
# ╠═ca77e96e-dda3-4c5a-b063-e782535045e5
# ╟─43b32242-7577-4593-b8ea-9f3fcd0cfda5
# ╠═fa1ea836-6750-43d0-b574-d1490ecd6ebf
# ╠═6b0fb54e-6d52-40ed-87d5-96aff5f3d93d
# ╟─46bd023d-d3d7-4b8f-b8da-ef1800a2b939
# ╠═858fd79e-d41e-4302-8197-c6c628919844
# ╟─078406b4-7546-4da7-a13a-57acbd1f6925
# ╠═7492f556-957c-48a3-8597-9c7c7ba20547
# ╠═81cd8c4b-bf8a-424d-840d-3834beae575b
# ╟─f1abf877-98ce-4d0f-8408-ad6e203c6fbd
# ╠═97a3ef6d-e7ff-4851-95b5-c8482e953b70
# ╟─0d0e0197-d826-4f62-943c-79783e5fa701
