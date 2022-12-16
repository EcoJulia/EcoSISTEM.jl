### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
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
	using Diversity
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
	grd = (numGrid, numGrid); area = 100.0*km^2;
	totalK = (4.5e11kJ/km^2, 192.0mm/km^2)
	abenv1 = simplehabitatAE(298.0K, grd, totalK[1], area)
	abenv2 = simplehabitatAE(298.0K, grd, totalK[2], area)
	bud = BudgetCollection2(abenv1.budget, abenv2.budget)
	abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, 		abenv1.active, bud, abenv1.names)
	
	# Species characteristics
	individuals = 100_000
	abun = rand(Multinomial(individuals, numSpecies))
	death = 0.15/ year
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
	req=(450000.0kJ/m^2, 192.0nm/m^2)
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
	times = 10years; timestep = 1month; interval = 1month
	lensim = length(0month:timestep:times)
	abuns = zeros(Int64, numSpecies, prod(grd), lensim)
	simulate_record!(abuns, eco, times, timestep, interval)
	
	# Plot abundances
	sumabun = sum(abuns, dims = 2)[:, 1, :]
	floor_plot = zeros(lensim)
	p1 = plot(sumabun[1,:], ribbon=(sumabun[1,:], floor_plot), 
		guidefontsize = 16, tickfontsize= 16, size = (1000, 1000), titlefontsize = 16, 
		grid= false, layout = (2, 1), label = "", titleloc = :left)
	for i in 2:(numSpecies)
		before_spp = sum(sumabun[1:(i-1), :], dims = 1)[1, :]
		p1 = plot!(sumabun[i, :] .+ before_spp, ribbon=(sumabun[i,:], floor_plot), 
			xlabel="Time", ylabel="Population", label="", 
			title="A", subplot = 1)
	end
	sumabun_space = sum(eco.abundances.matrix, dims = 1)[1, :]
	p1 = heatmap!(reshape(sumabun_space, grd[1], grd[2]), subplot = 2, 
		aspect_ratio = 1, axis = false, xlabel = "Grid location", 
		ylabel = "Grid location", guidefontsize = 16, tickfontsize= 16, 
		titlefontsize = 16, grid= false, right_margin = 20.0 * Plots.mm, label = "", title = "B", titleloc = :left)
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
	area_size = 100.0*km^2

	# The total water it will have available
	totalW = 200.0mm/km^2

	# Overall temperature it will be
	totalT = 298.0K

	# Perfect, now we can build a simple habitat!
	temp_env = simplehabitatAE(totalT, grid, totalW, area_size)

	# Let's plot it to see what it looks like
	heatmap(temp_env.habitat.matrix ./ K, clim = (278, 308), title = "Habitat", layout = 2)
	heatmap!(temp_env.budget.matrix ./ mm, title = "Resource", subplot = 2, clim = (100, 200))
	
end

# ╔═╡ d1602be3-2724-4441-8dc3-8fa7713ed249
md"Very boring! But we have succeeded in making a simple environment with a temperature of 298K (or 25 degrees C) with 200 mm of water in each square.

Let's try something slightly more adventurous:"

# ╔═╡ e1ca4d60-383f-4c2a-945b-c7665b51bff9
begin
	# A temperature gradient spanning 10 degrees either side of total temperature from above. We can also give it a rate over which to change over time
	temp_change_rate = 0.2K/month
	temp_grad_env = tempgradAE(totalT - 10.0K, totalT + 10.0K, grid, totalW, 
		area_size, temp_change_rate)

	# Let's plot it to see what it looks like now
	heatmap(temp_grad_env.habitat.matrix' ./ K, clim = (278, 308), title = "Habitat")
end

# ╔═╡ 8f7b85a1-1e64-4d1a-90f6-55c6307a03cf
md"That's better! How about something even fancier?"

# ╔═╡ 01a11afc-7c9b-41f5-b308-003303dfa72a
begin
	# A temperature peak spanning 10 degrees either side of total temperature from above. We can also give it a rate over which to change over time
	temp_peak_env = peakedgradAE(totalT - 10.0K, totalT + 10.0K, grid, totalW, 
		area_size, temp_change_rate)

	# Let's plot it to see what it looks like now
	heatmap(temp_peak_env.habitat.matrix' ./ K, clim = (278, 308), title = "Habitat")
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
	death_rates = 0.15/ year
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
	water_req=(100.0nm/m^2);
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
	p3 = plot(288:0.1:308, pdf(Normal(optima[1]/K, niche_width[1]/K), 288:0.1:308), 
		label = "", xlabel = "Temperature (K)", ylabel = "Match to environment")
	for j in 2:numSpp
		plot!(288:0.1:308, pdf(Normal(optima[j]/K, niche_width[j]/K), 288:0.1:308), 
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
	simulation_time = 10years; time_step = 1month; record_interval = 1month
	len_sim = length(0month:time_step:simulation_time)
	record_abuns = zeros(Int64, numSpp, prod(grid), len_sim)
	simulate_record!(record_abuns, example_eco, simulation_time, time_step, record_interval)
end

# ╔═╡ 0d0e0197-d826-4f62-943c-79783e5fa701
md"# _That's all folks!_

There are more examples that can be found in the examples folder. Be aware to check first how big the simulation is before running it on your laptop!
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Diversity = "d3d5718d-52de-57ab-b67a-eca7fd6175a4"
EcoSISTEM = "ed2dc23b-ada4-5fdb-a26f-56368a14ad8f"
OnlineStats = "a15396b6-48d5-5d58-9928-6d29437db91e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
SpatialEcology = "348f2d5d-71a3-5ad4-b565-8af070f99681"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
Distributions = "~0.24.18"
Diversity = "~0.5.6"
EcoSISTEM = "~0.1.2"
OnlineStats = "~1.5.11"
Plots = "~0.28.4"
SpatialEcology = "~0.9.7"
Unitful = "~1.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "345a14764e43fe927d6f5c250fe4c8e4664e6ee8"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "2.4.0"

[[ArchGDAL]]
deps = ["Dates", "DiskArrays", "GDAL", "GeoFormatTypes", "GeoInterface", "Tables"]
git-tree-sha1 = "ac9a3bdf0b1cc5bd276c530b645194bf57373ecf"
uuid = "c9ce4bd3-c3d5-55b8-8973-c0e20141b8c3"
version = "0.6.0"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[ArrayInterface]]
deps = ["IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "be3671b34caec1d28a7915ca59cf8ba5a89a34fb"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.20"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "a4d07a1c313392a77042855df46c5f534076fab9"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.0"

[[AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "d127d5e4d86c7680b20c35d40b503c74b9a39b5e"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.4"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

[[Blosc]]
deps = ["Blosc_jll"]
git-tree-sha1 = "84cf7d0f8fd46ca6f1b3e0305b4b4a37afe50fd6"
uuid = "a74b3585-a348-5f62-a45c-50e91977d574"
version = "0.7.0"

[[Blosc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Lz4_jll", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "e747dac84f39c62aff6956651ec359686490134e"
uuid = "0b7ba130-8d10-5ba8-a3d6-c5182647fed9"
version = "1.21.0+0"

[[CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[Calculus]]
deps = ["Compat"]
git-tree-sha1 = "f60954495a7afcee4136f78d1d60350abd37a409"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.4.1"

[[CategoricalArrays]]
deps = ["DataAPI", "Future", "JSON", "Missings", "Printf", "Statistics", "StructTypes", "Unicode"]
git-tree-sha1 = "2ac27f59196a68070e132b25713f9a5bbc5fa0d2"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.8.3"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f53ca8d41e4753c41cdafa6ec5f7ce914b34be54"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "0.10.13"

[[CodecZlib]]
deps = ["BinaryProvider", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "05916673a2627dd91b4969ff8ba6941bc85a960e"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.6.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b9de8dc6106e09c79f3f776c27c62360d30e5eb8"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.9.1"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "InteractiveUtils", "Printf", "Reexport"]
git-tree-sha1 = "177d8b959d3c103a6d57574c38ee79c81059c31b"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.11.2"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b0b7e8a0d054fada22b64095b46469627a138943"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "2.2.1"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[Dagger]]
deps = ["Distributed", "LinearAlgebra", "MemPool", "Profile", "Random", "Serialization", "SharedArrays", "SparseArrays", "Statistics", "StatsBase", "Test"]
git-tree-sha1 = "39dd9351c6637a71d1925e26e49adcd8a25ebf0b"
uuid = "d58978e5-989f-55fb-8d15-ea34adc7bf54"
version = "0.8.0"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataFrames]]
deps = ["CategoricalArrays", "Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "Missings", "PooledArrays", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ecd850f3d2b815431104252575e7307256121548"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "0.21.8"

[[DataFramesMeta]]
deps = ["DataFrames", "Reexport"]
git-tree-sha1 = "5dee98ae14cca827ead9d7ede7210f5eb8b34d42"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.6.1"

[[DataStructures]]
deps = ["InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "88d48e133e6d3dd68183309877eac74393daa7eb"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.17.20"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffEqDiffTools]]
deps = ["LinearAlgebra", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "2d4f49c1839c1f30e4820400d8c109c6b16e869a"
uuid = "01453d9d-ee7c-5054-8395-0335cb756afa"
version = "0.13.0"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "85d2d9e2524da988bffaf2a381864e20d2dae08d"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.2.1"

[[DiskArrays]]
git-tree-sha1 = "599dc32bae654fa78056b15fed9b2af36f04ee44"
uuid = "3c3547ce-8d99-4f5e-a174-61eb10b00ae3"
version = "0.2.11"

[[Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "abe4ad222b26af3337262b8afb28fab8d215e9f8"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.3"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "a837fdf80f333415b69684ba8e8ae6ba76de6aaa"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.24.18"

[[Diversity]]
deps = ["AxisArrays", "DataFrames", "EcoBase", "InteractiveUtils", "LinearAlgebra", "Missings", "RecipesBase", "Requires", "Statistics"]
git-tree-sha1 = "819eb004d95692407c7ab5f4638607a998c50118"
uuid = "d3d5718d-52de-57ab-b67a-eca7fd6175a4"
version = "0.5.6"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[DoubleFloats]]
deps = ["GenericLinearAlgebra", "LinearAlgebra", "Polynomials", "Printf", "Quadmath", "Random", "Requires", "SpecialFunctions"]
git-tree-sha1 = "1c962cf7e75c09a5f1fbf504df7d6a06447a1129"
uuid = "497a8b3b-efae-58df-a0af-a86822472b78"
version = "1.1.23"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EcoBase]]
deps = ["RecipesBase"]
git-tree-sha1 = "4cb846e6dd521a122d6be123ab77375f563581ea"
uuid = "a58aae7d-b440-5a11-b283-399458f99aac"
version = "0.1.3"

[[EcoSISTEM]]
deps = ["ArchGDAL", "AxisArrays", "Calculus", "DataFrames", "Dates", "Distributions", "Diversity", "EcoBase", "Feather", "HCubature", "IndexedTables", "Interpolations", "JLD", "JuliaDB", "LibGit2", "Libdl", "LinearAlgebra", "Logging", "MPI", "Markdown", "Measures", "Missings", "NetCDF", "OnlineStats", "Optim", "Phylo", "Plots", "Printf", "RCall", "REPL", "Random", "RecipesBase", "Requires", "SHA", "SpatialEcology", "Statistics", "StatsBase", "UUIDs", "Unitful"]
git-tree-sha1 = "5dbab521599ed59e6db8b0dbf1ce8cd5a7c56b28"
uuid = "ed2dc23b-ada4-5fdb-a26f-56368a14ad8f"
version = "0.1.2"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "8041575f021cba5a099a456b4163c9a08b566a02"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.1.0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFMPEG]]
deps = ["BinaryProvider", "Libdl"]
git-tree-sha1 = "9143266ba77d3313a4cf61d8333a1970e8c5d8b6"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.2.4"

[[Feather]]
deps = ["CategoricalArrays", "DataFrames", "Dates", "FlatBuffers", "Mmap", "Tables"]
git-tree-sha1 = "5a9c74bfef934bcd813570a513ccf38f25555352"
uuid = "becb17da-46f6-5d3c-ad1b-1c5fe96bc73c"
version = "0.5.9"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "256d8e6188f3f1ebfa1a5d17e072a0efafa8c5bf"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.10.1"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "693210145367e7685d8604aee33d9bfb85db8b31"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.11.9"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[FixedPointNumbers]]
git-tree-sha1 = "d14a6fa5890ea3a7e5dcab6811114f132fec2b4b"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.6.1"

[[FlatBuffers]]
deps = ["Parameters", "Pkg", "Test"]
git-tree-sha1 = "dd5c2460639eef7962178216fd150c60fa42a1c4"
uuid = "53afe959-3a16-52fa-a8da-cf864710bae9"
version = "0.5.3"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "b5e930ac60b613ef3406da6d4f42c35d8dc51419"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.19"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GDAL]]
deps = ["CEnum", "GDAL_jll", "NetworkOptions", "PROJ_jll"]
git-tree-sha1 = "899e159dab7953918c3615290743c401121f02f9"
uuid = "add2ef01-049f-52c4-9ee2-e494f65e021a"
version = "1.2.2"

[[GDAL_jll]]
deps = ["Artifacts", "Expat_jll", "GEOS_jll", "JLLWrappers", "LibCURL_jll", "LibSSH2_jll", "Libdl", "Libtiff_jll", "MbedTLS_jll", "OpenJpeg_jll", "PROJ_jll", "Pkg", "SQLite_jll", "Zlib_jll", "Zstd_jll", "libgeotiff_jll", "nghttp2_jll"]
git-tree-sha1 = "439c33eb4dfa74a43a1e96b1d758aeb3cbc33dc3"
uuid = "a7073274-a066-55f0-b90d-d619367d196c"
version = "300.202.100+0"

[[GEOS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "45d0ddfd29620ac9b2d1072801e90fb016c5f94c"
uuid = "d604d12d-fa86-5845-992e-78dc15976526"
version = "3.9.0+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "LinearAlgebra", "Printf", "Random", "Serialization", "Sockets", "Test"]
git-tree-sha1 = "c690c2ab22ac9ee323d9966deae61a089362b25c"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.44.0"

[[GenericLinearAlgebra]]
deps = ["LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "ff291c1827030ffaacaf53e3c83ed92d4d5e6fb6"
uuid = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
version = "0.2.5"

[[GeoFormatTypes]]
git-tree-sha1 = "bb75ce99c9d6fb2edd8ef8ee474991cdacf12221"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.3.0"

[[GeoInterface]]
deps = ["RecipesBase"]
git-tree-sha1 = "38a649e6a52d1bea9844b382343630ac754c931c"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "0.5.5"

[[GeometryTypes]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "78f0ce9d01993b637a8f28d84537d75dc0ce8eef"
uuid = "4d00f742-c7ba-57c2-abde-4428a4b178cb"
version = "0.7.10"

[[Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HCubature]]
deps = ["Combinatorics", "DataStructures", "LinearAlgebra", "QuadGK", "StaticArrays"]
git-tree-sha1 = "134af3b940d1ca25b19bc9740948157cee7ff8fa"
uuid = "19dc6840-f33b-545b-b366-655c7e3ffd49"
version = "1.5.0"

[[HDF5]]
deps = ["Blosc", "HDF5_jll", "Libdl", "Mmap", "Random"]
git-tree-sha1 = "0713cbabdf855852dfab3ce6447c87145f3d9ea8"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.13.6"

[[HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fd83fa0bde42e01952757f01149dd968c06c4dba"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.0+1"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[IndexedTables]]
deps = ["DataAPI", "DataValues", "Distributed", "IteratorInterfaceExtensions", "OnlineStatsBase", "PooledArrays", "SparseArrays", "Statistics", "StructArrays", "TableTraits", "TableTraitsUtils", "Tables", "WeakRefStrings"]
git-tree-sha1 = "e41ee5688e404b49795a85dcb1da2dafb4409645"
uuid = "6deec6e2-d858-57c5-ab9b-e6ca5bd20e43"
version = "0.12.6"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "1470c80592cf1f0a35566ee5e93c5f8221ebc33a"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.3"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[InvertedIndices]]
deps = ["Test"]
git-tree-sha1 = "15732c475062348b0165684ffe28e85ea8396afc"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.0.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IterableTables]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Requires", "TableTraits", "TableTraitsUtils"]
git-tree-sha1 = "70300b876b2cebde43ebc0df42bc8c94a144e1b4"
uuid = "1c8ee90f-4401-5389-894e-7a04a3dc0f4d"
version = "1.0.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLD]]
deps = ["FileIO", "HDF5", "Printf"]
git-tree-sha1 = "f6cf928214ae7c0e7550b2424a57f11875d7e49a"
uuid = "4138dd39-2aa7-5051-a626-17a0bb65d9c8"
version = "0.10.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[JuliaDB]]
deps = ["Dagger", "DataValues", "Distributed", "Glob", "IndexedTables", "MemPool", "Nullables", "OnlineStats", "OnlineStatsBase", "PooledArrays", "Printf", "Random", "RecipesBase", "Serialization", "Statistics", "StatsBase", "TextParse", "WeakRefStrings"]
git-tree-sha1 = "a19da3e8b971e65d4859ee13bfae8a06df5681f4"
uuid = "a93385a2-3734-596a-9a66-3cfbb77141e6"
version = "0.13.0"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "7bd5f6565d80b6bf753738d2bc40a5dfea072070"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.2.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[MPI]]
deps = ["Distributed", "DocStringExtensions", "Libdl", "MPICH_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "Pkg", "Random", "Requires", "Serialization", "Sockets"]
git-tree-sha1 = "6e8c30afdcbb6167cf5d470b6333f4db01cc366f"
uuid = "da04e1cc-30fd-572f-bb4f-1f8673147195"
version = "0.17.2"

[[MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6cafe3f9747c0a0740611e2dffc4d37248fb691"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "3.4.2+0"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[MemPool]]
deps = ["DataStructures", "Distributed", "Mmap", "Random", "Serialization", "Sockets", "Test"]
git-tree-sha1 = "d52799152697059353a8eac1000d32ba8d92aa25"
uuid = "f9f48841-c794-520a-933b-121f7ba6ed94"
version = "0.2.0"

[[MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5c90234b3967684c9c6f87b4a54549b4ce21836"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.3+0"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f8c673ccc215eb50fcadb285f522420e29e69e1c"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "0.4.5"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "748f6e1e4de814b101911e64cc12d83a6af66782"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.2"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "144bab5b1443545bc4e791536c9f1eacb4eed06a"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.1"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetCDF]]
deps = ["DiskArrays", "Formatting", "NetCDF_jll"]
git-tree-sha1 = "23b0e32fde256a4e2e497e678abcf956ed26204b"
uuid = "30363a11-5582-574a-97bb-aa9a979735b9"
version = "0.11.3"

[[NetCDF_jll]]
deps = ["Artifacts", "HDF5_jll", "JLLWrappers", "LibCURL_jll", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Pkg", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "0cf4d1bf2ef45156aed85c9ac5f8c7e697d9288c"
uuid = "7243133f-43d8-5620-bbf4-c2c921802cf3"
version = "400.702.400+0"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Nullables]]
deps = ["Compat"]
git-tree-sha1 = "ae1a63457e14554df2159b0b028f48536125092d"
uuid = "4d1e1d77-625e-5b40-9113-a560ec7a8ecd"
version = "0.0.8"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "4f825c6da64aebaa22cc058ecfceed1ab9af1c7e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.3"

[[OnlineStats]]
deps = ["AbstractTrees", "Dates", "LinearAlgebra", "OnlineStatsBase", "OrderedCollections", "Random", "RecipesBase", "Statistics", "StatsBase", "SweepOperator"]
git-tree-sha1 = "e8ac92503ddcbe9ea586bfbb4c8e8c6e2c14d5cc"
uuid = "a15396b6-48d5-5d58-9928-6d29437db91e"
version = "1.5.11"

[[OnlineStatsBase]]
deps = ["AbstractTrees", "Dates", "LinearAlgebra", "OrderedCollections", "Statistics", "StatsBase"]
git-tree-sha1 = "8efa5acf7af1623eabeebdc82ef54396adf16f71"
uuid = "925886fa-5bf2-5e8e-b522-a9147a512338"
version = "1.4.5"

[[OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "a784e5133fc7e204c900f2cf38ed37a92ff9248d"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "4.1.1+2"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Optim]]
deps = ["Calculus", "DiffEqDiffTools", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "Random", "SparseArrays", "StatsBase", "Test"]
git-tree-sha1 = "a626e09c1f7f019b8f3a30a8172c7b82d2f4810b"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "0.18.1"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[PROJ_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "LibSSH2_jll", "Libdl", "Libtiff_jll", "MbedTLS_jll", "Pkg", "SQLite_jll", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "2435e91710d7f97f53ef7a4872bf1f948dc8e5f8"
uuid = "58948b4f-47e0-5654-a9ad-f609743f8632"
version = "700.202.100+0"

[[Parameters]]
deps = ["Markdown", "OrderedCollections", "REPL", "Test"]
git-tree-sha1 = "70bdbfb2bceabb15345c0b54be4544813b3444e4"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.10.3"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "94bf17e83a0e4b20c8d77f6af8ffe8cc3b386c0a"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.1"

[[Phylo]]
deps = ["AxisArrays", "DataFrames", "DataStructures", "Distributions", "IterableTables", "LightGraphs", "Missings", "Printf", "Random", "RecipesBase", "Requires", "SimpleTraits", "Statistics", "Tokenize", "Unitful"]
git-tree-sha1 = "c3ee113b4696d033d83c480de66c6daf6a9b6050"
uuid = "aea672f4-3940-5932-aa44-993d1c3ff149"
version = "0.4.19"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "87a4ea7f8c350d87d3a8ca9052663b633c0b2722"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "1.0.3"

[[PlotUtils]]
deps = ["Colors", "Dates", "Printf", "Random", "Reexport"]
git-tree-sha1 = "51e742162c97d35f714f9611619db6975e19384b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "0.6.5"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "FFMPEG", "FixedPointNumbers", "GR", "GeometryTypes", "JSON", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "Reexport", "Requires", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "efbe466a790d7e8a5c4b5ee1601c0c8edc99780b"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "0.28.4"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "0bbfdcd8cda81b8144de4be8a67f5717e959a005"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.14"

[[PooledArrays]]
deps = ["DataAPI"]
git-tree-sha1 = "b1333d4eced1826e15adbdf01a4ecaccca9d353c"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "0.5.3"

[[PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[Quadmath]]
deps = ["Printf", "Random", "Requires"]
git-tree-sha1 = "5a8f74af8eae654086a1d058b4ec94ff192e3de0"
uuid = "be4d8f0f-7fa4-5f49-b795-2f01399ab2dd"
version = "0.5.5"

[[RCall]]
deps = ["CategoricalArrays", "Conda", "DataFrames", "DataStructures", "Dates", "Libdl", "Missings", "REPL", "Random", "Requires", "StatsModels", "WinReg"]
git-tree-sha1 = "80a056277142a340e646beea0e213f9aecb99caa"
uuid = "6f49c342-dc21-5d91-9882-a32aef131414"
version = "0.13.12"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RandomBooleanMatrices]]
deps = ["Random", "RandomNumbers", "SparseArrays", "StatsBase", "Test"]
git-tree-sha1 = "c13963a90c579b3bc241d858e4d063d778deeef5"
uuid = "9ae346a0-3d16-5633-ad70-ddb60ab77eac"
version = "0.1.1"

[[RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "a752043df7488ca8bcbe05fa82c831b5e2c67211"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.2"

[[RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[Ratios]]
git-tree-sha1 = "37d210f612d70f3f7d57d488cb3b6eff56ad4e41"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.0"

[[RecipesBase]]
git-tree-sha1 = "7bdce29bc9b2f5660a6e5e64d64d91ec941f6aa2"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "0.7.0"

[[Reexport]]
deps = ["Pkg"]
git-tree-sha1 = "7b1d07f411bc8ddb7977ec7f377b97b158514fe0"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "0.2.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SQLite_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "9a0e24b81e3ce02c4b2eb855476467c7b93b8a8f"
uuid = "76ed43ae-9a5d-5a62-8c75-30186b810ce8"
version = "3.36.0+0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[ShiftedArrays]]
git-tree-sha1 = "22395afdcf37d6709a5a0766cc4a5ca52cb85ea0"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "1.0.0"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "ee010d8f103468309b8afac4abb9be2e18ff1182"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "0.3.2"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures", "Random", "Test"]
git-tree-sha1 = "03f5898c9959f8115e30bc7226ada7d0df554ddd"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "0.3.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpatialEcology]]
deps = ["DataFrames", "DataFramesMeta", "Distances", "EcoBase", "Random", "RandomBooleanMatrices", "RandomNumbers", "RecipesBase", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "0210cdedd3276e89a3d4a6ade1e9aeb91df06a06"
uuid = "348f2d5d-71a3-5ad4-b565-8af070f99681"
version = "0.9.7"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "508822dca004bf62e210609148511ad03ce8f1d8"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.0"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "62701892d172a2fa41a1f829f66d2b0db94a9a63"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "da4cf579416c81994afd6322365d00916c79b8ae"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "0.12.5"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics"]
git-tree-sha1 = "19bfcb46245f69ff4013b3df3b977a289852c3a1"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.32.2"

[[StatsFuns]]
deps = ["LogExpFunctions", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "30cd8c360c54081f806b1ee14d2eecbef3c04c49"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.8"

[[StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "3db41a7e4ae7106a6bcff8aa41833a4567c04655"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.21"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "Tables"]
git-tree-sha1 = "8099ed9fb90b6e754d6ba8c6ed8670f010eadca0"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.4.4"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "e36adc471280e8b346ea24c5c87ba0571204be7a"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.7.2"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[SweepOperator]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "20da2784e79fc0bfdd70592d1b47d7a6034e82d1"
uuid = "7522ee7d-7047-56d0-94d9-4bc626e7058d"
version = "0.3.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[TableTraitsUtils]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Missings", "TableTraits"]
git-tree-sha1 = "8fc12ae66deac83e44454e61b02c37b326493233"
uuid = "382cd787-c1b6-5bf2-a167-d5b971a19bda"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "d0c690d37c73aeb5ca063056283fde5585a41710"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TextParse]]
deps = ["CodecZlib", "DataStructures", "Dates", "DoubleFloats", "Mmap", "Nullables", "PooledArrays", "WeakRefStrings"]
git-tree-sha1 = "26b43d6746b52cca13c4cdef90f89652273b413e"
uuid = "e0df1984-e451-5cb5-8b61-797a481e67e3"
version = "0.9.1"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "81753f400872e5074768c9a77d4c44e70d409ef0"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.6"

[[Tokenize]]
git-tree-sha1 = "eee92eda3cc8e104b7e56ff4c1fcf0d78ca37c89"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.18"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "7c53c35547de1c5b9d46a4797cf6d8253807108c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.5"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a981a8ef8714cba2fd9780b22fd7a469e7aaf56d"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.9.0"

[[VersionParsing]]
git-tree-sha1 = "80229be1f670524750d905f8fc8148e5a8c4537f"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.0"

[[WeakRefStrings]]
deps = ["DataAPI", "Random", "Test"]
git-tree-sha1 = "28807f85197eaad3cbd2330386fac1dcb9e7e11d"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "0.6.2"

[[WinReg]]
deps = ["Test"]
git-tree-sha1 = "808380e0a0483e134081cc54150be4177959b5f4"
uuid = "1b915085-20d7-51cf-bf83-8f477d6f5128"
version = "0.3.1"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "59e2ad8fd1591ea019a5259bd012d7aee15f995c"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libgeotiff_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "PROJ_jll", "Pkg"]
git-tree-sha1 = "a5cc2e3dd7b1c1e783a61b8ab7de03eebddfed60"
uuid = "06c338fa-64ff-565b-ac2f-249532af990e"
version = "1.6.0+1"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─6b280ae0-6af9-4d3c-91c8-f02a20d28d18
# ╟─496434ae-1bc8-4472-8981-2d9fc08b688a
# ╟─39be3104-92c0-40fe-ad10-905d6c889917
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
