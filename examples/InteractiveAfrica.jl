### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ 1e9bcd31-9146-494b-98f5-27569c54ab0f


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
	file = "Africa.tif"
	africa = readfile(file, -25.0°, 50.0°, -35.0°, 40.0°)
	active =  Array{Bool, 2}(.!isnan.(africa))
	# Set up initial parameters for ecosystem
	grd = size(africa); req= 10.0kJ; individuals=3*10^8; area = 64e6km^2; totalK = 1000.0kJ/km^2

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
	vars[end] = (1/niche_width) * K
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
	burnin = 1years; times = 10years; timestep = 1month; record_interval = 1month;
	lensim = length(0years:record_interval:times)
	abuns = zeros(numSpecies, prod(grd), lensim)
	
	# Run simulation for burnin and then add invasive species
	@time simulate!(eco, burnin, timestep)
	eco.abundances.grid[end, 50, 50] = 100
	@time simulate_record!(abuns, eco, times, timestep, record_interval);
	
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
	africa_temp = world.array[-25°.. 50°, -35° .. 40°, 1]
	plot(africa_temp)
end

# ╔═╡ ee925e21-b0b6-478e-a3a0-573e8497b9f6
begin
	temp = uconvert.(K, africa_temp .* °C)
	africa_new = Worldclim_bioclim(AxisArray(temp, AxisArrays.axes(africa_temp)))
	active_new =  Array{Bool, 2}(.!isnan.(africa))
	# Set up initial parameters for ecosystem
	grd_new = size(africa_new.array); req_new= 10.0kJ; individuals_new=3*10^8; area_new = 64e6km^2; totalK_new = 1000.0kJ/km^2

	# Set up how much energy each species consumes
	energy_vec_new = SolarRequirement(fill(req, numSpecies))


	# Set rates for birth and death
	birth_new = 0.6/year
	death_new = 0.6/year
	longevity_new = 1.0
	survival_new = 0.1
	boost_new = 1.0
	# Collect model parameters together
	param_new = EqualPop(birth, death, longevity, survival, boost)

	# Create kernel for movement
	kernel_new = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
	movement_new = AlwaysMovement(kernel, Torus())


	# Create species list, including their temperature preferences, seed abundance and native status
	opts_new = fill(meantemp, numSpecies)
	vars_new = fill(50.0K, numSpecies)
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
	burnin_new = 1years; times_new = 10years; timestep_new = 1year; record_interval_new = 1year;
	lensim_new = length(0years:record_interval_new:times_new)
	abuns_new = zeros(numSpecies, prod(grd_new), lensim_new)
	
	# Run simulation for burnin and then add invasive species
	@time simulate!(eco_new, burnin_new, timestep_new)
	@time simulate_record!(abuns_new, eco_new, times_new, timestep_new, record_interval_new);
	
end

# ╔═╡ 59cf0c00-8a0b-41de-b895-cc511b7ed07b
begin
	mean_abuns = reshape(mean(eco_new.abundances.matrix, dims = 1)[1, :], grd_new)
	mean_abuns[.!eco_new.abenv.active] .= NaN
	heatmap(mean_abuns)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AxisArrays = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Diversity = "d3d5718d-52de-57ab-b67a-eca7fd6175a4"
EcoSISTEM = "ed2dc23b-ada4-5fdb-a26f-56368a14ad8f"
OnlineStats = "a15396b6-48d5-5d58-9928-6d29437db91e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
RasterDataSources = "3cb90ccd-e1b6-4867-9617-4276c8b2ca36"
SpatialEcology = "348f2d5d-71a3-5ad4-b565-8af070f99681"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
AxisArrays = "~0.4.6"
Distributions = "~0.24.18"
Diversity = "~0.5.6"
EcoSISTEM = "~0.1.2"
OnlineStats = "~1.5.11"
Plots = "~0.28.4"
RasterDataSources = "~0.5.6"
SpatialEcology = "~0.9.7"
Unitful = "~1.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ASCIIrasters]]
git-tree-sha1 = "0cb0046798af8ac8561334c5a2a31f015e53c2b1"
uuid = "81770e7c-c736-4fa5-8129-46dd21831640"
version = "0.1.1"

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
version = "1.1.1"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[ArrayInterface]]
deps = ["IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "2e004e61f76874d153979effc832ae53b56c20ee"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.22"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "a4d07a1c313392a77042855df46c5f534076fab9"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.0"

[[AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

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
version = "0.5.2+0"

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

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[EcoBase]]
deps = ["RecipesBase"]
git-tree-sha1 = "4cb846e6dd521a122d6be123ab77375f563581ea"
uuid = "a58aae7d-b440-5a11-b283-399458f99aac"
version = "0.1.3"

[[EcoSISTEM]]
deps = ["ArchGDAL", "AxisArrays", "Calculus", "DataFrames", "Dates", "Distributions", "Diversity", "EcoBase", "Feather", "HCubature", "IndexedTables", "Interpolations", "JLD2", "LibGit2", "Libdl", "LinearAlgebra", "Logging", "MPI", "Markdown", "Measures", "Missings", "NetCDF", "OnlineStats", "Optim", "Phylo", "Plots", "Printf", "REPL", "Random", "RecipesBase", "Requires", "SHA", "SimpleSDMLayers", "SpatialEcology", "Statistics", "StatsBase", "UUIDs", "Unitful"]
git-tree-sha1 = "27ce83d3a837ab5f7a880014abed5a60bd9c9f74"
uuid = "ed2dc23b-ada4-5fdb-a26f-56368a14ad8f"
version = "0.1.3"

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

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

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

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HCubature]]
deps = ["Combinatorics", "DataStructures", "LinearAlgebra", "QuadGK", "StaticArrays"]
git-tree-sha1 = "134af3b940d1ca25b19bc9740948157cee7ff8fa"
uuid = "19dc6840-f33b-545b-b366-655c7e3ffd49"
version = "1.5.0"

[[HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fd83fa0bde42e01952757f01149dd968c06c4dba"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.0+1"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

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

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

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

[[JLD2]]
deps = ["CodecZlib", "DataStructures", "MacroTools", "Mmap", "Pkg", "Printf", "Requires", "UUIDs"]
git-tree-sha1 = "9f2f2f24e60305feb6ae293a617ddf06f429efc3"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.3"

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

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

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
deps = ["Libdl", "libblastrampoline_jll"]
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

[[MPI]]
deps = ["Distributed", "DocStringExtensions", "Libdl", "MPICH_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "Pkg", "Random", "Requires", "Serialization", "Sockets"]
git-tree-sha1 = "6e8c30afdcbb6167cf5d470b6333f4db01cc366f"
uuid = "da04e1cc-30fd-572f-bb4f-1f8673147195"
version = "0.17.2"

[[MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "480a8be51647c9dec19ff2211e26c9f30787840b"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "3.4.2+2"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "a16aa086d335ed7e0170c5265247db29172af2f9"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.3+2"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f8c673ccc215eb50fcadb285f522420e29e69e1c"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "0.4.5"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

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
version = "1.2.0"

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

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "08334c7cd6cb033274bc90918ece2174b5f84cc7"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "4.1.1+5"

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
version = "1.8.0"

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

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
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

[[RasterDataSources]]
deps = ["ASCIIrasters", "Dates", "DelimitedFiles", "HTTP", "JSON", "URIs", "ZipFile"]
git-tree-sha1 = "aeb786b2cc040bc320bc97630f91fc2c33612873"
uuid = "3cb90ccd-e1b6-4867-9617-4276c8b2ca36"
version = "0.5.6"

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
version = "0.7.0"

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

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "ee010d8f103468309b8afac4abb9be2e18ff1182"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "0.3.2"

[[SimpleSDMLayers]]
deps = ["ArchGDAL", "Downloads", "RecipesBase", "Requires", "Statistics", "ZipFile"]
git-tree-sha1 = "a5f2bd0c13ae434963aadb8d1e26bf3dba81d9cc"
uuid = "2c645270-77db-11e9-22c3-0f302a89c64c"
version = "0.7.0"

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
version = "1.0.0"

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
version = "1.10.0"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[Tokenize]]
git-tree-sha1 = "eee92eda3cc8e104b7e56ff4c1fcf0d78ca37c89"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.18"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "7c53c35547de1c5b9d46a4797cf6d8253807108c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.5"

[[URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

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

[[WeakRefStrings]]
deps = ["DataAPI", "Random", "Test"]
git-tree-sha1 = "28807f85197eaad3cbd2330386fac1dcb9e7e11d"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "0.6.2"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "59e2ad8fd1591ea019a5259bd012d7aee15f995c"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.3"

[[ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "c3a5637e27e914a7a445b8d0ad063d701931e9f7"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.9.3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

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
version = "1.48.0+0"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─39b75180-f384-11eb-3449-4f7c2ad25d99
# ╠═108951ec-3ecf-4f0c-b5e4-e79d00b1bfac
# ╠═5dd6e516-cfb5-4a4f-9f1c-50334c1f743f
# ╠═1e9bcd31-9146-494b-98f5-27569c54ab0f
# ╠═e9e9065b-cf7e-4d37-8162-88f0076ad1eb
# ╠═220e4af6-f228-4b8d-a77e-0ddbf5fc6705
# ╠═7e16f197-874b-482d-80b6-13a62ddda1f7
# ╠═03ade0dc-271d-407c-a0f1-583007b168d8
# ╠═ee925e21-b0b6-478e-a3a0-573e8497b9f6
# ╠═59cf0c00-8a0b-41de-b895-cc511b7ed07b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
