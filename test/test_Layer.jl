# SPDX-License-Identifier: LGPL-3.0-or-later

module TestLayer

using EcoSISTEM
using Diversity
# `materialise` is `public` rather than exported, so it has to be named.
using EcoSISTEM: materialise
using Unitful
using Unitful.DefaultSymbols
using Test
using EcoBase
using Plots
using EcoSISTEM.Units
using RasterDataSources
using DimensionalData
using Random
using DimensionalData: DimArray, X, Y, Ti

@testset "Habitats" begin
    grid = (10, 10)
    area = 25.0km^2
    totalK = 10000.0kJ / km^2 / day
    numNiches = 4
    active = fill(true, grid)

    # Grid decided first: 25 km^2 over 10 × 10 is 0.5 km cells - the same
    # `sqrt(area / prod(dimension))` the deprecated builders computed inside, which is what the
    # `xcellsize`/`ycellsize` assertions below still check against - now on the grid rather than
    # on the layer.
    studyarea = StudyArea(extent = (sqrt(area), sqrt(area)),
                          cellsize = sqrt(area) / grid[1], verbosity = :silent)
    supply = UniformSpec(totalK, axis = SolarRadiation)

    fillval = 0.0
    habitat = GridHabitat(regime = UniformSpec(fillval,
                                               axis = EcoSISTEM.NicheAxis),
                          supply = supply, area = studyarea)
    @test EcoSISTEM.iscontinuous(habitat.regime) == true
    @test eltype(habitat.regime) == typeof(habitat.regime.matrix[1])
    @test size(habitat.regime, 1) == grid[1]
    # **These moved off the layer and onto the grid**, and three of them changed answer in the
    # process. A layer holds values *on* a grid and has no origin or CRS of its own, so answering
    # from one means inventing an origin and stripping the cell size to a bare number.
    # See `test_EcoBaseInterface.jl` for the geographic case, which a layer could not express.
    gridof = EcoBase.getcoords(habitat)
    @test EcoBase.xmin(gridof) == 0.0km      # this grid really does start at the origin...
    @test EcoBase.ymin(gridof) == 0.0km
    @test EcoBase.xcellsize(gridof) == sqrt(area / prod(grid))   # ...and is unitful now
    @test EcoBase.ycellsize(gridof) == sqrt(area / prod(grid))
    # dimension 1 of regime.matrix is y (rows), dimension 2 is x (columns)
    @test EcoBase.xcells(gridof) == size(habitat.regime, 2)
    @test EcoBase.ycells(gridof) == size(habitat.regime, 1)
    # `coordinates` gives real coordinates, never the *indices* verbatim.
    @test EcoBase.coordinates(gridof) != EcoBase.indices(gridof)
    @test EcoBase.coordinates(gridof)[1, :] == [0.0km, 0.0km]
    # And column 1 is **x**, which is EcoBase's order and the opposite of this package's - so the
    # x index is the one that repeats in blocks, since the package's cell order runs y fastest.
    @test EcoBase.indices(gridof, 1) ==
          repeat(collect(1:grid[2]), inner = grid[1])
    @test EcoBase.indices(gridof, 2) ==
          repeat(collect(1:grid[1]), outer = grid[2])

    # A temperature gradient that warms: `GradientSpec` is the shape, `Varying(..., IncrementBy)`
    # the rate - the two stated separately rather than bundled into one builder's argument list.
    habitat = GridHabitat(regime = Varying(GradientSpec(-10.0K, 10.0K,
                                                        axis = Temperature),
                                           IncrementBy(0.01K /
                                                       month_mean_duration)),
                          supply = supply, area = studyarea)
    @test EcoSISTEM.iscontinuous(habitat.regime) == true
    @test eltype(habitat.regime) == typeof(habitat.regime.matrix[1])
    @test size(habitat.regime, 1) == grid[1]

    # A categorical niche regime.
    habitat = GridHabitat(regime = NicheSpec(numNiches,
                                             axis = EcoSISTEM.TypologyAxis),
                          supply = supply, area = studyarea)
    @test EcoSISTEM.iscontinuous(habitat.regime) == false
    @test eltype(habitat.regime) == typeof(habitat.regime.matrix[1])
    @test size(habitat.regime, 1) == grid[1]

    # Two more continuous regimes, built only to feed the `LayerCollection` assertions below - a
    # categorical member beside two continuous ones, which is the only such collection in the suite.
    # Their own per-layer properties are already asserted three times each above, so none of it is
    # repeated here.
    # Both are built on `studyarea`, so the members genuinely share a grid. The fixtures they
    # replace did **not**: one was on a degree grid and the other on a metre grid, so the
    # `_uniformcellside`/`getgridarea` assertions were reading one member's numbers for a collection that
    # could not have been built by any public route.
    temperature = GridHabitat(regime = UniformSpec(1.0K,
                                                   axis = Temperature),
                              supply = supply, area = studyarea)
    rainfall = GridHabitat(regime = UniformSpec(1.0mm / day,
                                                axis = Precipitation),
                           supply = supply, area = studyarea)

    # Test multi regimes
    regime = LayerCollection((temperature.regime, rainfall.regime))
    # **A collection has no `iscontinuous` and no `eltype` of its own** - its members may differ,
    # so both are asked of each member. That is what keeps a leaf and a collection-of-one from
    # answering differently, and it returns a `Tuple` rather than an allocating `Vector`.
    @test map(EcoSISTEM.iscontinuous, values(regime)) == (true, true)
    # **Built from a plain `Tuple`, yet named by AXIS** - the two members are on distinguishable
    # axes, so the collection derives `(:Temperature, :Precipitation)` rather than the old
    # `(:one, :two)`. A repeated axis, or a member with none, falls back to the positional names.
    @test keys(regime) == (:Temperature, :Precipitation)
    @test map(eltype, values(regime)) == (typeof(regime.Temperature.matrix[1]),
           typeof(regime.Precipitation.matrix[1]))
    # **The collection itself has no `eltype` method**, deliberately - see `src/collections.jl`.
    # That does **not** give a `MethodError`: `Base.eltype` has a universal fallback, so an
    # undefined type answers `Any`. What matters is that it is not the *wrong and specific*
    # answer - a generic method on `AbstractSpeciesRequirement` would hand back the whole backing
    # `NamedTuple` type, which is why that one had to move down onto the leaves.
    @test eltype(regime) === Any
    @test size(regime, 1) == grid[1]
    @test EcoSISTEM._uniformcellside(regime) == regime.Temperature.size
    @test isapprox(EcoSISTEM.getgridarea(regime),
                   regime.Temperature.size^2 * prod(grid),
                   rtol = 1e-5)
    @test EcoSISTEM.getgridshape(regime) == grid
    @test countsubcommunities(regime) == prod(grid)

    regime = LayerCollection((habitat.regime, temperature.regime,
                              rainfall.regime))
    @test map(EcoSISTEM.iscontinuous, values(regime)) == (false, true, true)
    # Named by axis here too, and the **categorical** member is what makes this worth asserting
    # separately: it is on `TypologyAxis`, so all three are distinguishable. It read
    # `(:one, :two, :three)` until `NicheSpec` was required to declare an axis - the positional
    # fallback was firing because the categorical member had none.
    @test keys(regime) ==
          (:TypologyAxis, :Temperature, :Precipitation)
    @test map(eltype, values(regime)) ==
          (typeof(regime.TypologyAxis.matrix[1]),
           typeof(regime.Temperature.matrix[1]),
           typeof(regime.Precipitation.matrix[1]))
    @test size(regime, 1) == grid[1]
    @test EcoSISTEM._uniformcellside(regime) == regime.TypologyAxis.size
    @test isapprox(EcoSISTEM.getgridarea(regime),
                   regime.TypologyAxis.size^2 * prod(grid),
                   rtol = 1e-5)
    @test EcoSISTEM.getgridshape(regime) == grid
    @test countsubcommunities(regime) == prod(grid)
end

@testset "type-level rejection of a wrongly-ordered regime/supply DimArray" begin
    # `ContinuousRegime`/`CategoricalRegime`/`Supply{A}` require their `matrix` to be a
    # `(Y, X)`-ordered `DimArray` - this is what actually prevents the historical X/Y mixup bug
    # from recurring (see `~/.claude/plans/fix-dimension-order-mismatch.md`, Step 2): a
    # wrongly-ordered array is rejected outright at construction (`MethodError`, not a runtime
    # check), never becoming a valid regime/supply in the first place.
    change = NoLayerChange()
    correct = DimArray(fill(274.0K, 3, 4), (Y((1:3)km), X((1:4)km)))
    wrong = DimArray(fill(274.0K, 3, 4), (X((1:3)km), Y((1:4)km)))

    @test ContinuousRegime(correct, change) isa ContinuousRegime
    @test_throws MethodError ContinuousRegime(wrong, change)

    dcorrect = DimArray(fill(1, 3, 4), (Y((1:3)km), X((1:4)km)))
    dwrong = DimArray(fill(1, 3, 4), (X((1:3)km), Y((1:4)km)))
    @test CategoricalRegime(dcorrect, change) isa CategoricalRegime
    @test_throws MethodError CategoricalRegime(dwrong, change)

    scorrect = DimArray(fill(1.0kJ / day, 3, 4), (Y((1:3)km), X((1:4)km)))
    swrong = DimArray(fill(1.0kJ / day, 3, 4), (X((1:3)km), Y((1:4)km)))
    @test Supply{SolarRadiation}(scorrect) isa Supply{SolarRadiation}
    @test_throws MethodError Supply{SolarRadiation}(swrong)
end

# **`materialise` returns a layer, so a layer has to plot** - otherwise making the array type
# private takes inspection with it, and inspection is what `materialise` is *for*.
#
# The regime recipes live in `src/Layer.jl`, not here: `RecipesBase` is a hard dependency, so
# the package draws its own layers and only the `ClimateRaster`/`ERA` recipes above need `Plots`.
# They are exercised from this file because it is the one that has `Plots` to render with.
@testset "a materialised layer plots" begin
    area = StudyArea(extent = (90.0km, 90.0km), cellsize = 10.0km,
                     verbosity = :silent)
    regime = materialise(GradientSpec(280.0K, 300.0K, axis = Temperature),
                         area)
    @test plot(regime).n == 1
    # **Titled from the layer's own axis, never from its unit.** The old recipes looked the value's
    # unit up in a four-entry `Dict`, which is the inference this subproject removes everywhere else.
    @test plot(regime).subplots[1][:title] == "Temperature (K)"
    # `(y, x)`: `matrix` is already in the orientation `heatmap` wants, so the plotted surface
    # must match it exactly - a transpose would pass a `size` check on a square grid and be wrong,
    # which is why the values are compared rather than the shape.
    @test plot(regime).series_list[1][:z].surf == ustrip.(parent(regime.matrix))

    # ...and a supply plots too, titled from *its* axis (`_resourcetitle`, `src/Demand.jl`).
    supply = materialise(UniformSpec(1.0kJ / (km^2 * day),
                                     axis = SolarRadiation),
                         area, role = EcoSISTEM.Resource)
    @test plot(supply).n == 1
    @test plot(supply).subplots[1][:title] == "Solar Radiation (kJ/day)"

    # **A categorical layer plots at all**, which it did not: its recipe read `unitdict` from a
    # *sibling recipe's* local scope, so this threw `UndefVarError` for every land-cover layer.
    niches = materialise(NicheSpec(4, axis = EcoSISTEM.TypologyAxis), area)
    @test plot(niches).n == 1
    @test plot(niches).subplots[1][:title] == "TypologyAxis"

    # ...and so does a layer whose unit was not one of the four the old `Dict` knew - a `KeyError`
    # before, and most axes are not one of those four.
    wind = materialise(UniformSpec(3.0m / s, axis = WindSpeed), area)
    @test plot(wind).subplots[1][:title] == "WindSpeed (m s⁻¹)"
end

# **The niche generator, on grids that are NOT square** - the shape that catches a y/x mix-up,
# and the only shape that does.
#
# **This is the shape of bug this file was otherwise blind to.** `_identifyclusters!` and
# `_fillin!` called `_getneighbours(M, y, x)` while its first argument has always been the *row*,
# left over from the package's x/y -> y/x switch. On a square grid the transposition is invisible -
# it silently reads the neighbours of the transposed cell - and every fixture in this file, and every
# `NicheSpec` fixture in the suite, was square: `(10, 10)`, `9 × 9`, `60km × 60km`. On any other
# shape it throws *"Coordinates outside grid"* from inside `_getneighbours`' own bounds check.
#
# So the sweep below is the test: **every** shape from 1 to 9 in each dimension, which is cheap
# (81 grids) and leaves the mistake nowhere to hide. Repeated, because the generator is
# stochastic - a percolation that happens to mark no cell exercises no neighbourhood at all.
@testset "the niche generator works on non-square grids" begin
    types = collect(1:3)
    weights = fill(1 / 3, 3)
    for ny in 1:9, nx in 1:9, _ in 1:3
        T = EcoSISTEM._nichefield((ny, nx), types, 0.5, weights)
        # The shape comes back as asked - `(y, x)`, rows then columns.
        @test size(T) == (ny, nx)
        # ...and every cell holds one of the classes asked for, none left unassigned.
        @test all(in(types), T)
    end

    # And through the builder that uses it, on a deliberately non-square grid. `extent` is `(y, x)`
    # like everything else, so 7 km by 11 km in 1 km cells is 7 rows and 11 columns -- which is
    # exactly the asymmetry under test, and is itself a place the order could be got wrong.
    habitat = GridHabitat(regime = NicheSpec(4, axis = EcoSISTEM.TypologyAxis),
                          supply = UniformSpec(10000.0kJ / km^2 / day,
                                               axis = SolarRadiation),
                          area = StudyArea(extent = (7.0km, 11.0km),
                                           cellsize = 1.0km,
                                           verbosity = :silent))
    @test size(habitat.regime.matrix) == (7, 11)
    @test EcoSISTEM.iscontinuous(habitat.regime) == false
    @test EcoBase.ycells(EcoBase.getcoords(habitat)) == 7
    @test EcoBase.xcells(EcoBase.getcoords(habitat)) == 11
end

# **The generator is DETERMINISTIC given a seed - and this is the assertion that has teeth.**
#
# The sweep above cannot catch what this catches. `_fillin!` used `isassigned(T, y, x)` on a
# `T = Array{Int64}(undef, ...)`, where it is **always `true`** (`isassigned` reports a filled
# *reference* slot, so it means nothing for a bits eltype). Every neighbour therefore counted,
# including never-written cells, and the tallies were taken over **uninitialised memory** - but the
# result is still `types[ind]`, so `all(in(types), T)` passed happily. Measured: the same seed
# gave **7 distinct maps in 12 runs**.
#
# **The churn between the two draws is the point, not decoration.** Garbage read from uninitialised
# memory comes from previously-freed pages, so a bare draw-twice can agree by luck; allocating and
# dropping a few megabytes in between makes the pages differ, which is exactly the condition under
# which such a defect shows itself.
# **The shapes and the clumpiness are CALIBRATED, not arbitrary** - detection is very sensitive to
# both. Measured against the pre-fix code, 10 seeds per cell: `clumpiness = 0.5` caught it **10/10**
# at every size, while `0.9` caught **0/10** on a small grid (a high clumpiness leaves almost no
# unlabelled cell for `_fillin!` to guess at, so the faulty branch barely runs). And it must be
# calibrated **in situ**: an earlier pair of shapes scored 10/10 in a fresh process but only **1/10**
# inside this file, where the 490-assertion sweep above churns the heap first.
# **A targeted poison was tried and is WORSE** - pre-filling the pool with same-sized `3`-valued
# matrices so the next `undef` would land on them scored **7/20** and **2/20** against generic
# churn's **10/20** and **9/20**. Measured, not assumed; do not "improve" it that way.
# **So the reliable lever is the number of draws, not a cleverer churn.** Detection per assertion
# is noisy (20-50% observed), so 30 of them is what makes a reintroduction essentially certain to be
# caught rather than merely likely; each is sub-millisecond, so the count is nearly free.
@testset "the niche generator is reproducible from a seed" begin
    types = collect(1:3)
    weights = fill(1 / 3, 3)
    for seed in 1:15, (ny, nx) in ((13, 17), (23, 29))
        Random.seed!(seed)
        first_pass = EcoSISTEM._nichefield((ny, nx), types, 0.5, weights)
        # Churn the heap so any read of uninitialised memory sees something different.
        for _ in 1:20
            sum(rand(Int64, 100_000))
        end
        Random.seed!(seed)
        @test EcoSISTEM._nichefield((ny, nx), types, 0.5, weights) == first_pass
    end
end

end
