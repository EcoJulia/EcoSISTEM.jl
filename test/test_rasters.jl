# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests for `src/rasters.jl` - reading, reprojecting and masking rasters: the layer of machinery that
# knows about rasters and grids, but not about layers, areas or habitats.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_rasters.jl"])'
#
# This is where the *memory* goes: reading real rasters can take a worker to 13 GB resident and
# trip Rasters' own memory guard. Keeping these testsets together bounds that to one worker rather
# than mixing it through a file that is otherwise cheap synthetic work.

module TestRasters

using EcoSISTEM
using Statistics: mean
# `[C7-VIS]` C: these are `public` rather than exported - a spec is what a user writes,
# and these are what it materialises into.
using EcoSISTEM: SeriesLayerChange, AbsoluteChange
using EcoSISTEM: hasdata, landcoverclass
using EcoSISTEM: materialise
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using DimensionalData: DimensionalData, DimArray, X, Y, dims
using Rasters
using RasterDataSources
using ArchGDAL
using Extents: Extent
using Test

include("rasterfixtures.jl")
include("buildfixtures.jl")

# A layer's footprint is the ground it *covers*, so it must be measured at its cells' outer
# edges, never at their reference corners. `extrema` of the lookup values gives corners: with the
# `Start` locus GDAL files declare, the largest value is where the last cell *begins*, so a whole
# cell is discarded off the top and right of every layer - WorldClim reporting -90°...89.8333° for a
# grid that reaches 90°. It matters because `Touches`, which
# `_applycut` selects with, works on intervals, so the footprint and the read disagreed.
@testset "a layer's footprint reaches its cells' outer edges" begin
    # The fixtures use the `Intervals(Start)` locus a reader produces, so a coordinate is its
    # cell's **lower** edge and the extent runs one whole cell beyond the last value.
    r = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5))
    e = EcoSISTEM._dimsextent(r.array, nothing)
    @test e.Y == (0.0°, 5.0°)       # values are 0°...4°, cells 1° wide
    @test e.X == (0.0°, 5.0°)
    # ...and asked for `Center` explicitly the same five values describe different ground, which is
    # the half-cell difference that made the old default a misdescription of real data.
    c = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5),
                    locus = DimensionalData.Lookups.Center())
    @test EcoSISTEM._dimsextent(c.array, nothing).Y == (-0.5°, 4.5°)
    # ...and the span is a whole number of cells, which the value-extrema version never gave.
    @test e.Y[2] - e.Y[1] == 5.0° && e.X[2] - e.X[1] == 5.0°

    # A projected layer is measured the same way, in its own units.
    bng = _bngraster(WorldClim{BioClim}, fill(290.0K, 9, 9))
    by = EcoSISTEM._dimsextent(bng.array, nothing)
    @test by.Y[2] - by.Y[1] == 9 * 2500.0m

    # A lookup with no `Regular` span falls back to its coordinate values - the behaviour this
    # function had throughout before. Both non-regular cases must be caught, and they fail
    # *differently*: an `Irregular` span with no recorded ends (what `Rasters.aggregate` leaves
    # behind) answers `(nothing, nothing)`, which would poison the extent arithmetic, while a
    # lookup still at `AutoSpan`/`AutoSampling` makes `Lookups.bounds` throw outright.
    auto = Y(Rasters.Projected([0.0, 1.0, 2.0] .* °,
                               crs = Rasters.EPSG(4326)))
    @test EcoSISTEM._axisbounds(auto) == (0.0°, 2.0°)
    irr = Y(Rasters.Projected([0.0, 1.0, 2.0] .* °,
                              crs = Rasters.EPSG(4326),
                              sampling = DimensionalData.Lookups.Intervals(DimensionalData.Lookups.Start()),
                              order = DimensionalData.Lookups.ForwardOrdered(),
                              span = DimensionalData.Lookups.Irregular((nothing,
                                                                        nothing))))
    @test EcoSISTEM._axisbounds(irr) == (0.0°, 2.0°)
end

# A layer the plan calls `:exact` must be **cropped, not warped**. A GDAL warp onto an
# identical grid does leave the values alone (bilinear's weights collapse to (1,0,0,0) on
# coincident cells - measured max difference 2e-9), but it does *not* leave the nodata mask
# alone: it poisons any output cell whose input stencil touches a `NaN`, so the coastline erodes,
# and which cells that catches depends on GDAL's internal source windowing - which the extent of
# the *read* changes. Reading a global layer and reading a window of the same layer therefore
# produced study areas differing in 138 coastal cells. Cropping removes the warp, so the answer
# is invariant by construction, and the report's "kept exactly" becomes true.
@testset "an exactly-aligned layer is cropped, never resampled" begin
    r = _testraster(WorldClim{BioClim},
                    reshape(collect(1.0:25.0), 5, 5) .* K)
    own = EcoSISTEM._owngrid(r)
    # The layer onto its own grid: values must come back bit-identical, which a warp does not give.
    whole = EcoSISTEM._sampledata(r, own, name = "l", categorical = false)
    @test !isnothing(EcoSISTEM._regrid(r, own, mean))
    @test parent(whole) == parent(r.array)

    # ...and onto a sub-grid of itself, it is exactly the corresponding block.
    sub = own[Y(2:4), X(3:5)]
    part = EcoSISTEM._sampledata(r, sub, name = "l", categorical = false)
    @test size(part) == (3, 3)
    @test parent(part) == parent(r.array)[2:4, 3:5]

    # A NaN must not spread into its neighbours, which is the whole point - a warp would take
    # the four cells around the hole with it.
    holed = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5))
    parent(holed.array)[3, 3] = NaN * K
    out = EcoSISTEM._sampledata(holed, EcoSISTEM._owngrid(holed),
                                categorical = false,
                                name = "l")
    @test count(isnan, ustrip.(parent(out))) == 1

    # A grid that does not sit on the layer's cell boundaries is not a tiling of it - a half-cell
    # offset here - so its cells are aggregated from the covering cells rather than selected. (A
    # coarser grid that *does* sit on them is an exact block aggregation; `test_materialise.jl` pins
    # that.)
    @test isnothing(EcoSISTEM._tilerange(collect(0.0:1.0:4.0) .* °,
                                         collect(0.5:1.0:4.5) .* °, 1))

    # Both routes must leave a layer with the *same* dims, or two layers of one collection
    # sampled differently are rejected downstream as being "on different grids" - which is
    # exactly what happened when the crop kept a `Regular` span and the warp inferred `Irregular`.
    L = DimensionalData.Lookups
    for d in (dims(whole, Y), dims(whole, X))
        @test L.span(d) isa L.Regular
    end
end

# The payoff of carrying layer codes through a read: whether a layer may be *interpolated* is a
# per-layer fact in the shipped `ValueType` column, and a source type alone cannot supply it -
# BioClimPlus holds all three value types at once. Before this, `_iscategorical` was a stub
# returning `false`, so a Köppen-Geiger layer was bilinearly interpolated *between class codes*.
# Note `:discrete` is **not** categorical: those are day-of-year and day-count layers, ordinary
# numbers that average perfectly well. Only `:categorical` forbids it.
@testset "a categorical layer is resampled by class, not interpolated" begin
    d = (Y(Rasters.Projected(collect(0.0:1.0:4.0) .* °,
                             crs = Rasters.EPSG(4326))),
         X(Rasters.Projected(collect(0.0:1.0:4.0) .* °,
                             crs = Rasters.EPSG(4326))))
    arr = DimArray(fill(1.0, 5, 5), d)
    cat = ClimateRaster(CHELSA{BioClimPlus}, arr, :kg0)   # Koppen-Geiger class codes
    disc = ClimateRaster(CHELSA{BioClimPlus}, arr, :gsl)  # growing-season length, in days
    cont = ClimateRaster(WorldClim{BioClim}, arr, :bio1)
    bare = ClimateRaster(WorldClim{BioClim}, arr)         # derived: no code to look up
    # With no axis given, a **coded** raster is classified from the catalogue - which asks
    # `_iscategorical` of the layer's own axis, so both routes end in the same place.
    @test EcoSISTEM._iscategorical(cat)
    for r in (disc, cont)
        @test !EcoSISTEM._iscategorical(r)
    end
    # ...and a declared axis is authoritative, needing no code at all: this is what a derived
    # layer uses, and `LandCoverTypology` is a `TypologyAxis`, so class codes.
    @test EcoSISTEM._iscategorical(bare, LandCoverTypology)
    @test !EcoSISTEM._iscategorical(bare, Temperature)
    # The axis wins over the catalogue where both could answer - it *is* the declaration.
    @test !EcoSISTEM._iscategorical(cat, Temperature)
    # ...and the catalogue can also be asked directly, without a raster at all.
    @test EcoSISTEM._iscategorical(CHELSA{BioClimPlus}, :kg0)
    @test !EcoSISTEM._iscategorical(CHELSA{BioClimPlus}, :gsl)

    # Whether a layer holds class codes is a property of its axis, which the raster does not carry;
    # the reducer that follows from it is `_reducer`'s to choose, and nothing is interpolated.
    @test EcoSISTEM._reducer(nothing, LandCoverTypology) ===
          EcoSISTEM._majorityclass
    @test EcoSISTEM._reducer(nothing, Temperature) === EcoSISTEM._meanpresent

    # A uniformly categorical stack is categorical.
    codes(cs...) = collect(EcoSISTEM.CODE_TYPE, cs)
    @test EcoSISTEM._iscategorical(ClimateRaster(CHELSA{BioClimPlus}, arr,
                                                 codes(:kg0, :kg1)))
    # A *mixed* stack has no correct answer - every caller picks one behaviour for the whole
    # array - so it errors rather than returning either. `read` calls the same method before
    # downloading anything (below), so reaching it here means the raster was assembled some other
    # way. This is the one method of `_iscategorical` that can throw.
    mixed = ClimateRaster(CHELSA{BioClimPlus}, arr, codes(:kg0, :bio3))
    @test_throws ErrorException EcoSISTEM._iscategorical(mixed)
    # No named axis and no code answers `false`: there is nothing to look up and nothing
    # claiming the values are class labels. A caller who means class codes says so with the axis.
    @test !EcoSISTEM._iscategorical(bare)
    # `NicheAxis` has its own method, and it defers to the raster's own codes rather than
    # answering `false`: it is the *absence* of a named axis, not a claim that the values are
    # continuous. This is what keeps a legitimate multi-axis stack - `SourceSpec(WorldClim{BioClim})`
    # derives `NicheAxis` from codes that disagree - from having its class codes interpolated.
    @test EcoSISTEM._iscategorical(cat, EcoSISTEM.NicheAxis)
    @test !EcoSISTEM._iscategorical(cont, EcoSISTEM.NicheAxis)
    @test !EcoSISTEM._iscategorical(bare, EcoSISTEM.NicheAxis)
    # ...while a *named* axis stays authoritative, including one that contradicts the catalogue.
    @test !EcoSISTEM._iscategorical(cat, Temperature)
    @test EcoSISTEM._iscategorical(bare, LandCoverTypology)
end

# An array has one eltype and gets one resample method, so its layers must agree on both unit and
# categorical-ness. Refused where a single array is actually demanded - not at construction,
# which would rule out `WorldClim{BioClim}` and three other shipped datasets outright.
@testset "layers that cannot share one array are refused" begin
    msg(spec) =
        try
            (read(spec); "")
        catch e
            sprint(showerror, e)
        end
    # Mixed units: one array cannot carry °C and mm and m s^-1.
    u = msg(SourceSpec(WorldClim{Climate}))
    @test occursin("different units", u)
    @test occursin("ConstructedRasterSpec", u)       # the error names the way out

    @test layerunit(CHELSA{BioClimPlus}, :kg0) ==
          layerunit(CHELSA{BioClimPlus}, :fcf)
    c = msg(SourceSpec(CHELSA{BioClimPlus}, [:kg0, :fcf]))
    @test occursin("class codes", c)
    @test occursin("no one resampling method", c)
    @test occursin("ConstructedRasterSpec", c)       # the same way out
end

# One rule, in one place. Reading a multi-layer `SourceSpec` and asking a raster about itself are
# the same question about the same codes, so they must be the same method: two copies of the filter
# are two wordings that drift apart.
@testset "the shared-stack rule is stated once" begin
    S = CHELSA{BioClimPlus}
    @test EcoSISTEM._iscategorical(S, [:kg0, :kg1])     # uniformly categorical
    @test !EcoSISTEM._iscategorical(S, [:bio1, :bio5])  # uniformly numeric
    @test !EcoSISTEM._iscategorical(S, [])              # nothing to disagree about
    # ...and the scalar method answers the same question of one layer, without throwing.
    @test EcoSISTEM._iscategorical(S, :kg0)
    @test !EcoSISTEM._iscategorical(S, :fcf)

    # The two entry points reach the same method, so they cannot disagree - assert it by
    # comparing the messages themselves, not merely that both threw.
    d = (Y(Rasters.Projected(collect(0.0:1.0:4.0) .* °,
                             crs = Rasters.EPSG(4326))),
         X(Rasters.Projected(collect(0.0:1.0:4.0) .* °,
                             crs = Rasters.EPSG(4326))))
    arr = DimArray(fill(1.0, 5, 5), d)
    bad = collect(EcoSISTEM.CODE_TYPE, [:kg0, :fcf])
    msgs = map([() -> read(SourceSpec(S, bad)),
                   () -> EcoSISTEM._iscategorical(ClimateRaster(S, arr, bad))]) do f
        try
            f()
            ""
        catch e
            sprint(showerror, e)
        end
    end
    @test !isempty(first(msgs))
    @test msgs[1] == msgs[2]
    @test occursin("class codes must not be interpolated", msgs[1])
    @test occursin("kg0", msgs[1]) && occursin("fcf", msgs[1])
end

# `compress_landcover` argmaxes EarthEnv's per-class % bands into class codes. Averaging
# *between* class codes is meaningless - and invisible, because GDAL rounds on an integer band,
# so the output still looks like tidy whole numbers. Measured over Scotland it invented five
# classes that are not there (Evergreen Broadleaf Trees, Deciduous Broadleaf Trees, Regularly
# Flooded Vegetation, Snow/Ice, Barren) and lost one that is.
#
# **There are now TWO independent guards against that, and this testset pins both.**
#   1. **Ordering** - a `ConstructedRasterSpec` combine runs *after* its children are sampled, so the
#      *percentages* are interpolated and the argmax taken afterwards, which keeps each cell's
#      sub-cell mix and cannot fabricate a class.
#   2. **Declaration** - `compress_landcover` marks its result `valuetype = :categorical`, so if
#      codes are resampled anyway (the `CombineOnSourceGrid` escape hatch) the resampler takes
#      the nearest class (`:mode`) instead of interpolating.
# Guard 2 is what makes the source-grid ordering safe too: with the declaration in place,
# collapsing first no longer invents classes. The assertion below therefore checks that the
# declaration is **load-bearing** - strip it and the five invented classes come straight back -
# rather than asserting a defect that
# has been fixed.
@testset "a land-cover collapse runs after resampling, not before" begin
    scot = EcoSISTEM.boundingbox("Scotland", coverage = AllTerritories())
    # `axis = LandCoverTypology` because the *result* is a class code, not the cover fraction
    # its inputs are (`SurfaceArea`). A raster carries no axis of its own, so a derived layer's
    # axis is declared on the spec - the same place `AvailableGround.jl` declares `SurfaceArea`.
    # `scale = 10` reads the source at 5 arcminutes, coarser than the 5 km grid, so the collapse
    # on the target loses no class and the ground-truth comparison below is exact; it also keeps
    # the whole-globe read this testset makes outside the study area to a few hundred megabytes.
    spec = ConstructedRasterSpec(SourceSpec(EarthEnv{LandCover}, scale = 10),
                                 axis = LandCoverTypology) do lc
        return compress_landcover(lc)
    end
    area = StudyArea(regime = spec, within = scot,
                     crs = Rasters.EPSG(27700), cellsize = 5km,
                     verbosity = :silent)
    target, cache = area.report.active, area.report.cache

    classes(a) = Set(Int.(filter(!isnan, ustrip.(vec(parent(a))))))
    # Ground truth is the collapse on its **own** grid, never resampled - not the other
    # ordering. Comparing the two orderings only shows they differ; it cannot show which one is
    # right, and a version that invented *fewer* classes would pass such a test just as happily.
    src = EcoSISTEM._asraster(spec, cache)
    truth = classes(EcoSISTEM._applycut(src.array, scot))
    fresh = classes(EcoSISTEM._materialiseon(spec, target, cache).array)
    stale = classes(EcoSISTEM._sampledata(src, target, name = "old",
                                          categorical = true))

    # The property that matters: resampling may *lose* a class that only ever wins in a handful
    # of fine cells, but it must never conjure one that occurs nowhere in the data.
    @test issubset(fresh, truth)
    @test fresh == truth        # here it recovers the ground-truth set exactly
    # ...and with `valuetype = :categorical` declared, collapsing *first* is safe too - `:mode`
    # picks a nearest class and so cannot invent one. Measured over Scotland at 5 km the two
    # orderings now agree exactly.
    @test issubset(stale, truth)

    # **The axis is what does that, so prove it is not vacuous**: put the identical raster on the
    # grid as though it held continuous values and its class codes are *averaged*, which yields
    # values that are no class at all.
    @test EcoSISTEM._iscategorical(src, LandCoverTypology)
    naive = EcoSISTEM._sampledata(src, target, name = "raw",
                                  categorical = false)
    @test any(v -> !isinteger(v), filter(!isnan, ustrip.(vec(parent(naive)))))
end

# **A derived layer says what it holds through its AXIS**, which is the only thing that can
# know: the shipped table has no row for it, and its combine is an arbitrary closure. There used
# to be a separate `valuetype` keyword for this; it is gone, because it could only ever repeat
# the axis or contradict it.
@testset "a derived layer's axis says whether it holds class codes" begin
    d = (Y(Rasters.Projected(collect(0.0:1.0:4.0) .* °,
                             crs = Rasters.EPSG(4326))),
         X(Rasters.Projected(collect(0.0:1.0:4.0) .* °,
                             crs = Rasters.EPSG(4326))))
    arr = DimArray(fill(1.0, 5, 5), d)
    bare = ClimateRaster(WorldClim{BioClim}, arr)
    @test EcoSISTEM._iscategorical(bare, LandCoverTypology)
    @test !EcoSISTEM._iscategorical(bare, Temperature)
    # A layer derived *from* a categorical source can say it is not one - a per-class fraction,
    # say - simply by declaring the axis it is actually on.
    @test !EcoSISTEM._iscategorical(ClimateRaster(CHELSA{BioClimPlus}, arr,
                                                  :kg0), Precipitation)
    # ...and the spec is where that declaration lives.
    spec = ConstructedRasterSpec(() -> ClimateRaster(WorldClim{BioClim}, arr),
                                 axis = LandCoverTypology)
    @test spec.axis === LandCoverTypology
    @test EcoSISTEM._iscategorical(spec.axis)
    @test !EcoSISTEM._iscategorical(Temperature)
end

# The step this whole subproject exists for: a monthly total becomes a rate by being divided by
# **its own** month, not by a fixed 30.4375-day one. Exercised on a synthetic array so it costs no
# download - what is being tested is the arithmetic and the slice-to-month mapping, not the reader.
@testset "a monthly read is divided by each slice's own month" begin
    spec = SourceSpec(WorldClim{Climate}, :prec)
    # 100 units in every cell of every month, so the divisor is the only thing that can vary.
    array = fill(100.0, 2, 2, 12)
    out = EcoSISTEM._applyperiod(array, spec, NamedTuple())

    # February is the case the plan is named after: 100/28.25 rather than 100/30.4375, which is
    # 7.7% higher - the error that was silently in every monthly climatology until now.
    @test out[1, 1, 2] ≈ 100 / 28.25
    @test out[1, 1, 2] / (100 / 30.4375) ≈ 1.0774 atol = 1e-4
    @test out[1, 1, 1] ≈ 100 / 31           # January
    @test out[1, 1, 4] ≈ 100 / 30           # April
    # Every slice divided by its own month, and nothing else touched.
    @test all(out[1, 1, m] ≈ 100 / d
              for (m, d) in enumerate([31, 28.25, 31, 30, 31, 30, 31, 31,
                                          30,
                                          31, 30, 31]))

    # A partial read is divided by the months actually requested. Under the old fixed divisor this
    # distinction did not exist; getting it wrong here would be silent.
    part = EcoSISTEM._applyperiod(fill(100.0, 2, 2, 3), spec,
                                  (month = 2:4,))
    @test part[1, 1, 1] ≈ 100 / 28.25       # February, not January
    @test part[1, 1, 3] ≈ 100 / 30          # April

    # A slice count that disagrees with the months asked for is refused rather than divided by
    # whatever happens to line up.
    @test_throws ErrorException EcoSISTEM._applyperiod(fill(100.0, 2, 2, 5),
                                                       spec, (month = 2:4,))

    # A layer whose canonical reading is the accumulated amount is returned untouched - a heat sum
    # must never be quietly turned into a temperature.
    heat = fill(100.0, 2, 2)
    @test EcoSISTEM._applyperiod(heat,
                                 SourceSpec(CHELSA{BioClimPlus}, :gdd0),
                                 NamedTuple()) === heat
end

@testset "attach units to a raster layer" begin
    # a bare (unitless) layer gets its unit, preserving axes and source type
    bare = _testraster(WorldClim{BioClim}, fill(290.0, 5, 5))
    united = EcoSISTEM._attachunit(bare, u"K")
    @test eltype(united.array) <: Unitful.Temperature
    @test united isa ClimateRaster{WorldClim{BioClim}}
    @test DimensionalData.name.(dims(united.array)) == (:Y, :X)
    # NoUnits is a no-op -> stays bare numbers (what discrete land cover wants)
    landcover = _testraster(EarthEnv{LandCover},
                            Float64.(repeat(1:5, 1, 5)))
    @test eltype(EcoSISTEM._attachunit(landcover, NoUnits).array) ==
          Float64
end

@testset "SourceSpec carries its own read keywords (e.g. month)" begin
    # Construction stays pure - the keywords are only recorded, nothing is read.
    s = SourceSpec(WorldClim{Climate}, :wind, month = 1:2)
    @test s.readkw == (month = 1:2,)
    # `axis` is a *spec* keyword, so it must not be swallowed into the read keywords.
    t = SourceSpec(WorldClim{BioClim}, 1, axis = Temperature)
    @test t.readkw == NamedTuple()
    @test EcoSISTEM._specaxis(t) === Temperature
    # the whole-dataset (no-code) form takes them too
    @test SourceSpec(EarthEnv{LandCover}, month = 3).readkw == (month = 3,)
    # The read options common to every source are fields, not read keywords, and a catalogued spec
    # resolves its own files.
    @test SourceSpec === RasterSpec
    opts = SourceSpec(WorldClim{BioClim}, 1, scale = 4, fn = maximum)
    @test opts isa RasterSpec{Temperature}
    @test opts.scale == 4 && opts.fn === maximum && isnothing(opts.cut)
    @test opts.readkw == NamedTuple()
    @test isnothing(opts.files)
    @test isnothing(t.scale)
    @test_throws ErrorException SourceSpec(WorldClim{BioClim}, 1, scale = 0)

    # ...and the keywords really do reach `read` through `GridHabitat`. This is the
    # regression test: A1 deleted the only method that forwarded read keywords, so `month` became
    # silently unreachable. Asking for 2 months must give a 2-slice series - not the source's
    # own default of all 12.
    wind = SourceSpec(WorldClim{Climate}, :wind, month = 1:2)
    env = _env(wind, SUP,
               within = EcoSISTEM.boundingbox("Scotland",
                                              coverage = LargestLandmass()))
    # The regime is a plain 2-D layer holding the current slice; the stack lives in its change,
    # which is what knows how elapsed time picks between the slices.
    @test env.regime isa EcoSISTEM.ContinuousRegime
    @test ndims(env.regime.matrix) == 2
    @test env.regime.change isa SeriesLayerChange{AbsoluteChange}
    @test size(env.regime.change.slices, 3) == 2
    # `isequal`, not `==`: the read carries `NaN` in the cells outside the mask, and `NaN` is
    # never `==` itself.
    @test all(isequal.(env.regime.matrix,
                       env.regime.change.slices[:, :, 1]))
    # the default active mask comes from the first time slice rather than erroring on a 3-D array
    @test 0 < count(env.active) <= length(env.active)
end

@testset "a shape covers each cell by the share of its area inside" begin
    polygon(points) = ArchGDAL.createpolygon([[points..., first(points)]])
    part(g) = (prepared = ArchGDAL.preparegeom(g),
               envelope = ArchGDAL.envelope(g), geometry = g)
    # Three rows from y = 0 and four columns from x = 0, as `(lo, hi)`. Cells are 2 tall and 1.5
    # wide, so a cell's area is 3: with cells of area one, dividing by the wrong side - or by
    # nothing at all - gives the same shares and the tests pass anyway.
    rows = [(y, y + 2.0) for y in 0.0:2.0:4.0]
    cols = [(x, x + 1.5) for x in 0.0:1.5:4.5]

    # A rectangle over x 0.75 to 3.375 and y 0 to 3, so each cell's share is the product of the
    # shares of its row and its column, written out by hand.
    rectangle = part(polygon([
                                 (0.75, 0.0),
                                 (3.375, 0.0),
                                 (3.375, 3.0),
                                 (0.75, 3.0)
                             ]))
    covered = EcoSISTEM._coveredfraction((rectangle,), rows, cols)
    @test covered ≈ [0.5 1.0 0.25 0.0
                     0.25 0.5 0.125 0.0
                     0.0 0.0 0.0 0.0]
    # A cell wholly inside counts exactly one, not an intersection's area divided back out.
    @test covered[1, 2] == 1.0

    # An L covering three whole cells touches the fourth, at x 1.5 to 3 and y 2 to 4, only along two
    # of its edges: it intersects that cell, so the zero-area case is reached, and shares no area.
    corners = [
        (0.0, 0.0),
        (3.0, 0.0),
        (3.0, 2.0),
        (1.5, 2.0),
        (1.5, 4.0),
        (0.0, 4.0)
    ]
    l = part(polygon(corners))
    touching = polygon([(1.5, 2.0), (3.0, 2.0), (3.0, 4.0), (1.5, 4.0)])
    @test ArchGDAL.intersects(l.prepared, touching)
    lshares = [1.0 1.0 0.0 0.0
               1.0 0.0 0.0 0.0
               0.0 0.0 0.0 0.0]
    @test EcoSISTEM._coveredfraction((l,), rows, cols) == lshares

    # Disjoint pieces add cell by cell, and a piece off the grid adds nothing.
    far = part(polygon([(4.5, 4.0), (6.0, 4.0), (6.0, 6.0), (4.5, 6.0)]))
    away = part(polygon([
                            (20.0, 20.0),
                            (21.0, 20.0),
                            (21.0, 21.0),
                            (20.0, 21.0)
                        ]))
    both = EcoSISTEM._coveredfraction((rectangle, far, away), rows, cols)
    @test both[3, 4] == 1.0
    both[3, 4] = 0.0
    @test both ≈ covered

    # A descending axis gives the same shares in reversed rows, and units are stripped whatever
    # they are.
    @test EcoSISTEM._coveredfraction((rectangle,), reverse(rows), cols) ≈
          reverse(covered, dims = 1)
    @test EcoSISTEM._coveredfraction((rectangle,),
                                     [(lo * °, hi * °) for (lo, hi) in rows],
                                     [(lo * °, hi * °) for (lo, hi) in cols]) ≈
          covered

    # What this does not catch, and cannot: windowing each piece to its envelope only saves visits,
    # since a cell outside it shares no area either way, so widening or dropping the window leaves
    # every share as it is. Nor does it reach the cap at one, which needs pieces that overlap, where
    # a shape's own pieces are disjoint.
end

@testset "each mask rule is one comparison against the covered share" begin
    tol = EcoSISTEM._COVERTOL
    @test EcoSISTEM._rulepasses(AnyOverlap(), 1.0)
    @test EcoSISTEM._rulepasses(AnyOverlap(), 10 * tol)
    # A cell touched only along an edge shares no area; nor may the reprojection's own error, which
    # can turn a shared edge into a sliver of real overlap, read as ground.
    @test !EcoSISTEM._rulepasses(AnyOverlap(), 0.0)
    @test !EcoSISTEM._rulepasses(AnyOverlap(), tol / 2)

    @test EcoSISTEM._rulepasses(FractionWithin(0.5), 0.5)
    @test EcoSISTEM._rulepasses(FractionWithin(0.5), 0.5 - tol / 2)
    @test !EcoSISTEM._rulepasses(FractionWithin(0.5), 0.5 - 10 * tol)

    # `FullyWithin` is `FractionWithin(1.0)` rather than a second test that could drift from it.
    for f in (0.0, 0.25, 1 - 10 * tol, 1 - tol / 2, 1.0)
        @test EcoSISTEM._rulepasses(FullyWithin(), f) ==
              EcoSISTEM._rulepasses(FractionWithin(1.0), f)
    end
    # A cell that is whole but for the reprojection's 5e-7 still counts as whole.
    @test EcoSISTEM._rulepasses(FullyWithin(), 1 - 5e-7)
    @test !EcoSISTEM._rulepasses(FullyWithin(), 0.999)

    # A threshold outside (0, 1] is refused where it was written, not where it is read.
    @test_throws ArgumentError FractionWithin(0.0)
    @test_throws ArgumentError FractionWithin(1.5)

    # What this does not catch: whether the shares themselves are right, which the covered-fraction
    # testset above pins, and whether `_rastermask` hands the rule the share rather than something
    # else, which the mask-pipeline tests in `test_NaturalEarth.jl` cover on real geometry.
end

@testset "a shape mask applies its own rule to each cell's covered share" begin
    # Four rows and five columns of 1 degree cells labelled by their lower corner, deliberately not
    # square. A rectangle covers x 0.5 to 2.1 and y 0 to 2, so the bottom two rows are covered by
    # 0.5, 1.0 and 0.1 of their width and nothing else is covered at all - shares that tell the
    # three area rules apart, where a fixture of whole cells would let any of them pass as another.
    grid = _testraster(WorldClim{BioClim}, zeros(4, 5),
                       lat = (0.0:1.0:3.0) .* °, long = (0.0:1.0:4.0) .* °).array
    rectangle = ArchGDAL.createpolygon([[(0.5, 0.0), (2.1, 0.0), (2.1, 2.0),
                                           (0.5, 2.0), (0.5, 0.0)]])
    piece = (prepared = ArchGDAL.preparegeom(rectangle),
             envelope = ArchGDAL.envelope(rectangle), geometry = rectangle)
    masked(rule) = EcoSISTEM._rastermask((parts = (piece,), rule = rule),
                                         nothing, grid)

    @test masked(AnyOverlap()) == [true true true false false
           true true true false false
           false false false false false
           false false false false false]
    @test masked(FractionWithin(0.5)) == [true true false false false
           true true false false false
           false false false false false
           false false false false false]
    @test masked(FullyWithin()) == [false true false false false
           false true false false false
           false false false false false
           false false false false false]
    @test masked(FullyWithin()) == masked(FractionWithin(1.0))

    # The three disagree on this fixture, so a rule discarded on the way through - which is what a
    # payload carrying the geometry alone would do - cannot pass as one of the others.
    @test count(masked(AnyOverlap())) == 6
    @test count(masked(FractionWithin(0.5))) == 4
    @test count(masked(FullyWithin())) == 2
end

@testset "a mask payload is cropped with the grid only when it is tied to its shape" begin
    # A recut re-rasterises the mask on the cropped grid, so a payload indexed by cell position has
    # to be cropped alongside it or it meets `_rastermask`'s size check. Geometry is placed by
    # coordinates, so it must travel unchanged.
    bits = trues(5, 7)
    bits[1, 1] = false
    @test size(EcoSISTEM._croppayload(bits, 2:4, 3:6)) == (3, 4)
    @test size(EcoSISTEM._croppayload(Matrix(bits), 2:4, 3:6)) == (3, 4)
    @test EcoSISTEM._croppayload(Matrix(bits), 2:4, 3:6) == bits[2:4, 3:6]
    geometry = (parts = (), rule = AnyOverlap())
    @test EcoSISTEM._croppayload(geometry, 2:4, 3:6) === geometry
    @test isnothing(EcoSISTEM._croppayload(nothing, 2:4, 3:6))
end

@testset "_axiswindow reads the axis direction rather than assuming it" begin
    up = collect(0.0:1.0:10.0)
    @test EcoSISTEM._axiswindow(up, 2.5, 5.5) == 4:6          # 3.0, 4.0, 5.0
    @test EcoSISTEM._axiswindow(reverse(up), 2.5, 5.5) == 6:8  # the same three values
    @test isempty(EcoSISTEM._axiswindow(up, 20.0, 30.0))
    @test isempty(EcoSISTEM._axiswindow(reverse(up), 20.0, 30.0))
    # Monotonic in neither direction, so every index is offered. A grid axis always is monotonic,
    # but the fallback must not drop cells silently if one is not.
    @test EcoSISTEM._axiswindow([0.0, 5.0, 1.0], 0.5, 2.0) == 1:3
end

end
