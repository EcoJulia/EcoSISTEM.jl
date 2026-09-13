# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDatasetread

using EcoSISTEM
using EcoSISTEM: materialise, in_memory_raster
using EcoSISTEM.Units
import Extents
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources
using Rasters
using ArchGDAL
# Loads Rasters' affine-transform extension. A raster whose geotransform is rotated or skewed cannot
# be opened without it, and Rasters says so rather than failing obscurely - but only at read time, so
# the import has to be here rather than discovered on a runner.
using CoordinateTransformations
using DimensionalData
using Statistics
using Dates: Dates
import NCDatasets
using Test

# CHELSA's `bio1` is a 43200 x 20880 global grid, and coarsening does NOT bound the cost of reading
# it: the aggregate consumes the whole file either way. Measured on a development machine -- 11.7 GB
# peak resident for the bare `Rasters.aggregate` after reading whole, 10.8 GB for the lazy
# window-by-window fallback that exists to be the cheap path, and 25 GB through this package's full
# read pipeline. None of that fits the 16 GB a GitHub runner has, even with the machine to itself:
# the cache-priming job, which runs nothing else at all, was killed attempting exactly this.
#
# So the whole-file CHELSA reads run only where there is memory for them. Same predicate and the same
# reasoning as `heavydata()` in `test/canonical/canonical.jl`, kept local because that lives in a
# module of its own for the canonical suite; `ECOSISTEM_HEAVY_DATA=true` forces them on anywhere.
#
# What CI gives up is stated rather than hidden: the extent assertion below still runs for WorldClim
# and EarthEnv, but CHELSA is the case that MOTIVATES it -- its origins sit ~1.7% of a cell off the
# lattice -- so on a runner that particular lattice is unchecked. It is checked locally, and by
# `test/canonical/`.
function bigrasters()
    return haskey(ENV, "ECOSISTEM_HEAVY_DATA") ?
           ENV["ECOSISTEM_HEAVY_DATA"] == "true" : !haskey(ENV, "RUNNER_OS")
end

if !Sys.iswindows()
    # `getraster` returns the full path(s) to the downloaded file(s), so use those directly rather
    # than reconstructing RasterDataSources' folder layout. Pre-fetching here (outside the
    # `@test_nowarn`s) also keeps download messages out of those tests on an empty cache.
    bio1 = getraster(WorldClim{BioClim}, :bio1)               # one tif path
    wind = getraster(WorldClim{Climate}, :wind, month = 1:12) # 12 monthly tif paths
    getraster(EarthEnv{LandCover})
    getraster(CHELSA{BioClim}, 1)
    # A directory holding exactly the 12 downloaded wind tifs, to exercise the directory readers
    # `read(CRUTS, ...)`/`read(CHELSA{Climate}, ...)` (the variable name only fixes the unit that
    # gets attached).
    winddir = dirname(first(wind))

    @testset "Reading functions" begin
        @test_nowarn read(SourceSpec(WorldClim{Climate}, :wind, month = 1:12))
        @test_nowarn read(SourceSpec(CRUTS, "tavg", directory = winddir))
        # Every whole-dataset read here is downsampled, and the reason is a hard ceiling rather
        # than tidiness: a GitHub runner has 16 GB, and one of these files reading at native
        # resolution took the process to 7.0 GB and the runner to a shutdown signal. Measured peak
        # RSS for a fresh process, baseline 1.0 GB: CHELSA bioclim is a 43200×20880 global grid and
        # allocates several ~7 GiB Float64 arrays whole; one `EarthEnv{LandCover}` band read at
        # its own resolution costs 26 GB, and aggregating cold costs 8.5 GB for one band and
        # 11.6 GB for all twelve whatever the scale, since the scale only decides what comes out.
        # What makes these reads affordable on a runner is the aggregate cache `primecache.jl`
        # fills: a primed read costs 0.3 GB. `WorldClim{BioClim}` costs 1.5 GB against 0.3 GB at 4.
        #
        # `scale` is safe for what these assert -- that the read emits no warning, and that what
        # comes back is unitless -- since neither is a property of the resolution. A test that does
        # depend on the grid asks for one layer rather than a whole dataset, so it never reaches
        # this size.
        bigrasters() &&
            @test_nowarn read(SourceSpec(CHELSA{BioClim}, 1, scale = 20))
        @test_nowarn read(SourceSpec(EarthEnv{LandCover}, scale = 40))
        @test_nowarn read(RasterFileSpec(bio1, axis = EcoSISTEM.NicheAxis))
    end

    @testset "Output data" begin
        # A read carries the layer's unit from the table; a bare file carries what it was told.
        bio1r = read(SourceSpec(WorldClim{BioClim}, 1, scale = 4))
        cr = read(SourceSpec(CRUTS, "tavg", directory = winddir))
        rf = read(RasterFileSpec(bio1, axis = EcoSISTEM.NicheAxis))

        @test unit(bio1r.array[1]) == °C
        @test unit(rf.array[1]) == NoUnits
        @test rf isa EcoSISTEM.ClimateRaster{EcoSISTEM.SyntheticData}
        # A whole dataset of differing units is not one array.
        @test_throws "different units" read(SourceSpec(WorldClim{BioClim},
                                                       scale = 4))
        if bigrasters()
            ch_b = read(SourceSpec(CHELSA{BioClim}, 1, scale = 20))
            @test unit(ch_b.array[1]) == °C
        end
    end

    @testset "Output data 2" begin
        landcover = read(SourceSpec(EarthEnv{LandCover}, scale = 40))
        @test unit(landcover.array[1]) == percent      # the table's unit for every class
    end

    # A partial monthly read is labelled with the months it actually holds. Reads only files the
    # `getraster` above already downloaded, so nothing new is fetched.
    #
    # The axis must carry the months actually read, not be rebuilt from the slice *count* - which
    # would label `month = 2:4` as 1-3. That is not a cosmetic mislabel: the plot recipe looks a
    # month up by coordinate
    # (`At(ind * month_mean_duration)`), so asking a partial read for February returned **March's**
    # grid under February's name. The last assertion here is that regression, and it is the one that
    # would have caught it - the others only check the labels.
    # Which axis a multi-file source stacks on decides whether it is a **time series** or a stack
    # of unrelated bands - and so whether it can drive a layer through time at all. It is read off
    # the shipped catalogue: a source whose layer table declares a temporal resolution is a series.
    # A source silently read as bands loads fine - twelve monthly files, simply as twelve bands - so
    # this is asserted against its siblings rather than alone: a test of one source cannot say
    # whether `Ti` is right or the fallback is wrong.
    @testset "monthly sources stack on time, band sources on layers" begin
        @test EcoSISTEM._stackaxis(WorldClim{Climate}) == Ti
        @test EcoSISTEM._stackaxis(CHELSA{Climate}) == Ti
        @test EcoSISTEM._stackaxis(EcoSISTEM.ERA) == Ti
        # Only the `Climate` layers are monthly. A source's *bioclim* variables are one file per
        # variable, so they must keep stacking on `Dim{:layer}`.
        @test EcoSISTEM._stackaxis(CHELSA{BioClim}) == Dim{:layer}
        @test EcoSISTEM._stackaxis(WorldClim{BioClim}) ==
              Dim{:layer}
        @test EcoSISTEM._stackaxis(EarthEnv{LandCover}) ==
              Dim{:layer}
        # A source with no table at all is a stack of bands.
        @test EcoSISTEM._stackaxis(Int) == Dim{:layer}
    end

    @testset "a partial monthly read knows which months it holds" begin
        full = read(SourceSpec(WorldClim{Climate}, :wind, month = 1:12))
        part = read(SourceSpec(WorldClim{Climate}, :wind, month = 2:4))

        @test collect(DimensionalData.lookup(full.array, Ti)) ==
              (1:12) .* month_mean_duration
        @test collect(DimensionalData.lookup(part.array, Ti)) ==
              (2:4) .* month_mean_duration

        # An uneven request stays uneven, which is what lets a series hold each slice until the next
        # rather than pretending the months are consecutive.
        sparse = read(SourceSpec(WorldClim{Climate}, :wind, month = [1, 6, 12]))
        @test collect(DimensionalData.lookup(sparse.array, Ti)) ==
              [1, 6, 12] .* month_mean_duration

        # A single month has no series in it, so it stays 2-D and carries no time axis at all -
        # a deliberate carve-out, since a length-1 `Ti` would change the result's dimensionality and
        # every `ndims == 2` static-vs-series branch downstream with it.
        one = read(SourceSpec(WorldClim{Climate}, :wind, month = 2))
        @test ndims(one.array) == 2
        @test isnothing(DimensionalData.dims(one.array, Ti))

        # The regression: the same coordinate must select the same month's data whether the read
        # was partial or complete.
        #
        # `isequal`, not `==`: these grids are mostly ocean, and `NaN == NaN` is `false`, so `==`
        # reports two identical rasters as different. `isequal` is the comparison that means "the
        # same values", which is what is being asserted.
        @test isequal(full.array[Ti(At(2month_mean_duration))],
                      part.array[Ti(At(2month_mean_duration))])
        @test isequal(full.array[Ti(At(4month_mean_duration))],
                      part.array[Ti(At(4month_mean_duration))])
    end

    # A coarsened read must keep bounded spatial dims. An `Irregular((nothing, nothing))` span has
    # no bounds, so anything needing them fails: `_applycut`'s `Touches` selector compares
    # `nothing < 60.86°` and throws a bare `MethodError`, which breaks
    # `read(EarthEnv{LandCover}, ..., cut = ...)` outright and forces a whole-globe read and a crop.
    # A coarsened EarthEnv read is what surfaced it, which made it look source-specific rather
    # than a general consequence of aggregating. Guarded because it needs the
    # real file: the vector-lookup condition that triggers it cannot be reproduced synthetically.
    @testset "a coarsened read can be cut (Regular span with real bounds)" begin
        L = DimensionalData.Lookups
        scotland = EcoSISTEM.boundingbox("Scotland",
                                         coverage = AllTerritories())
        whole = read(SourceSpec(EarthEnv{LandCover}, 7, scale = 10))
        for d in (Y, X)
            @test L.span(dims(whole.array, d)) isa L.Regular
            @test all(!isnothing, L.bounds(dims(whole.array, d)))
        end
        # The end the fix exists for - this threw a MethodError before it.
        cut = read(SourceSpec(EarthEnv{LandCover}, 7, scale = 10),
                   cut = scotland)
        @test size(cut.array, 1) < size(whole.array, 1)
        @test size(cut.array, 2) < size(whole.array, 2)
        # ...and it is a *window*, not a token crop: Scotland is a tiny share of a global layer.
        @test prod(size(cut.array)) < 0.01 * prod(size(whole.array))
    end

    # A read grid must land exactly where its source says it does. Two things get it there: the axis is
    # rebuilt as a range from (origin, step, length) rather than carried as the file's coordinate vector
    # - which accumulates rounding over 43200 entries - and angular origins/steps are snapped onto the
    # lattice the source intends. CHELSA is the case that needs the second: its `bio1` records a step of
    # 29.99999988 arcsec for 30, and origins ~0.5 arcsec (1.7% of a cell) off, so untouched it lands
    # *near* its documented -180...180 by -90...84 rather than on it. Guarded because it needs the real files.
    # CHELSA is read coarsened here for the same reason as in "Reading functions" above: at full
    # resolution one layer is a 6.7 GiB array whose read pipeline peaks near 39 GB, which starves the
    # other test workers. Coarsening costs the assertion nothing, because block aggregation cannot
    # clean a lattice - it carries the origin through unchanged and multiplies the step, so an
    # unsnapped CHELSA still lands on 179.99985967° by -90.00847° and every exact comparison below
    # still fails. Measured, both halves.
    @testset "a read grid lands exactly on its source's stated extent" begin
        sources = Any[(read(SourceSpec(WorldClim{BioClim}, :bio1)),
                       (-180°, 180°), (-90°, 90°)),
                      (read(SourceSpec(EarthEnv{LandCover}, 7, scale = 10)),
                       (-180°, 180°), (-56°, 90°))]
        bigrasters() && push!(sources,
              (read(SourceSpec(CHELSA{BioClim}, 1, scale = 20)),
               (-180°, 180°), (-90°, 84°)))
        for (a, xs, ys) in sources
            for (D, (lo, hi)) in ((X, xs), (Y, ys))
                v = parent(DimensionalData.lookup(a.array, D))
                # A range, not a vector: exact from three numbers, with no drift to accumulate.
                @test v isa AbstractRange
                @test first(v) == lo                  # exact equality is the point
                @test last(v) + step(v) == hi
            end
        end
    end

    # A `cut` is now pushed down into the *lazy* read, so only the window comes off disk instead of
    # the whole layer being read and then cropped. It must be a pure optimisation: the result has to
    # equal what cropping after the read gave, cell for cell.
    #
    # The aggregated case is the dangerous one. `Rasters.aggregate` blocks from index 1, so a crop
    # that does not start on a block boundary moves every coarse cell - read-extent variance
    # reintroduced exactly where it was just removed. EarthEnv read at scale 10 exercises it.
    @testset "a windowed read equals the whole read cropped" begin
        CP = EcoSISTEM
        scot = CP.boundingbox("Scotland", coverage = AllTerritories())
        for (src, code, scale) in ((WorldClim{BioClim}, :bio1, 1),
            (EarthEnv{LandCover}, 7, 10))      # block alignment
            whole = read(SourceSpec(src, code, scale = scale))
            windowed = read(SourceSpec(src, code, scale = scale,
                                       cut = scot))
            cropped = EcoSISTEM._applycut(whole.array, scot)
            @test size(windowed.array) == size(cropped)
            # Coordinates agree to within float noise; aggregating a cropped raster differs from
            # aggregating the whole one in the last bit, which is why `_inrange` works on integers.
            for D in (Y, X)
                w = parent(DimensionalData.lookup(windowed.array, D))
                c = parent(DimensionalData.lookup(cropped, D))
                @test w ≈ c
            end
            @test isequal(parent(windowed.array), parent(cropped))
        end
    end

    # `sourcecrs` answers from the file header alone - a lazy open, no pixels - which is what lets the
    # target CRS be settled *before* deciding how much of each layer to read.
    @testset "a source's CRS can be had without reading it" begin
        @test EcoSISTEM.sourcecrs(WorldClim{BioClim}, :bio1) ==
              EcoSISTEM._rastercrs(read(SourceSpec(WorldClim{BioClim},
                                                   :bio1)))
        # Read keywords a `SourceSpec` may carry are accepted and ignored - `scale` cannot alter a CRS.
        @test EcoSISTEM.sourcecrs(WorldClim{BioClim}, :bio1,
                                  scale = 20, cut = nothing) ==
              EcoSISTEM.sourcecrs(WorldClim{BioClim}, :bio1)
    end

    @testset "Output data 3" begin
        cr = read(SourceSpec(CRUTS, "tavg", directory = winddir))
        worldclim = read(SourceSpec(WorldClim{Climate}, :wind))

        # Every read carries the unit the layer table declares, however the files were named.
        @test unit(cr.array[1]) == °C
        @test unit(worldclim.array[1]) == m / s
        @test layerunit(WorldClim{Climate}, :wind) == m / s
    end
end

# `_snaparcsec` (ReadData) produces angular steps and `_arcsecs` (StudyArea) recognises them, in two
# different files - so the pairing has to be pinned, not assumed. A whole number of arcseconds has no
# exact `Float64` degree representation, and *which* neighbouring value you land on depends on how it
# was built: `uconvert` scales by a precomputed factor and rounds twice, a hand-written `n / 3600`
# rounds once, and the two differ in the last bit for 29 of the first 3600 counts. Both routes occur
# (`_snaparcsec` makes the first, a caller writing `cellsize = (30 / 3600)°` the second), so both must
# survive the trip. Synthetic, so it runs on every platform.
# Asserted in aggregate rather than one `@test` per count: this is a single invariant swept over a
# large input space, and emitting 10800 assertions for it would swamp the suite's totals and hide any
# real change in them.
@testset "every whole arcsecond survives snap -> recognise, however it was built" begin
    snap, arcsecs = EcoSISTEM._snaparcsec, EcoSISTEM._arcsecs
    counts = 1:3600
    # The two constructions that occur in practice: `uconvert` (what `_snaparcsec` itself emits) and a
    # hand-written `n / 3600` (what a caller writing `cellsize = (30 / 3600)°` gets). They disagree in
    # the last bit for 29 of these counts, so both must survive.
    viaunits = [uconvert(°, float(n) * arcsecond) for n in counts]
    manual = [(n / 3600) * ° for n in counts]
    # ...and a CHELSA-style badly-written step (~4e-7 relative error, its real worst case) must be pulled
    # back onto the lattice rather than rejected.
    sloppy = [((n - 4.0e-7 * n) / 3600) * ° for n in counts]
    for built in (viaunits, manual, sloppy)
        @test [arcsecs(snap(s)) for s in built] == collect(counts)
    end
    # Recognised directly too, not only after a snap - a user-supplied `cellsize` never passes through it.
    @test [arcsecs(s) for s in viaunits] == collect(counts)
    @test [arcsecs(s) for s in manual] == collect(counts)
    # And nothing *off* the lattice is accepted: a few ULP of slack must not become a tenth of an
    # arcsecond of slack. (0.1 arcsec is still 360× further out than the worst genuine value.)
    @test all(isnothing(arcsecs(((n + d) / 3600) * °))
              for n in counts, d in (0.4, -0.4, 0.1, -0.1))
    # Sub-arcsecond grids must never round away to nothing, and projected steps are not touched.
    @test snap((0.4 / 3600) * °) == (0.4 / 3600) * °
    @test snap(2500.0u"m") == 2500.0u"m"
    @test isnothing(arcsecs((0.4 / 3600) * °))
end

# `_blockrange` is what keeps a windowed read of a *coarsened* source honest. `Rasters.aggregate`
# groups cells in blocks of `scale` starting from index 1, so a crop beginning anywhere else shifts
# every coarse cell and the answer starts depending on how much was read. Widening the crop's start
# back to a block boundary makes the crop's blocks exactly the whole file's. Synthetic, so it runs on
# every platform.
@testset "a windowed crop is widened to whole aggregation blocks" begin
    CP = EcoSISTEM
    L = DimensionalData.Lookups
    # A **unitless** lookup, because that is what `_blockrange` actually sees: it runs on the raw
    # lazy raster, before `_rastertodimarray` attaches units, which is why `_lazycrop` `ustrip`s the
    # cut box first. A united fixture is not merely unrealistic here but silently useless -
    # `selectindices` answers a unit mismatch with an *empty range* rather than an error.
    d = X(Rasters.Projected(collect(0.0:1.0:99.0),
                            sampling = L.Intervals(L.Start()),
                            order = L.ForwardOrdered(), span = L.Regular(1.0),
                            crs = Rasters.EPSG(4326)))
    # scale 1: no widening. `Touches` is inclusive, so the cell *ending* at 12.0 counts too.
    @test EcoSISTEM._blockrange(d, 12.0, 15.0, 1) == 12:16
    # scale 10: the start is pulled back to a block boundary and the end pushed out to one
    r = EcoSISTEM._blockrange(d, 12.0, 15.0, 10)
    @test first(r) == 11 && last(r) == 20
    @test (first(r) - 1) % 10 == 0        # ...the invariant that keeps aggregation blocks aligned
    @test length(r) % 10 == 0
    # ...and it always contains what scale 1 would have selected, never less
    @test issubset(EcoSISTEM._blockrange(d, 12.0, 15.0, 1), r)
    # a box reaching the far end clamps to the axis rather than running past it
    @test last(EcoSISTEM._blockrange(d, 95.0, 99.0, 10)) == 100
    # a box off the axis entirely selects nothing, and `_lazycrop` then reads everything
    @test isnothing(EcoSISTEM._blockrange(d, 500.0, 600.0, 1))
end

# An origin is snapped onto the lattice implied by its **own cell size**, not onto whole arcseconds.
# The distinction is the whole correction rather than a detail: CHELSA `bio1`'s Y origin is 302399.4975
# arcsec, which at arcsecond granularity rounds *away* from the intended value to 83.999722°, but onto
# its own 30 arcsec lattice gives exactly 84° - CHELSA's stated northern limit. Testing against the
# arcsecond lattice is what makes these origins look irretrievably ambiguous. Synthetic, so it runs
# on every platform.
# The reducer a coarsening read applies is decided from the axis unless given, and the majority
# reducer must be reproducible: ties go to the smallest code, and missing or NaN cells do not vote.
@testset "aggregation reducer follows the axis" begin
    @test EcoSISTEM._reducer(nothing, Temperature) === EcoSISTEM._meanpresent
    @test EcoSISTEM._reducer(nothing, EcoSISTEM.NicheAxis) ===
          EcoSISTEM._meanpresent
    @test EcoSISTEM._reducer(nothing, LandCoverTypology) ===
          EcoSISTEM._majorityclass
    @test EcoSISTEM._reducer(nothing, ClimateTypology) ===
          EcoSISTEM._majorityclass
    @test EcoSISTEM._reducer(maximum, LandCoverTypology) === maximum
    maj = EcoSISTEM._majorityclass
    @test maj([1, 1, 2]) == 1
    @test maj([2, 3, 2, 3]) == 2                       # a tie goes to the smallest code
    @test maj([7.0, NaN, NaN, 7.0, 9.0]) == 7.0        # NaN does not vote
    @test maj([missing, 4, missing]) == 4              # nor does missing, however many
    @test ismissing(maj([missing, missing]))           # nothing present: absent
    @test ismissing(maj(Union{Missing, Float64}[NaN]))
    # The mean is over the cells present, and absent only where none is.
    mp = EcoSISTEM._meanpresent
    @test mp([1.0, 2.0, 3.0, 4.0]) == 2.5
    @test mp([1.0, NaN, 3.0, 5.0]) == 3.0
    @test mp([1.0, NaN, NaN, 5.0]) == 3.0
    @test isnan(mp([NaN, NaN]))
    @test mp(Union{Missing, Float64}[1.0, missing, missing]) == 1.0
    @test mp([1.0K, 3.0K, NaN * K]) == 2.0K            # units survive, NaN is absent
    # No source pins a reducer; the axis decides.
    @test isnothing(EcoSISTEM._defaultfn(WorldClim{BioClim}))
    # The axis a dataset read chooses by comes from the catalogue, `NicheAxis` where it cannot.
    @test EcoSISTEM._readaxis(WorldClim{BioClim}, :bio1) === Temperature
    @test EcoSISTEM._readaxis(WorldClim{BioClim}, [:bio1, :bio12]) ===
          EcoSISTEM.NicheAxis
    @test EcoSISTEM._readaxis(WorldClim{BioClim}, :nosuchlayer) ===
          EcoSISTEM.NicheAxis
end

@testset "origins snap to the cell lattice, not the arcsecond lattice" begin
    CP = EcoSISTEM
    arcsec(n) = (n / 3600)°

    # The real case, and the one nearest-arcsecond gets wrong.
    @test EcoSISTEM._snaporigin(arcsec(302399.4975), arcsec(30)) == 84°
    @test EcoSISTEM._snaporigin(arcsec(-648000.4999986), arcsec(30)) == -180°
    # An origin already on the lattice is returned untouched.
    @test EcoSISTEM._snaporigin(-180°, arcsec(600)) == -180°
    # A deliberate half-cell offset (centre- vs edge-registration) must survive: it is 50% of a cell,
    # far outside the tolerance, so it is a registration choice rather than a rounding error.
    off = arcsec(15)
    @test EcoSISTEM._snaporigin(off, arcsec(30)) == off
    # Projected axes have no such anchor - and whole metres are already exact - so they pass through.
    @test EcoSISTEM._snaporigin(245000.0u"m", 2500.0u"m") == 245000.0u"m"
end

# `_mask_int_fills` removes the raw integer-band fill sentinels (GDAL `typemax`/`typemin`) that a file's
# declared nodata misses - e.g. CHELSA's `0xffffffff`, which the default scaled read otherwise leaves as a
# spurious ~4.29e8. Synthetic (no download), so it runs on every platform.
@testset "integer-band fill-sentinel masking" begin
    CP = EcoSISTEM
    dims = (X(1:3), Y(1:3))
    A = Float64[1 2 3; 4 5 6; 7 8 9]
    r = Rasters.Raster(A, dims)

    # unsigned band: only typemax is a fill; ordinary cells are kept
    rawU = Rasters.Raster(UInt16[typemax(UInt16) 2 3; 4 5 6;
                                 7 8 typemax(UInt16)], dims)
    maskedU = EcoSISTEM._mask_int_fills(r, rawU)
    @test ismissing(maskedU[1, 1]) && ismissing(maskedU[3, 3])
    @test maskedU[2, 2] == 5.0
    @test count(ismissing, Array(maskedU)) == 2

    # signed band: both typemin and typemax are fills
    rawS = Rasters.Raster(Int16[typemin(Int16) 2 3; 4 5 6; 7 8 typemax(Int16)],
                          dims)
    maskedS = EcoSISTEM._mask_int_fills(r, rawS)
    @test ismissing(maskedS[1, 1]) && ismissing(maskedS[3, 3])

    # a float band carries no such sentinel -> returned unchanged
    @test EcoSISTEM._mask_int_fills(r, Rasters.Raster(Float32.(A), dims)) === r

    # the scale>1 path: aggregate(mean) must propagate the introduced missings, not error
    B = Rasters.Raster(Float64.(reshape(1:16, 4, 4)), (X(1:4), Y(1:4)))
    rawB = Rasters.Raster(reshape(UInt8.(vcat(typemax(UInt8), 2:16)), 4, 4),
                          (X(1:4), Y(1:4)))
    agg = Rasters.aggregate(mean, EcoSISTEM._mask_int_fills(B, rawB), 2)
    @test ismissing(agg[1, 1])                       # the fill-touching block is masked
    @test !ismissing(agg[2, 2])                      # a clean block survives
end

# A GeoTIFF written with no CRS comes back from Rasters as `WellKnownText("")`, not `nothing`, so the
# a `_crsunit(::Nothing)` fallback alone misses it and `ArchGDAL.importCRS("")` fails with the opaque
# "Failed to initialize SRS based on WKT string (Corrupt data.)" - which is what an empty-CRS file
# such as `data/Africa.tif` produces. Synthetic (no download), so it runs on every platform.
# A CF netCDF file shaped like an ERA5 download from the Copernicus Data Store: two variables on a
# 6 x 8 grid, twelve monthly `DateTime`s, a `units` attribute each, 0-360 longitudes, and no
# filename extension - which is what forces the backend to come from the catalogue rather than
# from the name. Values encode their own position so a roll of the longitude axis can be checked.
function _erafixture(dir)
    path = joinpath(dir, "era5_fixture")
    NCDatasets.NCDataset(path, "c") do ds
        NCDatasets.defDim(ds, "longitude", 8)
        NCDatasets.defDim(ds, "latitude", 6)
        NCDatasets.defDim(ds, "valid_time", 12)
        lon = NCDatasets.defVar(ds, "longitude", Float64, ("longitude",),
                                attrib = ["units" => "degrees_east"])
        lon[:] = 0.0:45.0:315.0
        lat = NCDatasets.defVar(ds, "latitude", Float64, ("latitude",),
                                attrib = ["units" => "degrees_north"])
        lat[:] = 75.0:-30.0:-75.0
        t = NCDatasets.defVar(ds, "valid_time", Int64, ("valid_time",),
                              attrib = ["units" => "seconds since 1970-01-01",
                                  "calendar" => "proleptic_gregorian",
                                  "standard_name" => "time"])
        t[:] = Dates.DateTime.(1990, 1:12, 1)
        v = NCDatasets.defVar(ds, "t2m", Float32,
                              ("longitude", "latitude", "valid_time"),
                              attrib = ["units" => "K",
                                  "long_name" => "2 metre temperature"])
        v[:, :, :] = [Float32(270 + i + 10j + 100k)
                      for i in 1:8, j in 1:6, k in 1:12]
        p = NCDatasets.defVar(ds, "tp", Float32,
                              ("longitude", "latitude", "valid_time"),
                              attrib = ["units" => "m"])
        p[:, :, :] = fill(0.001f0, 8, 6, 12)
        w = NCDatasets.defVar(ds, "swvl1", Float32,
                              ("longitude", "latitude", "valid_time"),
                              attrib = ["units" => "m**3 m**-3"])
        w[:, :, :] = fill(0.3f0, 8, 6, 12)
        return nothing
    end
    return path
end

@testset "an ERA netCDF file reads through its catalogue row" begin
    E = EcoSISTEM
    path = _erafixture(mktempdir())
    spec = SourceSpec(E.ERA, "t2m", file = path)
    @test spec.unit == K && E._specaxis(spec) === Temperature
    cr = read(spec)
    @test cr isa ClimateRaster{E.ERA}
    # The unit is the table's, the time axis the file's own dates.
    @test Unitful.unit(eltype(cr.array)) == K
    @test collect(lookup(cr.array, Ti)) == Dates.DateTime.(1990, 1:12, 1)
    # The row leaves the longitude convention to the file, whose 0-360 longitudes are rolled onto
    # (-180, 180] with the data columns following: the cell written at 180 degrees east is the one
    # now labelled -180, and the axis ends ascending.
    xs = collect(lookup(cr.array, X))
    @test xs == (-180.0:45.0:135.0) .* °
    @test collect(lookup(cr.array, Y)) == (-75.0:30.0:75.0) .* °
    @test ustrip(cr.array[Y(At(-75.0°)), X(At(-180.0°)), Ti(1)]) ==
          270 + 5 + 10 * 6 + 100
    @test ustrip(cr.array[Y(At(75.0°)), X(At(0.0°)), Ti(12)]) ==
          270 + 1 + 10 + 1200
    # A window on the rolled axis is applied after the roll, and a coarsened read aggregates the
    # slices' cells.
    cut = Extents.Extent(Y = (-80.0°, 0.0°), X = (-100.0°, 0.0°))
    windowed = read(spec, cut = cut)
    @test size(windowed.array)[1:2] == (4, 3)      # every cell the window touches
    @test size(read(spec, scale = 2).array) == (3, 4, 12)
    # `tp` accumulates per day, so the read is the rate the table declares, not the file's depth.
    @test Unitful.unit(eltype(read(SourceSpec(E.ERA, "tp", file = path)).array)) ==
          u"m/d"
    # A file stating its unit is converted to the table's; the fixture agrees with it, so the
    # conversion is exercised on the array directly.
    ras = Rasters.Raster(path, name = :t2m, source = Rasters.NCDsource())
    celsius = E._rastertodimarray(E._centresampled(ras), expressedin = °C)
    @test maximum(celsius) ≈ ustrip(°C, maximum(ras) * K)
    # `times` replaces the file's dates, one per slice, and refuses any other count.
    timed = read(SourceSpec(E.ERA, "t2m", file = path,
                            times = collect((1:12) .* month_mean_duration)))
    @test collect(lookup(timed.array, Ti)) == (1:12) .* month_mean_duration
    @test_throws "one per slice" read(SourceSpec(E.ERA, "t2m", file = path,
                                                 times = [1month_mean_duration]))
    # A directory is its files in name order, joined along time.
    two = mktempdir()
    cp(path, joinpath(two, "era5_a"))
    cp(path, joinpath(two, "era5_b"))
    write(joinpath(two, "era5_a.provenance.toml"), "written = true\n")
    dirspec = SourceSpec(E.ERA, "t2m", directory = two)
    # A netCDF file may carry no extension, as a CDS download does; the sidecar is not one.
    @test length(dirspec.files) == 2
    joined = read(dirspec)
    @test size(joined.array, 3) == 24
    @test joined.array == read(SourceSpec(E.ERA, "t2m",
                          files = [joinpath(two, "era5_a"),
                              joinpath(two, "era5_b")])).array
    # The dataset-typed readers are the same reads spelled the old way.
    @test read(E.ERA, path, "t2m").array == cr.array
    @test collect(lookup(read(E.ERA, path, "t2m",
                              collect((1:12) .* month_mean_duration)).array,
                         Ti)) == (1:12) .* month_mean_duration
    both = read(E.ERA, dirname(path), "era5_fixture", "t2m",
                [collect((1:12) .* month_mean_duration)])
    @test size(both.array, 3) == 12
    # A source that does not fetch its own files must be told where they are.
    @test_throws "file = path" SourceSpec(E.ERA, "t2m")
    @test_throws "not $(2) of them" SourceSpec(E.ERA, "t2m", file = path,
                                               files = [path])
end

@testset "a dated series carries its policy from the spec to the layer" begin
    E = EcoSISTEM
    path = _erafixture(mktempdir())
    area = StudyArea(regime = SourceSpec(E.ERA, "t2m", file = path,
                                         atend = HoldAtEnd()),
                     cellsize = 45.0°, verbosity = :silent)
    # Real month starts are unevenly spaced, so the default `RepeatAtEnd` cannot derive a period
    # from them and says so.
    @test_throws "evenly spaced" materialise(SourceSpec(E.ERA, "t2m",
                                                        file = path), area)
    layer = materialise(SourceSpec(E.ERA, "t2m", file = path,
                                   atend = HoldAtEnd()), area)
    @test layer.change isa E.SeriesLayerChange
    @test layer.change.atend isa HoldAtEnd
    @test layer.change.calendar isa DatedSeries
    # A volumetric fraction over a 7 cm layer reads as a depth of water - 21 mm - so it is a
    # regime in the axis's canonical unit and, times the cell's area, a supply of cubic metres.
    water = SourceSpec(E.ERA, "swvl1", file = path, atend = HoldAtEnd())
    @test Unitful.unit(eltype(read(water).array)) == cm
    @test ustrip(read(water).array[1, 1, 1]) ≈ 0.3 * 7 rtol=1e-6
    regime = materialise(water, area, role = EcoSISTEM.Condition)
    @test Unitful.unit(eltype(regime.matrix)) == mm
    @test all(v -> isapprox(v, 21.0mm, rtol = 1e-6), regime.matrix)
    supply = materialise(water, area, role = EcoSISTEM.Resource)
    @test supply isa EcoSISTEM.Supply{SoilWaterVolume}
    @test Unitful.unit(eltype(supply.matrix)) == m^3
    @test isapprox(supply.matrix[1, 1],
                   0.021m * EcoSISTEM.getcellareas(m^2, area)[1, 1],
                   rtol = 1e-6)
    # The same policy reaches a derived stack through `in_memory_raster`.
    derived = materialise(in_memory_raster(read(SourceSpec(E.ERA, "t2m",
                                                           file = path)),
                                           axis = Temperature,
                                           atend = ErrorAtEnd()), area)
    @test derived.change.atend isa ErrorAtEnd
end

# A CF netCDF file shaped like a 20CRv3 monthly-means file from the NOAA Physical Sciences
# Laboratory: the five catalogued variables on a 1 degree grid running right round the globe in
# the 0 to 360 convention (a rolled longitude axis is only contiguous when the file is), six
# rows of latitude, three monthly `DateTime`s from 1806 in `hours since 1800-1-1`, PSL's unit
# spellings (`degK`, `kg/m^2/s`, `W/m^2`, `frac.`), and `soilw` on a `level` axis of four layer
# tops in cm. Values encode their position and, for `soilw`, their level, so the roll and the
# level selection can both be checked.
function _twentycrfixture(dir)
    path = joinpath(dir, "twentycr_fixture.nc")
    NCDatasets.NCDataset(path, "c") do ds
        NCDatasets.defDim(ds, "lon", 360)
        NCDatasets.defDim(ds, "lat", 6)
        NCDatasets.defDim(ds, "level", 4)
        NCDatasets.defDim(ds, "time", 3)
        lon = NCDatasets.defVar(ds, "lon", Float32, ("lon",),
                                attrib = ["units" => "degrees_east",
                                    "axis" => "X"])
        lon[:] = 0.0f0:1.0f0:359.0f0
        lat = NCDatasets.defVar(ds, "lat", Float32, ("lat",),
                                attrib = ["units" => "degrees_north",
                                    "axis" => "Y"])
        lat[:] = -3.0f0:1.0f0:2.0f0
        lev = NCDatasets.defVar(ds, "level", Float32, ("level",),
                                attrib = ["units" => "cm", "axis" => "Z",
                                    "coordinate_defines" => "top of layer"])
        lev[:] = Float32[0, 10, 40, 100]
        t = NCDatasets.defVar(ds, "time", Float64, ("time",),
                              attrib = [
                                  "units" => "hours since 1800-1-1 00:00:0.0",
                                  "standard_name" => "time", "axis" => "T"])
        t[:] = Dates.DateTime.(1806, 1:3, 1)
        surface(name, units, f) = begin
            v = NCDatasets.defVar(ds, name, Float32, ("lon", "lat", "time"),
                                  attrib = ["units" => units])
            v[:, :, :] = [Float32(f(i, j, k))
                          for i in 1:360, j in 1:6, k in 1:3]
        end
        surface("air", "degK", (i, j, k) -> 270 + i + 10j + 100k)
        surface("prate", "kg/m^2/s", (i, j, k) -> 2e-5)
        surface("dswrf", "W/m^2", (i, j, k) -> 200 + i)
        surface("uswrf", "W/m^2", (i, j, k) -> 50 + i)
        w = NCDatasets.defVar(ds, "soilw", Float32,
                              ("lon", "lat", "level", "time"),
                              attrib = ["units" => "frac."])
        w[:, :, :, :] = [Float32(0.1 * l)
                         for i in 1:360, j in 1:6, l in 1:4, k in 1:3]
        return nothing
    end
    return path
end

@testset "a 20CRv3 netCDF file reads through its catalogue row" begin
    E = EcoSISTEM
    path = _twentycrfixture(mktempdir())
    # A named file is read in place of the download the row would otherwise fetch, and the
    # catalogue's header facts - 1 degree cells, 0 to 360 longitudes - agree with it.
    air = read(SourceSpec(TwentyCR, "air", file = path))
    @test air isa ClimateRaster{TwentyCR}
    @test size(air.array) == (6, 360, 3)
    @test unit(eltype(air.array)) == K
    @test collect(lookup(air.array, Ti)) == Dates.DateTime.(1806, 1:3, 1)
    # Rolled onto -180 to 180: the columns east of 180 now lead, in ascending order.
    lons = ustrip.(parent(lookup(air.array, X)))
    @test lons == collect(-180.0:1.0:179.0)
    # The value at (lat -3, lon 180) was written at i = 181, j = 1, k = 1.
    @test ustrip(air.array[Y(At(-3.0°)), X(At(-180.0°)), Ti(1)]) ≈
          270 + 181 + 10 + 100
    # A precipitation rate stated as a mass of water per area per second reads as a depth per
    # second by the density of water, and converts to the axis's canonical unit from there.
    rain = read(SourceSpec(TwentyCR, "prate", file = path))
    @test unit(eltype(rain.array)) == m / s
    @test ustrip(rain.array[1, 1, 1]) ≈ 2e-8 rtol=1e-6
    @test uconvert(mm / day, rain.array[1, 1, 1]) ≈ 1.728mm / day rtol=1e-6
    # A flux already per unit time is divided by nothing.
    down = read(SourceSpec(TwentyCR, "dswrf", file = path))
    @test unit(eltype(down.array)) == W / m^2
    @test ustrip(down.array[Y(1), X(At(176.0°)), Ti(1)]) ≈ 200 + 177
    # The soil moisture layer selects the file's top level - the row's vertical extent runs
    # -10 cm to 0 cm, so the level at 0 cm - and reads as a depth of water, the fraction times
    # the 10 cm thickness. Level 1 holds 0.1, so every cell is 1 cm of water.
    soil = read(SourceSpec(TwentyCR, "soilw", file = path))
    @test size(soil.array) == (6, 360, 3)
    @test unit(eltype(soil.array)) == cm
    @test all(v -> isapprox(v, 1.0cm, rtol = 1e-6), soil.array)
    # Selecting a level that the file does not hold, or none on a file that has them, is refused
    # rather than read as a series in time.
    r = Rasters.Raster(path, lazy = true, source = Rasters.NCDsource(),
                       name = :soilw)
    @test_throws "no level was selected" E._selectlevel(r, nothing)
    @test_throws "no level was selected" E._lazyopen(path,
                                                     source = Rasters.NCDsource(),
                                                     name = :soilw)
    @test_throws "none of which" E._selectlevel(r, 5cm)
    @test size(E._selectlevel(r, 0.4m)) == (360, 6, 3)
    @test isnothing(E._openkw(SourceSpec(TwentyCR, "air", file = path)).level)
    @test E._openkw(SourceSpec(TwentyCR, "soilw", file = path)).level == 0cm
    # The whole series materialises on a matching grid as a dated regime, and the soil water
    # as a stock supply.
    area = StudyArea(regime = SourceSpec(TwentyCR, "air", file = path,
                                         atend = HoldAtEnd()),
                     cellsize = 1.0°, verbosity = :silent)
    # A cell-centred global layer reaches half a cell past 180 either way, so the grid keeps the
    # 358 whole-degree columns it covers completely and the two half columns at the antimeridian
    # go; each grid cell holds the layer cell it is labelled by.
    @test E.getgridshape(area) == (6, 358)
    layer = materialise(SourceSpec(TwentyCR, "air", file = path,
                                   atend = HoldAtEnd()), area)
    @test layer.change isa E.SeriesLayerChange
    @test layer.change.calendar isa DatedSeries
    @test ustrip.(parent(lookup(layer.matrix, X)))[1:3] ==
          [-179.0, -178.0, -177.0]
    @test ustrip.(layer.matrix[1, 1:3]) ≈ [270 + i + 10 + 100 for i in 182:184]
    supply = materialise(SourceSpec(TwentyCR, "soilw", file = path,
                                    atend = HoldAtEnd()), area,
                         role = E.Resource)
    @test supply isa E.Supply{SoilWaterVolume}
    @test unit(eltype(supply.matrix)) == m^3
end

@testset "a 20CRv3 layer names its own download" begin
    E = EcoSISTEM
    # Nothing is fetched at construction: the spec's one file is the layer's URL, owned by the
    # source so its cache directory is its own.
    spec = SourceSpec(TwentyCR, "air")
    @test length(spec.files) == 1
    asset = only(spec.files)
    @test asset isa E.CachedAsset && asset.owner === TwentyCR
    @test asset.url ==
          "https://downloads.psl.noaa.gov/Datasets/20thC_ReanV3/Monthlies/2mSI-MO/air.2m.mon.mean.nc"
    @test E.layerinfo(TwentyCR, "soilw").file ==
          "https://downloads.psl.noaa.gov/Datasets/20thC_ReanV3/Monthlies/subsfcSI-MO/soilw.mon.mean.nc"
    @test isnothing(E.layerinfo(E.ERA, "t2m").file)
    @test E.layerinfo(TwentyCR, "soilw").verticalextent == (-10cm, 0cm)
    # A named copy takes precedence over the download.
    @test only(SourceSpec(TwentyCR, "air", file = "air.nc").files) == "air.nc"
end

# The real 20CRv3 files, where a checkout holds them: the directory `ECOSISTEM_TWENTYCR_DIR`
# names, or the dependent package's own download beside this one. Never on a runner, and never
# fetched here - each file is several hundred megabytes.
function twentycrdir()
    haskey(ENV, "RUNNER_OS") && return nothing
    dir = get(ENV, "ECOSISTEM_TWENTYCR_DIR",
              joinpath(dirname(pkgdir(EcoSISTEM)), "Africa_plants", "data",
                       "twentycr", "SI-MO"))
    return isfile(joinpath(dir, "soilw.mon.mean.nc")) ? dir : nothing
end

if !isnothing(twentycrdir())
    @testset "the real 20CRv3 soil moisture file reads as catalogued" begin
        E = EcoSISTEM
        dir = twentycrdir()
        spec = SourceSpec(TwentyCR, "soilw",
                          file = joinpath(dir, "soilw.mon.mean.nc"),
                          atend = HoldAtEnd())
        # The header check passes against the row, and the level view keeps the read to one
        # layer of the four.
        cr = read(spec)
        @test size(cr.array) == (181, 360, 2520)
        @test unit(eltype(cr.array)) == cm
        ts = lookup(cr.array, Ti)
        @test first(ts) == Dates.DateTime(1806, 1, 1) &&
              last(ts) == Dates.DateTime(2015, 12, 1)
        lons = ustrip.(parent(lookup(cr.array, X)))
        @test first(lons) == -180.0 && last(lons) == 179.0
        # A fraction of at most one over a 10 cm layer is at most 10 cm of water; the sea is
        # missing and reads as NaN.
        finite = filter(isfinite, ustrip.(cr.array[:, :, 1]))
        @test !isempty(finite) && maximum(finite) <= 10.05 &&
              minimum(finite) >= 0
    end
end

@testset "CF unit spellings" begin
    E = EcoSISTEM
    @test E._parsecfunit("J m**-2") == u"J*m^-2"
    @test E._parsecfunit("m**3 m**-3") == u"m^3*m^-3"
    @test E._parsecfunit("") == NoUnits
    # The two spellings NOAA PSL uses that Unitful does not: `degK` and `Kg`. Named here because
    # the alias table is what makes them parse; nothing else does.
    @test E._parsecfunit("degK") == u"K"
    @test E._parsecfunit("Kg/m^2/s") == u"kg/m^2/s"
    @test E._parsecfunit("kg/m^2/s") == u"kg/m^2/s"
    @test E._parsecfunit("frac.") == NoUnits
end

@testset "a blank CRS is treated as an absent one" begin
    CP = EcoSISTEM
    blank = Rasters.GeoFormatTypes.WellKnownText(Rasters.GeoFormatTypes.CRS(),
                                                 "")

    @test CP._isblankcrs(blank)
    @test CP._isblankcrs(Rasters.GeoFormatTypes.WellKnownText(Rasters.GeoFormatTypes.CRS(),
                                                              "   "))
    # A real CRS is not blank, whatever form it arrives in.
    @test !CP._isblankcrs(Rasters.EPSG(4326))
    @test !CP._isblankcrs(Rasters.GeoFormatTypes.WellKnownText(Rasters.GeoFormatTypes.CRS(),
                                                               ArchGDAL.toWKT(ArchGDAL.importEPSG(4326))))

    # ...so a blank CRS gets the same WGS84 assumption an absent one does, rather than throwing.
    @test CP._crsunit(blank) == CP._crsunit(nothing) == °
    @test CP._crsunit(Rasters.EPSG(27700)) == u"m"

    # And the whole read normalises it away, so nothing downstream ever sees the blank text: a
    # CRS-less file comes back on ° coordinates carrying `nothing`, which `_samecrs`/`_dimsextent`
    # already handle but `importCRS` would not.
    path = joinpath(mktempdir(), "nocrs.tif")
    ArchGDAL.create(path, driver = ArchGDAL.getdriver("GTiff"), width = 4,
                    height = 4, nbands = 1,
                    dtype = Float32) do ds
        ArchGDAL.write!(ds, fill(1.0f0, 4, 4), 1)
        return ArchGDAL.setgeotransform!(ds, [0.0, 1.0, 0.0, 4.0, 0.0, -1.0])
    end
    @test CP._isblankcrs(Rasters.crs(Rasters.Raster(path)))
    r = read(RasterFileSpec(path, axis = EcoSISTEM.NicheAxis))
    @test r isa EcoSISTEM.ClimateRaster{EcoSISTEM.SyntheticData}
    # A named source is recorded as given - provenance, not a claim the file is a layer of that
    # dataset, so the dataset's catalogue row is not checked against it.
    @test read(RasterFileSpec(path, axis = EcoSISTEM.NicheAxis,
                              source = WorldClim{BioClim})) isa
          EcoSISTEM.ClimateRaster{WorldClim{BioClim}}
    a = r.array
    @test size(a) == (4, 4)
    @test unit(eltype(parent(DimensionalData.lookup(a, Y)))) == °
    @test isnothing(Rasters.crs(DimensionalData.dims(a, Y)))
end

end
