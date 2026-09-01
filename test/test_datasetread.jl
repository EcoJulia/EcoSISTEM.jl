# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDatasetread

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources
using Rasters
using ArchGDAL
using DimensionalData
using Statistics
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
        @test_nowarn read(WorldClim{Climate}, :wind, month = 1:12)
        @test_nowarn read(CRUTS, winddir, "tavg")
        @test_nowarn read(CHELSA{Climate}, winddir, "wind")
        # Every whole-dataset read here is downsampled, and the reason is a hard ceiling rather
        # than tidiness: a GitHub runner has 16 GB, and one of these files reading at native
        # resolution took the process to 7.0 GB and the runner to a shutdown signal. Measured peak
        # RSS for a fresh process, baseline 1.0 GB: CHELSA bioclim is a 43200×20880 global grid and
        # allocates several ~7 GiB Float64 arrays whole; `EarthEnv{LandCover}` costs 2.2 GB over
        # baseline at its default scale against 0.3 GB at 40, and is read twice in this file;
        # `WorldClim{BioClim}` costs 1.5 GB against 0.3 GB at 4.
        #
        # `scale` is safe for what these assert -- that the read emits no warning, and that what
        # comes back is unitless -- since neither is a property of the resolution. A test that does
        # depend on the grid asks for one layer rather than a whole dataset, so it never reaches
        # this size.
        bigrasters() && @test_nowarn read(CHELSA{BioClim}, 1, scale = 20)
        @test_nowarn read(EarthEnv{LandCover}, scale = 40)
        @test_nowarn readfile(bio1)
    end

    @testset "Output data" begin
        bioclim = read(WorldClim{BioClim}, scale = 4)
        cr = read(CRUTS, winddir, "tavg")
        rf = readfile(bio1)

        @test unit(bioclim.array[1]) == unit(rf[1]) == NoUnits
        if bigrasters()
            ch_b = read(CHELSA{BioClim}, 1, scale = 20)
            @test unit(ch_b.array[1]) == NoUnits
        end
    end

    @testset "Output data 2" begin
        landcover = read(EarthEnv{LandCover}, scale = 40)
        @test unit(landcover.array[1]) == NoUnits
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
    # of unrelated bands - and so whether it can drive a layer through time at all.
    #
    # `CHELSA{Climate}` had no method and fell through to the `Dim{:layer}` fallback, so a
    # `SourceSpec(CHELSA{Climate}, ...)` read through the generic path was **not recognised as a
    # series**. It failed silently: twelve monthly files loaded fine, simply as twelve bands. This is
    # asserted against its siblings rather than alone, because the bug was an omission - a test that
    # only checked CHELSA would not have said whether `Ti` was right or the fallback was wrong.
    @testset "monthly sources stack on time, band sources on layers" begin
        @test EcoSISTEM._stackaxis(WorldClim{Climate}) == Ti
        @test EcoSISTEM._stackaxis(CHELSA{Climate}) == Ti
        # Only the `Climate` layers are monthly. A source's *bioclim* variables are one file per
        # variable, so they must keep stacking on `Dim{:layer}`.
        @test EcoSISTEM._stackaxis(CHELSA{BioClim}) == Dim{:layer}
        @test EcoSISTEM._stackaxis(WorldClim{BioClim}) ==
              Dim{:layer}
        @test EcoSISTEM._stackaxis(EarthEnv{LandCover}) ==
              Dim{:layer}
    end

    @testset "a partial monthly read knows which months it holds" begin
        full = read(WorldClim{Climate}, :wind, month = 1:12)
        part = read(WorldClim{Climate}, :wind, month = 2:4)

        @test collect(DimensionalData.lookup(full.array, Ti)) ==
              (1:12) .* month_mean_duration
        @test collect(DimensionalData.lookup(part.array, Ti)) ==
              (2:4) .* month_mean_duration

        # An uneven request stays uneven, which is what lets a series hold each slice until the next
        # rather than pretending the months are consecutive.
        sparse = read(WorldClim{Climate}, :wind, month = [1, 6, 12])
        @test collect(DimensionalData.lookup(sparse.array, Ti)) ==
              [1, 6, 12] .* month_mean_duration

        # A single month has no series in it, so it stays 2-D and carries no time axis at all -
        # a deliberate carve-out, since a length-1 `Ti` would change the result's dimensionality and
        # every `ndims == 2` static-vs-series branch downstream with it.
        one = read(WorldClim{Climate}, :wind, month = 2)
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
    # EarthEnv is the only shipped source with `_defaultscale` > 1, which is what made it look
    # source-specific rather than a general consequence of aggregating. Guarded because it needs the
    # real file: the vector-lookup condition that triggers it cannot be reproduced synthetically.
    @testset "a coarsened read can be cut (Regular span with real bounds)" begin
        L = DimensionalData.Lookups
        scotland = EcoSISTEM.boundingbox("Scotland", islands = true)
        whole = read(EarthEnv{LandCover}, 7)
        for d in (Y, X)
            @test L.span(dims(whole.array, d)) isa L.Regular
            @test all(!isnothing, L.bounds(dims(whole.array, d)))
        end
        # The end the fix exists for - this threw a MethodError before it.
        cut = read(EarthEnv{LandCover}, 7, cut = scotland)
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
        sources = Any[(read(WorldClim{BioClim}, :bio1), (-180°, 180°),
                       (-90°, 90°)),
                      (read(EarthEnv{LandCover}, 7), (-180°, 180°),
                       (-56°, 90°))]
        bigrasters() && push!(sources,
              (read(CHELSA{BioClim}, 1, scale = 20),
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
    # reintroduced exactly where it was just removed. EarthEnv (`_defaultscale` 10) is the only
    # shipped source that exercises it.
    @testset "a windowed read equals the whole read cropped" begin
        CP = EcoSISTEM
        scot = CP.boundingbox("Scotland", islands = true)
        for (src, code) in ((WorldClim{BioClim}, :bio1),   # scale 1
            (EarthEnv{LandCover}, 7))      # scale 10 - block alignment
            whole = EcoSISTEM._read(SourceSpec(src, code))
            windowed = EcoSISTEM._read(SourceSpec(src, code, cut = scot))
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
              EcoSISTEM._rastercrs(EcoSISTEM._read(SourceSpec(WorldClim{BioClim},
                                                              :bio1)))
        # Read keywords a `SourceSpec` may carry are accepted and ignored - `scale` cannot alter a CRS.
        @test EcoSISTEM.sourcecrs(WorldClim{BioClim}, :bio1,
                                  scale = 20, cut = nothing) ==
              EcoSISTEM.sourcecrs(WorldClim{BioClim}, :bio1)
    end

    @testset "Output data 3" begin
        cr = read(CRUTS, winddir, "tavg")
        worldclim = read(WorldClim{Climate}, :wind)
        ch_m = read(CHELSA{Climate}, winddir, "wind")

        # the directory readers self-attach the actual unit (temperature now in °C, no hidden K shift)
        @test unit(cr.array[1]) == °C
        @test unit(ch_m.array[1]) == m / s
        # `read` returns bare magnitudes - the unit lives in the layer table (`layerunit`)
        @test unit(worldclim.array[1]) == NoUnits
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
    a = readfile(path)
    @test size(a) == (4, 4)
    @test unit(eltype(parent(DimensionalData.lookup(a, Y)))) == °
    @test isnothing(Rasters.crs(DimensionalData.dims(a, Y)))
end

end
