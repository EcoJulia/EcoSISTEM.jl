# SPDX-License-Identifier: LGPL-3.0-or-later
#
# READING A CATALOGUED DATASET. The pipeline between a set of files on disk and a `ClimateRaster`
# with its coordinates, units and layer axis in place:
#
# axes      `_stackaxis`, `_mkstackaxis`, `_stackcoords`, `_stackdim`, `_spatialdim` -- what the
# extra dimension of a multi-layer or monthly read is, and how it is labelled
# geometry  `_snaparcsec`, `_snaporigin`, `_sourcestep`, `_locus` -- putting a file's coordinates
# onto the lattice its provider intends
# reading   `_readraster`, `_readbytes`, `_readsource`, `_readmultilayer`, `_readraw`, `readfile`
# caching   `_cachedlayer`, `_aggcachepath`, `_fnid`
#
# It sits in the main module rather than in an extension because it names no climate-data-source
# package at all: every function here takes a type `T` and file paths,
# and the dataset-specific methods (`_defaultscale(::Type{<:EarthEnv{<:LandCover}}) = 10` and the
# rest) live in `EcoSISTEMRasterDataSourcesExt`, which supplies them over these generic helpers.
# Putting the generic half in an extension would have made an extension that names nothing from its
# own trigger, which stops precompiling the moment that dependency is weakened.
#
# `ClimatePref` keeps what genuinely is about climate data: the `ERA`/`CERA` NetCDF readers,
# `boundingbox`, and its deprecations. It imports `_applycut`, `_locus` and `_rastertodimarray` back
# from here, so the dependency runs one way.

using Unitful

using Unitful.DefaultSymbols

using EcoSISTEM.Units

using Dates: Dates

using CSV

using Rasters

using DimensionalData

using DimensionalData.Lookups: Sampled, Categorical, Intervals, Center,
                               locus, ForwardOrdered, Irregular, Regular

import Rasters: Projected

using JLD2: jldsave, jldopen

import Rasters: X, Y, Ti

import Unitful.°, Unitful.°C, Unitful.mm

# ArchGDAL is needed to register the GDAL backend Rasters uses to read GeoTIFFs.
import ArchGDAL

import Base.read

# How far a step may sit from a whole arcsecond and still be taken as one. The worst real offender is
# CHELSA's 89.999964 arcsec for 90 - a relative error of 4e-7 - while the smallest *genuine* step
# difference anyone could mean is a whole arcsecond out of the finest grids we see (1 in 30, i.e.
# 3e-2). Five orders of magnitude of clear water either side, so this discriminates rounding noise
# from real structure rather than papering over a comparison.
const _ARCSEC_RTOL = 1e-5

# How far an origin may sit from the cell lattice and still be taken as belonging to it, as a fraction
# of one cell. The worst real offender is CHELSA, whose `bio1` origins are ~1.7% of a cell out; the
# smallest *deliberate* offset anyone uses is a half cell (centre- vs edge-registration, 50%). An
# order of magnitude of clear water either side.
const _ORIGIN_CELLFRAC = 0.05

# A raster is read whole (then aggregated in memory) when it fits in this fraction of *total* system
# RAM; anything larger (e.g. the multi-GB CHELSA bioclim file on a small machine) stays lazy.
const _READ_WHOLE_FRACTION = 0.5

# Content hash of this file, folded into the aggregate cache key so any change to the reading /
# aggregation machinery here invalidates the cache - a cached result is only valid for the code that
# produced it. Evaluated at precompile, so it tracks the source automatically (a change to this file
# recompiles the module and updates the hash).
const _AGGCODEHASH = hash(read(@__FILE__, String))

"""
    readfile(file::String; cut = nothing)

Import a raster file from a path string into a `DimArray`.

# Arguments

  - `file`: path to the raster, in any format GDAL reads.
  - `cut`: a region to restrict the read to as an `Extents.Extent` of `°` intervals - from
    [`boundingbox`](@ref), say - or `nothing` for the whole file. Applied **lazily**, before the
    pixels are fetched, so cutting a global layer to one country costs the country rather than the
    globe.
  - `xmin`, `xmax`, `ymin`, `ymax`: **deprecated**, and warn when used. All four together mean
    `cut = Extent(Y = (ymin, ymax), X = (xmin, xmax))`; passing some of them, or passing them
    alongside `cut`, is an error. Pass `cut` instead.
"""
function readfile(file::String; cut = nothing, xmin = nothing, xmax = nothing,
                  ymin = nothing, ymax = nothing)
    n = count(!isnothing, (xmin, xmax, ymin, ymax))
    if n == 4
        isnothing(cut) ||
            error("`readfile`: pass either `cut` or the `xmin`/`xmax`/`ymin`/`ymax` extent, not both.")
        Base.depwarn("The `xmin`/`xmax`/`ymin`/`ymax` keywords are deprecated; they are being used " *
                     "as `cut = Extent(Y = (ymin, ymax), X = (xmin, xmax))`. Pass `cut` (e.g. from " *
                     "`boundingbox`) instead.", :readfile)
        cut = Extents.Extent(Y = (ymin, ymax), X = (xmin, xmax))
    elseif n != 0
        error("`readfile` needs all four of `xmin`/`xmax`/`ymin`/`ymax` (deprecated) or none of " *
              "them; got $n.")
    end
    return _applycut(_rastertodimarray(_readraster(file)), cut)
end

# Whether a Rasters dimension is one of the two spatial axes, so `_rastertodimarray` can
# filter to the non-spatial (band/time) ones.
_isspatialdim(::Union{X, Y}) = true

_isspatialdim(::Rasters.Dimension) = false

# The sampling locus (`Center`/`Start`/`End`) to build a Y/X lookup with: whatever the source
# raster's own dimension declares, preserved rather than assumed - confirmed real GDAL files vary
# (the NatureScot HLCM/BNG files declare `Start`; a WGS84 GeoTIFF commonly declares plain
# `Points`, i.e. no locus at all). A locus mismatch silently shifts every cell's implied bounds by
# half a pixel, which matters for crop/selection accuracy at the edges.
#
# **`Points` is refused rather than defaulted**, and the refusal is the honest form: a point
# coordinate has no extent, so there is no cell for it to label, and falling back to `Center` would
# be exactly the assumption it looks like it avoids. No shipped source reaches this - every one
# declares `Intervals` - so it guards data yet to arrive.
function _locus(sourcedim)
    s = DimensionalData.sampling(sourcedim)
    s isa Intervals ||
        error("this raster's `$(nameof(typeof(sourcedim)))` axis declares `$(nameof(typeof(s)))` " *
              "sampling, but a spatial axis must declare `Intervals` - a point coordinate has no " *
              "extent, so there is no cell for it to label and no honest way to guess one. " *
              "Re-declare the source's sampling, and its locus, to say what ground each " *
              "coordinate covers.")
    return locus(s)
end

# Kind of axis stacked when combining multiple files of a single layer/variable: a monthly series
# (`Ti`) or a band/variable index (`Dim{:layer}`). Returns the *type itself*, not a `:time`/
# `:layer` Symbol - so `_mkstackaxis` below is genuine multiple dispatch, resolved at compile
# time from the statically-known source type `T`, rather than a runtime branch on a value (the
# `Val(Symbol(...))`-shaped anti-pattern this project avoids, just the value-vs-dispatch mirror of
# it: a compile-time-knowable choice was being funnelled through a runtime `Symbol` comparison).
# The per-dataset methods (WorldClim's and CHELSA's monthly climate stack on `Ti`) are in
# `EcoSISTEMRasterDataSourcesExt`; this is the fallback they specialise.
_stackaxis(::Type) = Dim{:layer}

# Build the axis stacked when combining `n` files of a single layer/variable - one method per
# axis type (see `_stackaxis`), not a runtime branch on a value.
_mkstackaxis(::Type{Ti}, n::Integer) = Ti((1:n) .* month_mean_duration)

# **`Categorical`, matching what a real multi-layer read builds.** `_stacklayers` labels a real read
# `Dim{:layer}(layernames)`, a `Vector{Symbol}`, which DimensionalData formats as `Categorical`. A
# bare `1:n` here would be formatted `Sampled`/`Points` instead, making the axis a **different kind**
# depending on provenance: `code = 2` would resolve as a *label* against a real read and as a
# *coordinate on a numeric axis* against this fallback - indistinguishable only because label and
# position coincide when the labels are `1:n`. Declaring it `Categorical` makes a layer code mean one
# thing everywhere.
_mkstackaxis(::Type{Dim{:layer}}, n::Integer) = Dim{:layer}(Categorical(1:n))

# The coordinates of a stacked axis, given the *slices actually read* rather than just how many.
#
# A monthly `Ti` is labelled with the months it really holds, so `month = 2:4` reads as
# February-April rather than as the first three months. Falling back to `1:n` there is not a
# harmless mislabel: the plot recipe looks a month up by coordinate (`At(ind * month_mean_duration)`),
# so a partial read handed `February` the *next* month's grid under February's own name.
#
# Only a time axis has this identity - a band index genuinely is `1:n`, so `Dim{:layer}` and any
# future axis take the count-only form. Dispatched on the axis type, like `_mkstackaxis` itself.
#
# Uneven requests stay uneven (`month = [1, 6, 12]` really is spaced 1, 6, 12), which is what lets
# a series hold each slice until the next and lets `RepeatAtEnd` refuse a series with no derivable
# period - a check that a relabelled-to-`1:n` axis silently passed.
function _stackcoords(::Type{Ti}, n::Integer, slices)
    length(slices) == n ||
        error("read $n file(s) but was asked for $(length(slices)) month(s) ($slices) - the two " *
              "must agree, or the time axis would name each slice after the wrong month.")
    return Ti(collect(slices) .* month_mean_duration)
end

# `(Ti, n, nothing)` needs its own method, or it is genuinely **ambiguous**: the `Ti` method above
# is more specific in the first argument and the `::Nothing` one below in the third, so neither wins.
# A caller reaching `_readsource` without going through `read` - the deprecated `readworldclim`, which
# is handed its own file list - passes exactly that, and an ambiguity is a call-time error, so it
# surfaced only in the suite rather than at load.
_stackcoords(::Type{Ti}, n::Integer, ::Nothing) = _mkstackaxis(Ti, n)

_stackcoords(A, n::Integer, ::Nothing) = _mkstackaxis(A, n)

_stackcoords(A, n::Integer, _) = _mkstackaxis(A, n)

# The axis stacked onto a `_rastertodimarray` result: preserve the source raster's own
# non-spatial dimension when it carries real information worth keeping - confirmed empirically
# that `Rasters`/`NCDatasets` already decodes a netCDF file's CF `units`/`calendar` time metadata
# into genuine `Dates.DateTime` values (not approximated fixed-duration ticks) when read via
# `source = NCDsource()`, so that must not be thrown away and rebuilt as a synthetic
# `(1:n)*month_mean_duration` sequence. Falls back to a synthetic sequential axis (`_mkstackaxis`) when there is
# no source dimension (`nothing`, e.g. stacking bare 2-D files) or it carries no real information
# of its own (e.g. a plain GDAL `Band` index, itself just `1:n`) - dispatched on the source
# dimension's own type (a `Ti` genuinely carrying `DateTime`s vs. anything else), not an `isa`
# check, since which case applies is only known from the file's actual runtime structure.
#
# Tested on `Dates.TimeType`, never `Dates.AbstractDateTime` - see `_istimeaxis`. The narrower name
# would discard a real `Date` axis and every non-Gregorian `CFTime` one, rebuilding a synthetic
# `1:n` in their place, which is the silent version of the failure this branch exists to prevent.
_stackdim(::Nothing, axistype, n) = _mkstackaxis(axistype, n)

function _stackdim(sourcedim::Ti, axistype, n)
    return eltype(sourcedim.val) <: Dates.TimeType ? sourcedim :
           _mkstackaxis(axistype, n)
end

_stackdim(_, axistype, n) = _mkstackaxis(axistype, n)

# Convert a Rasters raster (2-D `X`/`Y`, or 3-D `X`/`Y` with a band/time dimension) into the
# `DimArray` EcoSISTEM uses internally: dim 1 = `Y` (latitude/northing, ascending), dim 2 = `X`
# (longitude/easting, ascending), and an optional dim 3 = `Ti` or `Dim{:layer}` (see `_stackdim`).
# Coordinates come from the raster's own `Y`/`X` lookups - `Rasters.jl` already computes their
# `Regular` span from the file's geotransform, so it's kept as-is (via
# `parent(Rasters.lookup(...))`, not `collect`-ed into a bare `Vector`) rather than rebuilt by
# hand. The CRS is preserved (a `Projected` lookup, not stripped to a bare `Sampled`), and the
# axis unit (`_crsunit`) and sampling locus (`_locus`) are read from the source's own CRS/
# metadata rather than assumed `°`/`Center` - a British National Grid file's axes are in metres
# with a `Start` locus, not degrees with a `Center` locus, and reading it as though it were WGS84
# would silently mislabel every coordinate. `unit` (the *data* unit, not the axis unit) is
# attached; pass `NoUnits` for a dimensionless layer.
#
# `:latitude` is `Y` and `:longitude` is `X`, which is the geographically correct pairing and the
# order used throughout the package.
# One spatial dimension of a read raster, unit-tagged and with an explicitly **`Regular`** span.
#
# The span must be stated, not inferred. `parent(lookup(...))` is a `Vector`, and `AutoSpan` infers
# `Irregular((nothing, nothing))` from a vector - a span with *no bounds*. Anything that needs real
# bounds then fails: `_applycut`'s `Touches` selector compares `nothing < 60.86°` and throws a bare
# `MethodError`, which is why `read(EarthEnv{LandCover}, ..., cut = ...)` was broken and the Scotland
# example had to crop after the read instead. `Rasters.aggregate` is what drops the span in the first
# place, so a coarsened read hits this and a plain one does not, making it look source-specific.
#
# Recomputing from the values is safe and is what `_striplength` already does for the same reason:
# every raster this codebase reads is a regular grid, so a step taken from the first two coordinates
# describes all of them.
#
# The axis is **rebuilt** from `(origin, step, length)` rather than passed through as the file's
# coordinate vector. Two reasons. The vector accumulates rounding across its length - 43200 additions
# of an inexact step - so `vals[end]` is not where the grid says it ends; a range is exact by
# construction from three numbers. And it is the natural place to correct a source's own metadata
# errors, which `_snaparcsec`/`_snaporigin` do for angular axes.
function _spatialdim(D, r, axisunit, crs)
    vals = parent(Rasters.lookup(r, D)) .* axisunit
    step = _sourcestep(r, D, axisunit, vals)
    lattice = range(_snaporigin(first(vals), step), step = step,
                    length = length(vals))
    return D(Projected(lattice, sampling = Intervals(_locus(dims(r, D))),
                       crs = crs, order = DimensionalData.order(dims(r, D)),
                       span = Regular(step)))
end

# The cell step to declare: the source's **own** where it states one, and only otherwise derived from
# the coordinates.
#
# Preferring the source's matters beyond tidiness. Differencing the first two coordinates gives a
# subtly different `Float64` depending on *where* in the grid you difference - a global WorldClim read
# yields 0.1666666666666714° and the same layer cut to Scotland 0.1666666666666643°, because
# `(a + 2s) - (a + s)` rounds differently for different `a`. Two reads of one file would then disagree
# on their resolution, and the grids built from them would differ by a cell here and there. The
# geotransform's own step has no such ambiguity. Only an aggregated read (whose span Rasters drops)
# needs the fallback.
function _sourcestep(r, D, axisunit, vals)
    sp = DimensionalData.Lookups.span(dims(r, D))
    sp isa Regular && return _snaparcsec(sp.step * axisunit)
    return _snaparcsec(length(vals) > 1 ? vals[2] - vals[1] :
                       oneunit(eltype(vals)))
end

# A step snapped to a whole arcsecond, where that is plainly what was meant.
#
# Every geographic grid this package reads is *intended* to sit on a whole-arcsecond lattice - 43200
# columns of "0.0083333333°" is exactly 360° of 30 arcsec cells - but some files record it badly:
# CHELSA writes 30 arcsec as 29.99999988 and 90 as 89.999964. Left alone those make every resolution
# question a matter of tolerance (0.025 / 0.0083333333 is 3.0000000048, not 3); snapped, they divide
# exactly, so "is this grid a whole multiple of that one?" has a straight yes/no answer.
#
# Two deliberate limits. Only *angular* axes: a projected CRS's step is a length, where arcseconds
# mean nothing. And only where the value is already within `_ARCSEC_RTOL` of a whole arcsecond, so a
# genuinely finer-than-arcsecond or non-integral grid is left exactly as found rather than silently
# moved - the `n == 0` case covers sub-arcsecond steps, which must never round to nothing.
#
# `round(::Quantity)` without a target type is a trap for angles - Unitful rounds via radians, so
# `round(600.0arcsecond)` is `0.0`, a bare `Float64`. The unit must be named.
#
# `_arcsecs` (`StudyArea.jl`) is what recognises the result, and the two are only safe as a pair.
# A whole number of arcseconds has no exact `Float64` degree representation, and *which* neighbouring
# value you land on depends on the route: `uconvert` scales by a precomputed factor and rounds twice,
# a hand-written `n / 3600` rounds once, and the two differ in the last bit for 29 of the first 3600
# counts. `_arcsecs` therefore accepts anything within a few ULP rather than testing equality against
# one particular reconstruction - see its comment. `test_datasetread.jl` pins the round trip for all
# 3600 counts by both routes; changing either function without re-running that sweep is how the
# integer-arcsecond path would silently switch off.
function _snaparcsec(step)
    _isangle(step) || return step
    n = round(arcsecond, step)
    (!iszero(n) && isapprox(step, n, rtol = _ARCSEC_RTOL)) || return step
    # Returned in the unit it arrived in, so a caller's choice of `°`/`′`/`″` survives the snap.
    return uconvert(unit(step), n)
end

# An origin snapped onto the lattice its own cell size implies, where it plainly belongs there.
#
# Geographic grids are anchored on the equator and the prime meridian, so a grid of `step`-sized cells
# should start at a whole multiple of `step` - CHELSA's `bio1` is meant to run -180°...180° by 84°...-90°,
# and with the step snapped those extents come out exact. What its file actually records is ~0.5
# arcsec adrift on each axis, which is 1.7% of a 30 arcsec cell (~15 m of a ~900 m cell).
#
# This corrects the *extent* as well as the origin: `_spatialdim` builds the axis as
# `range(origin, step, length)`, so a snapped origin and an exact step put both ends of the grid where
# the source intends them, rather than wherever accumulated rounding lands.
#
# Snap to a multiple of the **step**, not to a whole arcsecond. That distinction is the whole
# correction: `bio1`'s Y origin is 302399.4975 arcsec, which rounds to *83.999722°* at arcsecond
# granularity - the wrong way, and not a value anyone intended - but to exactly **84°** on its own 30
# arcsec lattice. Testing against the arcsecond lattice is what previously made these origins look
# irretrievably ambiguous.
#
# Angular axes only, and only within `_ORIGIN_CELLFRAC` of the lattice. A projected grid has no
# comparable anchor (and no representation error to fix - whole metres are exact in `Float64`), and a
# grid genuinely registered off-lattice must keep the offset it was given. CHELSA's `clt` stays odd
# whatever we do - 14401 columns where 14400 would be global, so 360.025° wide - but snapping at least
# makes that a visible, reproducible overshoot instead of float noise.
function _snaporigin(origin, step)
    (_isangle(origin) && _isangle(step)) || return origin
    cells = NoUnits(origin / step)      # a count of cells, so genuinely dimensionless
    n = round(cells)
    abs(cells - n) < _ORIGIN_CELLFRAC || return origin
    return n * step
end

# A `Rasters.Raster` as the plain `DimArray` this package stores, with missing values as `NaN`,
# latitude ascending and the unit attached. This is the one boundary where Rasters' representation is
# converted to ours, so a change in what a layer holds is a change here and nowhere else.
function _rastertodimarray(ras::Rasters.AbstractRaster; unit = NoUnits,
                           stackaxis::Type = Ti)
    r = Rasters.replace_missing(ras, NaN)
    # latitude/northing ascending (GeoTIFF `Y` is usually stored north->south)
    (first(Rasters.lookup(r, Y)) < last(Rasters.lookup(r, Y))) ||
        (r = reverse(r, dims = Y))
    # reorder axes to (Y, X[, band/time])
    others = filter(!_isspatialdim, Rasters.dims(r))
    r = permutedims(r, (Y, X, others...))
    # Normalised here, at the one point a file's CRS enters the package, rather than defended against
    # at each use: a blank CRS *is* an absent one, and every downstream consumer (`_samecrs`,
    # `_dimsextent`, `_crsunit`, `_isprojectedcrs`) already handles `nothing`, while none of them
    # survives `ArchGDAL.importCRS("")`.
    crs = _isblankcrs(Rasters.crs(r)) ? nothing : Rasters.crs(r)
    axisunit = _crsunit(crs)
    ydim = _spatialdim(Y, r, axisunit, crs)
    xdim = _spatialdim(X, r, axisunit, crs)
    data = Float64.(Array(r)) .* unit
    ndims(data) == 2 && return DimArray(data, (ydim, xdim))
    stackdim = _stackdim(isempty(others) ? nothing : first(others), stackaxis,
                         size(data, 3))
    return DimArray(data, (ydim, xdim, stackdim))
end

# Apply the optional `cut` (an `Extents.Extent` of `°` bounds) to a lat×long[×z] world.
# `Touches`, not `Between`/`..`: a crop must be overlap-inclusive (every cell touching the box is
# kept, e.g. cropping to a region's bounding box must not lose its edge cells) - confirmed
# empirically that `Between` under `Intervals` sampling instead requires the *whole* cell to lie
# inside the box, silently dropping edge cells. Indexing only `Y`/`X` leaves any third dimension
# (time/band) untouched, so no separate 2-D/3-D branch is needed.
_applycut(world, ::Nothing) = world

function _applycut(world, cut)
    return world[Y(Touches(cut.Y[1], cut.Y[2])),
                 X(Touches(cut.X[1], cut.X[2]))]
end

# Estimated bytes to hold a raster whole in memory. Rasters' own guard compares against *free*
# memory, which is over-conservative - idle memory is reclaimable / pageable - so we instead gate
# on a fraction of *total* RAM (below).
function _readbytes(r)
    elbytes = try
        Base.aligned_sizeof(eltype(r))   # handles isbits `Union{Missing, T}` (selector byte incl.)
    catch
        sizeof(Float64)                  # conservative fallback
    end
    return prod(size(r)) * elbytes
end

# Stable identifier for the aggregation reducer `fn`, used in the aggregate cache key. Anonymous
# closures (whose `nameof` starts with `#`, and isn't stable across sessions) return `nothing`, so
# they are not cached.
function _fnid(fn)
    n = string(nameof(fn))
    return startswith(n, "#") ? nothing : Symbol(n)
end

# Path in EcoSISTEM's own scratch subdir (`assetdir()`, not RasterDataSources') where the
# `scale`-aggregated `unit`-tagged form of source file `f` (with reducer `fn`) is cached as a JLD2
# `DimArray`, or `nothing` if `fn` is not cacheable. Keyed on the source's path/size/mtime - a
# `stat`, deliberately not a read, so a cache hit never touches the multi-GB source - plus `scale`,
# the reducer id, the unit and `_AGGCODEHASH` (so a machinery change invalidates it).
function _aggcachepath(f, scale, fn, u)
    id = _fnid(fn)
    isnothing(id) && return nothing
    key = string(hash((abspath(f), filesize(f), mtime(f), scale, id, string(u),
                       _AGGCODEHASH)), base = 16)
    return joinpath(assetdir(), "aggregates", key * ".jld2")
end

# Mask raw integer-band fill sentinels. GDAL sources (e.g. CHELSA) store a fill at the raw band's
# `typemax` (or `typemin` for signed types) that the file's *declared* nodata may not capture - CHELSA
# declares `-99999`, unrepresentable in its `UInt` bands, so the real fill (`0xffffffff` etc.) survives the
# default scaled read and becomes a spurious huge value (gsp's `0xffffffff` -> 4.29e8). Locate those cells in
# the raw band `raw` and set them `missing` on the scaled raster `r`, lazily (so a disk-backed raster stays
# lazy). Only integer bands carry such sentinels; float bands rely on the file's declared nodata and pass
# through unchanged.
function _mask_int_fills(r, raw)
    T = nonmissingtype(eltype(raw))
    T <: Integer || return r
    fills = T <: Signed ? (typemin(T), typemax(T)) : (typemax(T),)
    return Rasters.rebuild(r,
                           broadcast((value,
                                      rawvalue) -> (!ismissing(rawvalue) &&
                                                    rawvalue in fills) ?
                                                   missing : value, r, raw))
end

# The index range of axis `d` that `[lo, hi]` touches, widened outward to whole `scale`-sized
# aggregation blocks - or `nothing` when the box selects nothing at all.
#
# The widening is the whole difficulty. `Rasters.aggregate` blocks from index 1, so cropping a
# raster before aggregating it moves every block boundary and the coarsened cells land somewhere
# else - read-extent variance reintroduced at the very point it was just eliminated. Starting the
# crop on a block boundary (`i ≡ 1 mod scale`) makes the crop's blocks exactly the uncropped ones, so
# the aggregated coordinates are identical whatever window was asked for. The end is clamped to the
# axis length, which leaves the same partial final block a whole-file read would have had.
function _blockrange(d, lo, hi, scale::Integer)
    idx = DimensionalData.Lookups.selectindices(DimensionalData.lookup(d),
                                                Rasters.Touches(lo, hi))
    isempty(idx) && return nothing
    i, j = first(idx), last(idx)
    scale > 1 || return i:j
    return (1 + scale * fld(i - 1, scale)):min(length(d), scale * cld(j, scale))
end

# Crop a *lazy* raster to `cut` before anything materialises, so only the window is read off disk.
#
# This is purely an optimisation, and deliberately so: `_readsource` still applies `_applycut` to the
# assembled result, which is idempotent (every surviving cell touches the box, so re-selecting keeps
# them all). The lazy crop may therefore be conservative without changing any answer - a property
# worth keeping, given how much trouble read-extent variance has caused.
#
# Skipped for a projected source: `cut` is a `LatLong` in degrees and such a raster's coordinates
# are metres, so there is nothing to compare without reprojecting. That matches the existing
# post-read `_applycut`, which would reject the same case; this must not quietly start supporting it.
_lazycrop(r, ::Nothing, ::Integer) = r

function _lazycrop(r, cut, scale::Integer)
    _crsunit(Rasters.crs(r)) === ° || return r
    rows = _blockrange(dims(r, Y), ustrip(°, cut.Y[1]),
                       ustrip(°, cut.Y[2]), scale)
    cols = _blockrange(dims(r, X), ustrip(°, cut.X[1]),
                       ustrip(°, cut.X[2]), scale)
    (isnothing(rows) || isnothing(cols)) && return r
    return r[Y(rows), X(cols)]
end

# Read a raster file, optionally coarsening it by an integer `scale` factor, aggregating each block
# with `fn`. The opens run under a `NullLogger` to drop
# Rasters' benign metadata warnings (e.g. the `-99999` nodata CHELSA declares for its `UInt` bands, which
# the eltype can't hold). `_mask_int_fills` then removes the raw integer fill sentinels that declared nodata
# misses. Aggregating a *lazy* (disk-backed) raster reads it window-by-window - slow and hugely allocating -
# so when it comfortably fits in RAM we read it whole first (~6× faster and far fewer allocations); larger
# files fall back to the lazy aggregate to stay within memory. `_rastertodimarray` then materialises the
# (small) aggregated result.
function _readraster(f::AbstractString; scale::Integer = 1, fn = mean,
                     cut = nothing)
    r, raw = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        return Raster(f, lazy = true), Raster(f, lazy = true, scaled = false)
    end
    r = _mask_int_fills(r, raw)
    r = _lazycrop(r, cut, scale)
    scale > 1 || return r
    fits = _readbytes(r) < _READ_WHOLE_FRACTION * Sys.total_memory()
    return Rasters.aggregate(fn, fits ? read(r, checkmem = false) : r, scale)
end

# --- per-source read traits -------------------------------------------------
# After the Rasters migration the RasterDataSources-backed readers differ in only a few small,
# type-keyed ways (plus `_stackaxis`, above); one generic `read` (below) drives them all by
# dispatching these traits on the source type `T`. Each has a safe default and a per-source
# override where it actually differs.

# Physical unit attached to the data by the reader: none - `read` returns the raster's values as bare
# magnitudes in their actual physical unit (e.g. temperature in °C), and the unit is supplied from the
# shipped layer table (`layerunit`) when a layer is built. The directory readers below
# (`_readmonthlydir`) take a bare directory with no `RasterDataSource` attached, but their `var_name`
# uses the same short codes as the shipped tables (`tavg`, `wind`, `prec`, ...), so they self-attach a
# unit by deferring to `layerunit` too - `read(CHELSA{Climate}, dir, var_name)` against
# `CHELSA{Climate}`, `read(CRUTS, dir, var_name)` against `WorldClim{Climate}` (CRU TS has no table
# of its own; its variable names are WorldClim's).
_layerunit(::Type, files) = NoUnits

# Default read-time block-aggregation factor + reducer (land cover is coarsened 10× by default).
# Land cover's 10× default is a fact about `EarthEnv` and lives in the extension, with the note
# below on why the reducer is a plain `mean`.
_defaultscale(::Type) = 1

_defaultfn(::Type) = mean

# **Land cover aggregates with a plain mean, and must not round.** Rounding each of the twelve
# per-class bands back to an integer independently means a stack that summed to 100 before
# aggregation need not afterwards: measured over Scotland at the default 10x, the RMS departure from
# 100 is 0.466 with a worst case of 3.0 percentage points, against 0.027 and 0.18 for a plain mean.
#
# The metric that flatters rounding - more cells summing to exactly 100 - is an artefact of
# **quantisation** rather than accuracy: rounding snaps most cells onto 100 and throws the rest
# further off. It matters more than it looks, because the available-ground combine sums **eight** of
# these bands, so the per-band errors accumulate. And a mean of percentages is legitimately
# fractional, all the more so now that land cover canonicalises to a fraction of 1.

# Default keywords forwarded to `getraster` (WorldClim monthly climate needs a month range).
_getrasterkw(::Type) = (;)

# Normalise `getraster`'s return (a single path, a per-time vector, or a keyed collection of
# per-layer paths) to a plain `Vector{String}`. Does *not* handle the multi-layer *and*
# multi-time case (`Vector{<:NamedTuple}`, one NamedTuple of per-layer paths per time step) -
# that shape needs per-layer grouping first, see `_readmultilayer`/`Base.read`'s dispatch on it.
_filelist(x::AbstractString) = [String(x)]

_filelist(x) = String[String(f) for f in values(x)]

# The aggregated per-file layer is deterministic in (file, scale, fn, unit); memoise the resulting
# `DimArray` to disk (JLD2) in EcoSISTEM's scratch subdir so later reads skip the (slow) `aggregate`.
# Only caches when `scale > 1` (scale-1 reads are already fast) and `fn` is cacheable (see
# `_aggcachepath`). Caching the `DimArray` - not the `Raster` - preserves the exact lat/long axes,
# which a GeoTIFF round-trip perturbs.
#
# A windowed read is never cached. The key is (file, scale, fn, unit) - deliberately a `stat`, not
# a read - and has no room for the window; caching under it would serve one window's data for
# another's request. Folding the window in was the alternative and is not worth it: the cache exists
# to skip a slow whole-file `aggregate`, which a window does not do anyway, and a per-window key
# would fill the cache with near-duplicates.
function _cachedlayer(f, scale, fn, u; cut = nothing)
    (scale > 1 && isnothing(cut)) ||
        return _rastertodimarray(_readraster(f, scale = scale, fn = fn,
                                             cut = cut), unit = u)
    path = _aggcachepath(f, scale, fn, u)
    if !isnothing(path) && isfile(path)
        return jldopen(path, "r") do io
            return io["layer"]
        end
    end
    layer = _rastertodimarray(_readraster(f, scale = scale, fn = fn), unit = u)
    if !isnothing(path)
        mkpath(dirname(path))
        jldsave(path, layer = layer)
    end
    return layer
end

# The first file `getraster` resolves a source to, whatever shape it returned - a single path, a
# per-time vector, or a keyed collection of named layers. Mirrors `_readraw`'s dispatch on the same
# three shapes; every file of one source shares a CRS, so the first is representative.
_firstfile(raw::Vector{<:NamedTuple}) = String(first(values(first(raw))))

_firstfile(raw::NamedTuple) = String(first(values(raw)))

_firstfile(raw) = first(_filelist(raw))

# Read a resolved set of raster file paths into a `ClimateRaster` of source `T`. Values are returned in
# their actual physical unit as bare magnitudes (`_layerunit` is `NoUnits` for every source); the stacked
# axis (bands or a monthly series) comes from `_stackaxis`. Shared by `read` and the deprecated `readworldclim`.
function _readsource(T::Type, files::Vector{String};
                     cut = nothing, scale = _defaultscale(T),
                     fn = _defaultfn(T), slices = nothing)
    u = _layerunit(T, files)
    aas = map(f -> _cachedlayer(f, scale, fn, u, cut = cut), files)
    world = _stacklayers(aas,
                         _stackcoords(_stackaxis(T), length(aas), slices))
    return ClimateRaster(T, _applycut(world, cut))
end

# Read a `getraster` result shaped `Vector{<:NamedTuple}` - multiple named layers, each with its
# own per-time file list (the shape `getraster` returns for e.g. `layers = [:tmin, :tavg]` with
# `month = 1:12`) - into one multi-layer `ClimateRaster`. Each layer is read via the ordinary
# per-layer stacking (so it keeps its own `Ti`/time axis as usual), then the per-layer arrays are
# combined along a new `Dim{:layer}` axis (the generalised `_stacklayers`) labelled with the
# layers' own names - `getraster` already returns genuinely canonical `Symbol` keys (the shipped
# layer table's own codes), no further name resolution needed here.
function _readmultilayer(T::Type,
                         raw::Vector{<:NamedTuple};
                         cut = nothing, scale = _defaultscale(T),
                         fn = _defaultfn(T), slices = nothing)
    layernames = collect(keys(first(raw)))
    perlayer = map(layernames) do name
        files = String[String(nt[name]) for nt in raw]
        u = _layerunit(T, files)
        aas = map(f -> _cachedlayer(f, scale, fn, u, cut = cut), files)
        return _stacklayers(aas,
                            _stackcoords(_stackaxis(T), length(aas), slices))
    end
    world = _stacklayers(perlayer, Dim{:layer}(layernames))
    return ClimateRaster(T, _applycut(world, cut))
end

# Dispatch on the shape `getraster` actually returned - a `NamedTuple` (multiple named layers,
# one time step) or `Vector{<:NamedTuple}` (multiple named layers, each with its own per-time
# file list) both need canonical-name per-layer grouping before they can be stacked, so both
# route through `_readmultilayer` (wrapping a bare `NamedTuple` as a length-1 "time" vector);
# every other shape (a single path, or a per-time vector for one already-named layer) has no
# layer names to preserve and flattens directly via `_filelist`.
function _readraw(T, raw::Vector{<:NamedTuple}; cut, scale, fn,
                  slices = nothing)
    return _readmultilayer(T, raw, cut = cut,
                           scale = scale, fn = fn, slices = slices)
end

# A bare `NamedTuple` is *one* time step wrapped as a length-1 vector, so it has no slice identity
# to carry however the read was requested - `slices` is deliberately not forwarded.
function _readraw(T, raw::NamedTuple; cut, scale, fn, slices = nothing)
    return _readmultilayer(T, [raw], cut = cut, scale = scale,
                           fn = fn)
end

function _readraw(T, raw; cut, scale, fn, slices = nothing)
    return _readsource(T, _filelist(raw), cut = cut, scale = scale,
                       fn = fn, slices = slices)
end

# Divide out any `PublishedScaleFactor` - the layers a provider publishes at a known multiple of its
# own documented unit (see `LayerCatalogue.jl`). **The check runs on the RAW values and the division
# after**, and that order is the whole of it: correcting first would leave `_checkupstreamscale`
# looking at values it had just fixed, so it could never fire.
#
# Here rather than deeper in the read path because this is the only point that still knows which
# layer *codes* were asked for - below it, `getraster` has turned them into file paths.
function _rescalepublished(T::Type, layers,
                           cr::ClimateRaster)
    a = cr.array
    # A multi-layer read stacks along `Dim{:layer}` and each member may differ, so each is checked
    # and scaled by its own factor; a single-layer read has no such dim and takes the code asked for.
    if hasdim(a, Dim{:layer})
        codes = collect(DimensionalData.lookup(a, Dim{:layer}))
        factors = map(c -> something(_publishedscale(T, c), 1 // 1), codes)
        all(isone, factors) && return cr
        for (c, f) in zip(codes, factors)
            isone(f) || _checkupstreamscale(T, c, view(a, Dim{:layer}(At(c))))
        end
        shape = ntuple(i -> dimnum(a, Dim{:layer}) == i ? length(factors) : 1,
                       ndims(a))
        return ClimateRaster(T, a ./ reshape(float.(factors), shape), cr.code)
    end
    codes = layers isa AbstractVector ? layers : [layers]
    length(codes) == 1 || return cr        # nothing to attribute a factor to
    factor = _publishedscale(T, only(codes))
    isnothing(factor) && return cr
    _checkupstreamscale(T, only(codes), a)
    return ClimateRaster(T, a ./ float(factor), cr.code)
end

# Assemble a monthly time series from every `.tif` in `dir`, tagged with the already-resolved unit
# `u`. Shared by the directory-based `CRUTS`/`CHELSA{Climate}` readers, which differ only in their
# wrapper type and which shipped layer table `u` (from `layerunit`) comes from.
function _readmonthlydir(dir::String, u; scale = 1, fn = mean, cut = nothing)
    files = joinpath.(dir, _searchdir(dir, ".tif"))
    aas = map(f -> _rastertodimarray(_readraster(f, scale = scale, fn = fn),
                                     unit = u), files)
    world = _stacklayers(aas, _mkstackaxis(Ti, length(aas)))
    return _applycut(world, cut)
end
