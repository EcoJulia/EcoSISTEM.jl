# SPDX-License-Identifier: LGPL-3.0-or-later
#
# GETTING RASTER DATA ONTO A GRID. Everything between a file on disk and a value in the right cell:
#
# reading        `_read`, `_asraster`, `_attachunit`, `_applyperiod` - a spec to a `ClimateRaster`
# reprojection   `_reproject`, `_resample2d`, `_targetcrs`, `_cropto`, `_cellintervals` - the bulk
# masks          `_shapegeoms`, `_circle`, `_coverage`, `_rastermask` - which cells are active
#
# `materialise.jl` sits directly on top of this and adds the cache and the layer wrapper;
# `StudyArea.jl` uses it to decide a grid before anything is built. Neither is here: this file knows
# about rasters and grids, not about layers, areas or habitats.

using Unitful

using Unitful.DefaultSymbols

using EcoSISTEM.Units

using Distributions

using Random

using Dates: Dates

using DimensionalData

using Rasters

import Rasters: Projected

# **Nothing here is imported from `EcoSISTEM.ClimatePref`, and that is the point.** `_crsunit`,
# `_isangle` and `_stacklayers` are raster geometry rather than climate data, so they live in the
# parent module. The dependency runs one way, and the submodule is the dependent.
# Load-bearing with no referenced symbol: it registers the GDAL backend that `Rasters` needs for
# reprojection. Do not remove in a dead-code pass.
import ArchGDAL

# Every projected EPSG CRS that PROJ's own database knows a usable extent for, read once and kept.
# Deprecated entries and those without a declared bounding box are dropped, as is anything whose box
# crosses the antimeridian (`west > east`) - the extents this engine builds never wrap, so such a CRS
# could never be compared against one sensibly.
#
# This is `proj.db`, reached through GDAL's OSR bindings, so it needs no dependency beyond the
# ArchGDAL already in use - deliberately not a hand-maintained table, which would go stale and could
# never cover the ~5000 CRSs available here.
const _CRS_CANDIDATES = Ref{Union{Nothing, Vector{NamedTuple}}}(nothing)

# --- Active-area masks -----------------------------------------------------
# Data-driven active-area masks are composed with `ConstructedSpec` from a data source plus a
# combine rule; `CircleMaskSpec`/`ShapeSpec` are the synthetic/vector mask specs. The two public
# rules below are reusable building blocks for writing your own combine (`_circle`/`_shapegeoms`+
# `_shape` remain the private geometry helpers for the synthetic/vector masks).

"""
    hasdata(layer)

A [`ConstructedSpec`](@ref) combine returning a `Bool` mask of the cells of `layer` that hold data
(are not missing/`NaN`) - the canonical combine-rule example. Pass it a data source to mask that
source's coverage: `ConstructedSpec(hasdata, WorldClim{BioClim}, 1)`.

This is about **data coverage**, not about whether a cell is active: it takes a raw
[`ClimateRaster`](@ref), so it runs *before* any active mask exists - it is one of the rules that
**produces** one.

Two limits worth knowing. It tests `NaN` only, so a source's own nodata sentinel counts as data
unless reading has already converted it (which it does). And on a multi-band raster it reports the
**first band** alone, so a twelve-class land-cover stack is masked by class 1's coverage rather than
the stack's.
"""
function hasdata(layer::ClimateRaster{S}) where {S}
    A = Array(layer.array)
    # A 2-D layer masks by plain broadcast, so it stays a raster and keeps its own code - a mask
    # of one layer is still identifiably that layer's.
    ndims(A) == 2 && return .!isnan.(layer)
    # A multi-band read collapses to its **first** band, which is a different layer from the stack
    # it came from, so the code is deliberately not carried across.
    return ClimateRaster(S,
                         DimArray(Matrix{Bool}(.!isnan.(A[:, :, 1])),
                                  (dims(layer.array, Y), dims(layer.array, X))))
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Public builders
# ---------------------------------------------------------------------------

# `_canonical(value, axis)` (axis-driven, in `src/Layer.jl`) canonicalises a regime value to its axis's
# unit; `_reaxis(layer, axis)` (also `src/Layer.jl`) tags the layer with that axis.

# **There is deliberately only one `materialise`, and it takes a `StudyArea`.** A second family
# taking a bare shape - `materialise(spec, dim, size)` - would be a second implementation of what the
# `StudyArea` form already does through `_syntheticsupplyfield`, and the two would have to agree by
# inspection. They give identical results for `UniformSpec`, `GradientSpec` and `PeakedSpec`,
# including the affine °C to K case, where canonicalising the endpoints
# and canonicalising the whole field agree because the offset cancels out of the interpolation), and
# `NicheSpec` differs only because `_randomniches` is stochastic. Its two private generators
# (`_gradient_regime`, `_peaked_regime`) went with it, having no other caller.

# ---------------------------------------------------------------------------
# Data-driven environments: reconcile one or more rasters onto a common grid
# and build regime + supply + active mask from them.
# ---------------------------------------------------------------------------

# Convert a step in degrees to a physical distance (≈111.32 km per degree, as the existing `*AE`
# constructors assume). Only the degrees->length direction survives: it is still needed for the
# *geographic* branch of `_circle`. The old length->degrees `_deg` is gone - a physical cell side is
# now only meaningful on a projected target, where it needs no conversion at all.
# Renamed from `_side` 2026-08-24: that name is ALSO a three-argument helper in
# `collections.jl` meaning *a side of a pairing check*. One generic function, two unrelated
# meanings, colliding only by arity - the same confusion `_inrange`/`_indexrange` had.
_degreelength(step_deg) = abs(uconvert(NoUnits, step_deg / °)) * 111.32km

# --- Combination-time reprojection ------------------------------------------
# `Rasters.resample` (GDAL `warp`-backed) needs (confirmed empirically against real British
# National Grid CRS data, not assumed): (a) an actual `Rasters.Raster`, not a bare
# `DimensionalData.AbstractDimArray`; (b) unitless coordinate values on `Y`/`X` - a genuine
# `Unitful.Quantity`-valued dim breaks GDAL's own geotransform construction; (c) unitless data
# values, for the same reason (GDAL needs a concrete numeric `GDALDataType`). `_reproject` strips
# both, resamples, and reattaches (data unit as given; axis unit from the target's own CRS via
# `_crsunit`, defined in this file rather than re-derived).

# `d` (a `Y`/`X` dim) with its coordinate values stripped of their Unitful unit - CRS, locus and
# order preserved exactly; the span is always rebuilt as `Regular` from the actual (stripped) step
# rather than copied, since a dim that has been index-cropped comes
# back with an `Irregular` span (computed min/max bounds, not a step) even though the underlying
# grid is still perfectly regular - `Rasters.resample` requires a genuinely `Regular` span on its
# `to =` target. Every raster this codebase handles is a regular grid, so recomputing is safe.
#
# **This is the only place a coordinate loses its unit**, and it is the only
# place that may: every grid this package builds or reads carries one, and GDAL is the single
# consumer that refuses them (`Rasters.resample` raises `DimensionError` on a unitful `to =`;
# everything else - `crs`, `bounds`, `aggregate`, `extent`, `selectindices` - is indifferent).
# **The unit is named rather than assumed.** `ustrip(u, x)` converts *then* strips, so a `km` dim
# and an `m` dim against the same CRS come out on the same scale, where a bare `ustrip.(x)` would
# take the numbers as written and silently mix the two.
# The non-`Quantity` branch is not dead: a caller may still hand over an already-bare dim (a hand
# built test fixture, or a second pass over a stripped raster), and `ustrip(m, 5.0)` would throw.
function _striplength(d)
    lk = DimensionalData.lookup(d)
    raw = parent(lk)
    vals = eltype(raw) <: Unitful.Quantity ?
           ustrip.(_crsunit(Rasters.crs(d)), raw) : raw
    step = length(vals) > 1 ? vals[2] - vals[1] : oneunit(eltype(vals))
    return DimensionalData.rebuild(d,
                                   Projected(vals,
                                             sampling = DimensionalData.sampling(d),
                                             crs = Rasters.crs(d),
                                             order = DimensionalData.order(lk),
                                             span = DimensionalData.Lookups.Regular(step)))
end

# `r` with unitless `Y`/`X` dims; any other dimension is untouched.
function _unitlessyx(r)
    yd, xd = _striplength(dims(r, Y)), _striplength(dims(r, X))
    return DimensionalData.rebuild(r,
                                   dims = map(d -> d isa Y ? yd :
                                                   d isa X ? xd : d, dims(r)))
end

# Resample one 2-D `(Y, X)` slice `r2d` onto `target` (already unitless `Y`/`X`) with GDAL method
# `method`. `resample`'s own missing-data convention is `missing`; converted back to `NaN` here to
# match this codebase's established sentinel (`_rastertodimarray`'s
# `Rasters.replace_missing(ras, NaN)`), rather than introducing a second missing-data convention.
function _resample2d(r2d, target, method)
    out = Rasters.resample(r2d, to = target, method = method)
    return Rasters.rebuild(out, data = Float64.(coalesce.(parent(out), NaN)))
end

# Resample `r` - 2-D, or N-D with exactly one extra non-spatial dimension (e.g. `Ti`/
# `Dim{:layer}`) - onto `target` (a unitless `(Y, X)` `Rasters.Raster` in the destination CRS/grid),
# reprojecting if `r`'s CRS differs. `method` is a GDAL resampling method `Symbol` (`:near`/
# `:mode`/`:bilinear`/`:cubic`/..., chosen by the caller - see `_resamplemethod`, which will choose it
# from the shipped `ValueType` column once the `Band` redesign gives a materialised layer its code). A 3-D+ array is resampled one 2-D slice at a time - `Rasters.resample` errors directly
# on 3-D input (`GDALError: Too few arguments for '-te'`, confirmed) - and restacked along the
# original axis values, so a real `Ti{DateTime}`/`Dim{:layer}` axis survives unchanged. Returns a
# plain `DimArray` (not a `Rasters.Raster`), matching `ClimateRaster.array`'s existing convention.
function _reproject(r::Rasters.AbstractRaster, target; method = :bilinear)
    dataunit = eltype(r) <: Unitful.Quantity ? unit(eltype(r)) : NoUnits
    stripped = dataunit === NoUnits ? r :
               Rasters.rebuild(r, data = ustrip.(parent(r)))
    u = _unitlessyx(stripped)
    # **The target is stripped too, and this is the whole of GDAL's boundary.** Since
    # the study area's grid carries its coordinate unit like everything else, and
    # `Rasters.resample` is the one consumer that refuses one - so it is undressed here, once, rather
    # than the grid being kept bare everywhere on its account. Done **outside** the slice loop: an
    # N-D raster resamples one 2-D slice at a time, and stripping per slice would repeat identical
    # work for every month of a monthly stack.
    to = _unitlessyx(target)
    extras = otherdims(u, (Y, X))
    resampled = if isempty(extras)
        _resample2d(u, to, method)
    else
        length(extras) == 1 ||
            error("_reproject only supports one extra (non-spatial) dimension; got $(length(extras))")
        ed = extras[1]
        bt = DimensionalData.basetypeof(ed)
        slices = [_resample2d(u[bt(i)], to, method) for i in eachindex(ed)]
        cat(slices..., dims = bt(DimensionalData.lookup(ed)))
    end
    yd, xd = _targetyx(target)
    newdims = map(d -> d isa Y ? yd : d isa X ? xd : d, dims(resampled))
    data = dataunit === NoUnits ? parent(resampled) :
           parent(resampled) .* dataunit
    return DimArray(data, newdims)
end

# GDAL resampling method for putting `raster`'s data on a common grid: a class-code layer must never
# be interpolated *between* classes, so it takes the nearest class by frequency (`:mode`); anything
# numeric may be interpolated (`:bilinear`). Shares `iscategorical` with the regime-type choice
# rather than classifying twice - it is one question.
#
# **Takes the answer, not the raster.** Whether a layer holds class codes is decided by its axis,
# which the raster does not carry - the caller does, from the spec. A **mask** passes `true`
# unconditionally: it is `Bool`-valued and has no niche axis at all, and interpolating between `true`
# and `false` is meaningless in exactly the way interpolating between class codes is.
_resamplemethod(categorical::Bool) = categorical ? :mode : :bilinear

# The `(Y, X)` dims any layer sampled onto `target` carries: the target's own coordinates, re-united
# in the target CRS's coordinate unit. Shared by both sampling routes so that a cropped layer and a
# warped one are indistinguishable afterwards.
#
# The span must be restated, not left to inference. `Projected` over a plain `Vector` infers
# `Irregular((nothing, nothing))` - a span with no bounds - which is the same trap `_spatialdim`
# guards against on the read path. It bites differently here: a cropped layer keeps its source's
# `Regular` span while a warped one would come back `Irregular`, and two layers of one collection
# sampled by different routes then compare as being on *different grids*, which `_yx` rejects
# outright.
#
# **Simply the target's own dims.** A grid stored bare would need its coordinate unit put back here,
# and again in `_unitedyx` for the synthetic path - and two copies of one re-attachment, disagreeing
# about whether to run at all, is a drift that costs cells silently. The grid carries its unit, so
# there is nothing to put back.
# The function stays because the *span* restatement above is still doing work: the resampled
# output's dims come back from GDAL, so they are replaced wholesale by the target's, which are known
# to be `Regular`.
_targetyx(target) = dims(target, (Y, X))

# The contiguous run of `sourcevals` equal to `targetvals`, or `nothing` if the target's coordinates
# are not simply a stretch of the source's. Exact equality is the right test and not a fragile one:
# both sides come off the same corrected lattice, and the alignment layer's target is built by
# cropping that layer's *own* grid, so the values are the identical `Float64`s. Anything else falls
# back to resampling, which is what happened before.
function _subgridrange(sourcevals, targetvals)
    isempty(targetvals) && return nothing
    i = findfirst(==(first(targetvals)), sourcevals)
    isnothing(i) && return nothing
    j = i + length(targetvals) - 1
    j <= length(sourcevals) || return nothing
    return view(sourcevals, i:j) == targetvals ? (i:j) : nothing
end

# `raster` cut to `target` by plain indexing, when `target` is literally a sub-grid of it - same CRS,
# same cells, cell for cell - or `nothing` when it is not and a real resample is needed.
#
# This exists because a GDAL warp is **not** a no-op even onto an identical grid, and the
# difference is not the interpolation itself. Values do pass through untouched (bilinear's weights
# collapse to (1,0,0,0) on coincident cells; measured max difference 2e-9). What changes is *nodata*:
# the warp poisons any output cell whose input stencil touches a `NaN`, so the coastline erodes - and
# precisely which cells that catches depends on GDAL's internal source windowing, which the extent of
# the read changes. Reading a global layer and reading a window of the same layer therefore produced
# study areas differing in 138 coastal cells (650 vs 614 active), every one of them adjacent to
# nodata. Cropping removes the warp entirely, so the answer is invariant *by construction* rather
# than by GDAL behaving consistently.
#
# It also makes the report honest. `_resamplecost` already classifies this case as `:exact` and
# `StudyAreaReport` prints "kept exactly" - which was false while a warp ran over the layer. The
# geometric verdict and what is actually done to the data now agree.
function _cropto(raster::ClimateRaster, target)
    A = raster.array
    yd, xd = dims(A, Y), dims(A, X)
    tcrs = Rasters.crs(target)
    _samecrs(Rasters.crs(yd), tcrs) || return nothing
    # **Compared with units on, on both sides.** Source and target both carry their coordinate unit
    # and `_subgridrange` only ever tests `==`, which Unitful evaluates across
    # scales, so a `km`-labelled source matches an `m`-labelled target of the same ground.
    # **Stripping only the source here was a real bug for exactly one commit**: the target became
    # unitful while this still bared its own side, so every `==` failed, `_cropto` returned `nothing`,
    # and an *exact crop* silently became a **resample** - which moved the grid by a row.
    # Caught by `test_StudyArea.jl:90` (*"windowing the reads does not change the answer"*),
    # `(37, 47)` against `(38, 47)`. A reminder that the crop/resample fork is invisible in the
    # result's *values* and only shows in its shape.
    rows = _subgridrange(parent(DimensionalData.lookup(yd)),
                         parent(DimensionalData.lookup(target, Y)))
    cols = _subgridrange(parent(DimensionalData.lookup(xd)),
                         parent(DimensionalData.lookup(target, X)))
    (isnothing(rows) || isnothing(cols)) && return nothing
    cropped = A[Y(rows), X(cols)]
    yd, xd = _targetyx(target)
    newdims = map(d -> d isa Y ? yd : d isa X ? xd : d, dims(cropped))
    return DimArray(parent(cropped), newdims)
end

# Reproject `raster` onto `target` (real interpolation/reprojection via `_reproject`, replacing a
# hand-rolled nearest-neighbour lookup), erroring clearly if `raster` turns out to have no real
# coverage of `target` at all (rather than silently building an all-inactive environment).
#
# Deliberately says nothing about a change of resolution: what the target grid costs each layer -
# copied exactly, aggregated by a whole factor, or genuinely resampled, and why - is classified per
# layer by `_analyse` and reported by the `StudyArea` that decided the grid. The blanket warning this
# would fire on exactly the case that classification calls an *exact* aggregation, so the
# two contradicted each other.
function _sampledata(raster::ClimateRaster, target; name = "raster",
                     categorical::Bool)
    # `_reproject` already returns a real `(Y, X[, extra])` `DimArray` - kept as-is (not stripped
    # to a bare `Array`) so the regime/supply built from it carries real CRS provenance, per
    # the `(y, x)` order used throughout.
    # Not `something(...)`: that would evaluate the resample eagerly, which is the whole cost
    # this is here to avoid.
    cropped = _cropto(raster, target)
    out = isnothing(cropped) ?
          _reproject(Rasters.Raster(raster.array), target,
                     method = _resamplemethod(categorical)) : cropped
    any(!isnan, out) ||
        error("$name does not overlap the study area's grid at all - check the layer's real extent " *
              "against the `StudyArea`, and the area's `within`/`crs` if given.")
    return out
end

# `target`'s own Y/X coordinate values, in its axis unit (° for a geographic target, a real length
# for a projected one).
#
# **Simply read**, because the grid carries its own unit. Multiplying by `_crsunit(crs)` would be
# needed only for a grid stored bare, and would cost twice over: a **synthetic** grid has no CRS but
# does have units, and `_crsunit(nothing)` answers `°`, so a `km` grid would come back as `° km` - a
# dimensionally absurd value that is merely *accidentally* harmless where it is compared against `°`.
# It would also leave a caller unable to tell a length grid from an angular one without knowing the
# rule.
#
# **Handing out `parent(lookup(dim))` instead - raw values with no stated meaning - leaves every
# caller to decide what they are, and they will not agree.** One reads them as cell **centres**,
# another as centres it then reconstructs edges from, a third as a span; so a change of *locus*, a
# labelling choice that moves no cell at all, propagates into all of them at once.
#
# **`DimensionalData.Lookups.intervalbounds` is locus-blind by construction**: measured, the
# same three cells described as `Start` and as `Center` give **identical** per-cell `(lo, hi)`. So a
# caller asks for what it means - the interval, the midpoint, or the span - and the convention stops
# being load-bearing anywhere outside the one function that builds the grid.
# A hand-rolled `± step/2` would have been the same mistake in a new place: it is itself a locus
# assumption.
#
# Both come in a two-argument form (one axis) and a no-argument-dimension form returning a named
# `(lat, long)` - the pairing `[UNIT-DUP]` wanted, since nearly every caller needs both.
function _cellintervals(target, D)
    return DimensionalData.Lookups.intervalbounds(DimensionalData.lookup(target,
                                                                         D))
end

function _cellintervals(target)
    return (lat = _cellintervals(target, Y), long = _cellintervals(target, X))
end

# The representative point of each cell - its midpoint, whatever the lookup labels its cells by.
function _cellcentres(target, D)
    return [(lo + hi) / 2 for (lo, hi) in _cellintervals(target, D)]
end

function _cellcentres(target)
    return (lat = _cellcentres(target, Y), long = _cellcentres(target, X))
end

# A CRS as an `ArchGDAL` spatial reference in *traditional* (long/easting-first) axis order. The
# `order = :trad` matters: some authority definitions - notably EPSG:4326 - declare (lat, long) order,
# which would silently swap the coordinates against the `createpoint(x, y)` convention used throughout
# this file.
_gdalcrs(crs) = ArchGDAL.importCRS(crs, order = :trad)

# The four bounds of an extent, in the order the old `(ylo, yhi, xlo, xhi)` convention used. An
# escape hatch for the handful of places that genuinely need the numbers -- a message naming all
# four, or a component-wise comparison -- not a licence to go back to passing them around.
_extentvalues(e::Extents.Extent) = (e.Y[1], e.Y[2], e.X[1], e.X[2])

# Build one from the same four, for the boundary where a foreign API hands them over separately.
_extentof(ylo, yhi, xlo, xhi) = Extents.Extent(Y = (ylo, yhi), X = (xlo, xhi))

# A geographic place expressed in `crs`'s own coordinates. A `LatLong` in, a `SpatialLocation` out --
# and the types carry the frame, which is the whole reason they exist: the input is degrees, the
# result is whatever `crs` measures in, metres for a projected one.
#
# A geographic `crs` needs no transform at all, its coordinates *being* lat/long, so the place comes
# back unchanged and still typed as a `LatLong`.
function _pointin(crs, place::LatLong)
    _isprojectedcrs(crs) || return place
    point = ArchGDAL.createpoint(ustrip(°, getlong(place)),
                                 ustrip(°, getlat(place)))
    ArchGDAL.createcoordtrans(_gdalcrs(Rasters.EPSG(4326)), _gdalcrs(crs)) do ct
        return ArchGDAL.transform!(point, ct)
    end
    u = _crsunit(crs)
    return SpatialLocation(ArchGDAL.gety(point, 0) * u,
                           ArchGDAL.getx(point, 0) * u)
end

# Each cell's angular `(y, x)` extent on a **projected** grid: the inverse of the transform every
# other function here applies, run once per cell.
#
# A projected cell has a real angular size -- it covers so many degrees of latitude and longitude --
# but getting it means transforming its coordinates back through the CRS, which is what this does.
# The extent VARIES across the grid, and on an equal-area projection it is precisely the quantity
# that varies while the area does not, so a single answer would be wrong wherever it was not taken.
#
# The transform is created ONCE and reused across every cell. `_bboxin` would have been the obvious
# reuse, but it builds a fresh `createcoordtrans` per call, which on a million-cell grid is a million
# GDAL context set-ups; the loop below pays that cost once.
#
# Each cell is transformed by its two opposite corners rather than all four. That is exact for the
# axis-aligned case and slightly under-states a rotated one, which is why the caller documents the
# result as the cell's extent along each axis rather than its bounding box.
function _angularextents(yx, crs)
    ys, xs = parent(DimensionalData.lookup(yx[1])),
             parent(DimensionalData.lookup(yx[2]))
    dy, dx = _axisstep(yx[1]), _axisstep(yx[2])
    ny, nx = length(ys), length(xs)
    ey, ex = Matrix{typeof(1.0°)}(undef, ny, nx),
             Matrix{typeof(1.0°)}(undef, ny, nx)
    ArchGDAL.createcoordtrans(_gdalcrs(crs),
                              _gdalcrs(Rasters.EPSG(4326))) do ct
        for i in 1:ny, j in 1:nx
            lo = ArchGDAL.createpoint(ustrip(xs[j]), ustrip(ys[i]))
            hi = ArchGDAL.createpoint(ustrip(xs[j] + dx), ustrip(ys[i] + dy))
            ArchGDAL.transform!(lo, ct)
            ArchGDAL.transform!(hi, ct)
            ey[i, j] = abs(ArchGDAL.gety(hi, 0) - ArchGDAL.gety(lo, 0)) * °
            ex[i, j] = abs(ArchGDAL.getx(hi, 0) - ArchGDAL.getx(lo, 0)) * °
        end
    end
    return (y = ey, x = ex)
end

# `src`'s extent re-expressed in `dst`, as the axis-aligned
# envelope of its four transformed corners (a rectangle in one CRS is not one in another, so the
# envelope is the conservative box that still encloses it). `nothing` for either CRS, or two equal
# CRSs, is a no-op.
function _bboxin(src, dst, e::Extents.Extent)
    (isnothing(src) || isnothing(dst) || _samecrs(src, dst)) && return e
    u = _crsunit(dst)
    ys, xs = Float64[], Float64[]
    ArchGDAL.createcoordtrans(_gdalcrs(src), _gdalcrs(dst)) do ct
        for y in e.Y, x in e.X
            point = ArchGDAL.createpoint(ustrip(x), ustrip(y))
            ArchGDAL.transform!(point, ct)
            push!(ys, ArchGDAL.gety(point, 0))
            push!(xs, ArchGDAL.getx(point, 0))
        end
    end
    return _extentof(minimum(ys) * u, maximum(ys) * u,
                     minimum(xs) * u, maximum(xs) * u)
end

# The bounds of an already-gridded array's own `Y`/`X`, expressed in `tcrs`. `nothing` when the array
# carries no CRS (a synthetic/`NoLookup` grid has no real-world position to offer).
# The ground a dimension actually covers: the outer *edges* of its first and last cells.
#
# Not the extrema of its coordinate values. A lookup value is a
# cell's reference corner - with the `Start` locus GDAL files declare, its lower-left one - so the
# largest value is where the last cell *begins*, and the layer reaches a further whole cell beyond it.
# WorldClim's Y values run -90°...89.8333°, but the layer covers -90°...90°. Measuring at the reference
# corner therefore threw away a full cell off the top and right of **every** layer, and since
# `Touches` (which `_applycut` selects with) works on intervals, a windowed read came back slightly
# inside what was asked for and then clipped the mask in `_planbounds`.
#
# Only a `Regular` span can be widened to cell edges, and asking anything else is not merely
# useless but unsafe. An `Irregular` span with no recorded ends - what `Rasters.aggregate` produces -
# answers `(nothing, nothing)`, which would propagate straight into `_bboxin`'s comparisons; and a
# lookup left at `AutoSpan`/`AutoSampling`, which a hand-built `Projected` mask can be, makes
# `Lookups.bounds` **throw** a `MethodError` rather than return anything at all. Both fall back to
# the coordinate values rather than the interval bounds, so the fallback is the
# previous behaviour rather than a new guess.
function _axisbounds(d)
    DimensionalData.Lookups.span(d) isa DimensionalData.Lookups.Regular ||
        return extrema(parent(DimensionalData.lookup(d)))
    b = DimensionalData.Lookups.bounds(d)
    any(isnothing, b) && return extrema(parent(DimensionalData.lookup(d)))
    return extrema(b)
end

# An array's extent expressed in `tcrs`, or `nothing` when it carries no CRS to convert from. Read
# from the dims rather than from any stored metadata, so it describes where the cells actually are.
function _dimsextent(A, tcrs)
    yd, xd = dims(A, Y), dims(A, X)
    acrs = Rasters.crs(yd)
    isnothing(acrs) && return nothing
    ylo, yhi, xlo, xhi = _extrema2(_axisbounds(yd), _axisbounds(xd))
    return _bboxin(acrs, tcrs, _extentof(ylo, yhi, xlo, xhi))
end

# Which cells of `target` lie **wholly inside** the ground `raster` covers - `simulate_safely`'s test.
#
# Purely geometric: a per-axis interval test and their outer product, so it needs no values and no
# sampling. That is deliberate - what it answers is *"does the data describe all of this cell"*, which
# is a question about the footprint, not about any particular pixel.
#
# **"Covered" means inside the footprint, NOT "every pixel has a value."** The second reading would
# drop every coastal cell, since a 10 km cell on a coastline is part sea and so partly `NaN`. Absence
# *within* the extent stays as it always was - a cell whose centre is `NaN` is already inactive.
#
# **Asked in the LAYER's own coordinates whenever the CRSs differ, never in the target's.** This
# is `_clipto`'s rule and for `_clipto`'s reason, and it is not a refinement:
# EarthEnv's global WGS84 footprint (-56°...90°, -180°...180°) re-expressed in British National Grid comes
# back as eastings **275 286 m...400 000 m** - a 125 km strip. Scotland's own grid runs 0...500 km, so
# most of it read as *outside the data* and `ScottishCultivatedLand.jl` lost **1859 of 3168** active
# cells. Re-expressing a *global* layer's four corners in a local projected CRS is meaningless; the
# small box into the big CRS is the well-conditioned direction, and answers `true` for the same case.
#
# A raster with no CRS - or a target with none - cannot be positioned at all, so it claims full
# coverage rather than inventing an answer: a synthetic grid has no footprint to fall outside of.
_fullycovered(r::ClimateRaster, target) = _fullycovered(r.array, target)

function _fullycovered(A, target)
    everywhere = trues(Base.size(target)[1:2])
    tcrs = Rasters.crs(target)
    yd, xd = dims(A, Y), dims(A, X)
    acrs = Rasters.crs(yd)
    (isnothing(tcrs) || isnothing(acrs)) && return everywhere
    native = _extentof(_extrema2(_axisbounds(yd), _axisbounds(xd))...)
    yb, xb = _cellintervals(target, Y), _cellintervals(target, X)
    _samecrs(acrs, tcrs) && return _cellsinbox(yb, xb, native)
    # One box into the layer's CRS first: when the layer simply contains the area - the common case,
    # and the one the target-CRS envelope got wrong - that settles every cell at the cost of four
    # point transforms, with no per-cell work at all.
    whole = _bboxin(tcrs, acrs,
                    _extentof(minimum(minimum, yb), maximum(maximum, yb),
                              minimum(minimum, xb), maximum(maximum, xb)))
    _boxwithin(whole, native) && return everywhere
    return _cellsinbox(_reprojectedcells(yb, xb, tcrs, acrs)..., native)
end

# The per-cell answer once both sides are in one CRS: a per-axis interval test, outer-producted.
function _cellsinbox(yb, xb, e::Extents.Extent)
    fy = [_intervalwithin(iv, e.Y[1], e.Y[2]) for iv in yb]
    fx = [_intervalwithin(iv, e.X[1], e.X[2]) for iv in xb]
    return fy .& permutedims(fx)
end

# Each of the target's cells re-expressed in `crs` as an axis-aligned envelope, given per axis.
#
# **An envelope, so the answer is conservative in the safe direction**: a reprojected cell is not a
# rectangle, and its envelope contains it - so "the envelope is inside the data" implies the cell is.
# The cost is that a cell may occasionally be dropped that was in fact covered, never kept when it was
# not, which is the way round `simulate_safely` wants to be wrong.
#
# The whole corner lattice goes through **one** coordinate transform: `(ny + 1) × (nx + 1)` points
# for `ny × nx` cells, rather than four per cell through `_bboxin` and a `createcoordtrans` with them.
function _reprojectedcells(yb, xb, tcrs, crs)
    yedges, xedges = _edgesequence(yb), _edgesequence(xb)
    u = _crsunit(crs)
    lat = Matrix{Float64}(undef, length(yedges), length(xedges))
    long = similar(lat)
    ArchGDAL.createcoordtrans(_gdalcrs(tcrs), _gdalcrs(crs)) do ct
        for (i, y) in enumerate(yedges), (j, x) in enumerate(xedges)
            point = ArchGDAL.createpoint(ustrip(x), ustrip(y))
            ArchGDAL.transform!(point, ct)
            lat[i, j] = ArchGDAL.gety(point, 0)
            long[i, j] = ArchGDAL.getx(point, 0)
        end
    end
    # A cell's envelope spans its own four corners, so both axes have to be reduced together - hence
    # the per-cell minima and maxima rather than one interval per axis.
    ny, nx = length(yb), length(xb)
    ys = [extrema((lat[i, j], lat[i + 1, j], lat[i, j + 1], lat[i + 1, j + 1]))
          for i in 1:ny, j in 1:nx]
    xs = [extrema((long[i, j], long[i + 1, j], long[i, j + 1],
                   long[i + 1, j + 1])) for i in 1:ny, j in 1:nx]
    # Reduced across the other axis, so the per-axis test still applies: a cell is covered iff its
    # widest reach on each axis is inside the data's.
    yiv = [(minimum(first, view(ys, i, :)) * u,
            maximum(last, view(ys, i, :)) * u) for i in 1:ny]
    xiv = [(minimum(first, view(xs, :, j)) * u,
            maximum(last, view(xs, :, j)) * u) for j in 1:nx]
    return yiv, xiv
end

# The `n + 1` edge coordinates of a run of `n` contiguous cell intervals, **in the lookup's own
# order**, so `(edges[i], edges[i + 1])` is cell `i`'s extent whichever way the axis runs. A
# descending axis (a north-up raster's latitude) walks its cells' upper edges downwards.
function _edgesequence(ivs)
    ascending = last(last(ivs)) >= first(first(ivs))
    edges = [ascending ? minimum(iv) : maximum(iv) for iv in ivs]
    push!(edges, ascending ? maximum(last(ivs)) : minimum(last(ivs)))
    return edges
end

# Whether one extent sits wholly inside another. **Strictly, with no tolerance**,
# unlike the per-cell test: this is the fast path's "then every cell is covered" shortcut, so slack
# here would claim full coverage on the strength of a rounding error rather than merely keep one
# borderline cell.
function _boxwithin(box::Extents.Extent, outer::Extents.Extent)
    return Extents.covers(_ordered(outer), _ordered(box))
end

# An extent with each axis low-to-high. The normalisation is NOT decorative: a `Y` lookup is
# routinely `ReverseOrdered` -- GDAL files run north to south -- so bounds taken from one arrive
# high-first, and `Extents`' predicates reasonably assume otherwise.
function _ordered(e::Extents.Extent)
    return _extentof(minmax(e.Y...)..., minmax(e.X...)...)
end

# Whether a cell's interval sits inside `[a, b]`, either axis given in either order.
#
# **The tolerance is not optional.** Where the data's span divides exactly by the cell size, the
# last cell's far edge is `lo + (n-1)*step + step` while the data's is `lo + n*step` - arithmetic that
# agrees mathematically and can differ by an ULP, which without slack silently drops a fully-covered
# row. A whisker of the cell's own width absorbs that and nothing else: at 10 km it is 0.15 mm.
function _intervalwithin((lo, hi), a, b)
    tol = sqrt(eps(Float64)) * abs(hi - lo)
    return min(lo, hi) >= min(a, b) - tol && max(lo, hi) <= max(a, b) + tol
end

# A circle's own bounds, `centre ± radius`, in `tcrs`. Only available when the centre is given *and*
# the target is projected: with the default centre the circle is defined *as* the grid's midpoint, so
# it cannot lead the extent (that would be circular), and on a geographic target converting a radius
# to degrees would need exactly the fixed-km-per-degree approximation this engine removed.
function _circleextent(cm::CircleMaskSpec, tcrs)
    (isnothing(cm.centre) || !_isprojectedcrs(tcrs)) && return nothing
    centre = _pointin(tcrs, cm.centre)
    r = uconvert(_crsunit(tcrs), cm.radius)
    return _extentof(centre.y - r, centre.y + r, centre.x - r,
                     centre.x + r)
end

# Resolve a mask as far as it can be *before* the target grid exists: the `payload` `_rastermask` will
# need once there is a grid, plus the `extent` the mask implies (`nothing` when it cannot state one, in
# which case it simply follows the data). Doing this once, here, is what stops a `ShapeSpec`'s file read
# or a `ConstructedSpec`'s data read happening twice - once to learn the extent and again to rasterise.
_preparemask(active::Nothing, tcrs) = (payload = nothing, extent = nothing)

# A plain `Matrix{Bool}` carries no coordinates at all, so it can only follow the data.
function _preparemask(active::AbstractMatrix{Bool}, tcrs)
    return (payload = active, extent = nothing)
end

function _preparemask(active::DimensionalData.AbstractDimArray, tcrs)
    return (payload = active, extent = _dimsextent(active, tcrs))
end

function _preparemask(active::CircleMaskSpec, tcrs)
    return (payload = active, extent = _circleextent(active, tcrs))
end

function _preparemask(active::ShapeSpec, tcrs)
    geoms, extent = _shapegeoms(active, tcrs)
    return (payload = geoms, extent = extent)
end

# A geographic `Extents.Extent` is *pure* extent: it says where the grid goes and nothing else, so once the grid has
# been cut to it every cell is inside and there is no payload left to rasterise. Always given in WGS84
# whatever the target CRS (a lat/long rectangle is not a rectangle in a projected one, so `_bboxin`
# takes the envelope of its transformed corners).
function _preparemask(active::Extents.Extent, tcrs)
    _checkgeographicextent(active)
    return (payload = nothing,
            extent = _bboxin(Rasters.EPSG(4326), tcrs, active))
end

# Materialising a `ConstructedSpec` mask yields a `Bool` array on its own real grid, so its extent comes
# free from that array's dims - and materialising here rather than later avoids reading the (possibly
# global) source data a second time.
function _preparemask(active::ConstructedSpec, tcrs)
    # Unwrapped **here and only here**: the combine's contract is that it hands back a raster, and
    # everything downstream of this point works in plain arrays. The package owns the wrapper, so
    # reaching through it internally is free; what the contract buys is that *user* code never has to.
    mask = _materialiseconstructed(active).array
    return (payload = mask, extent = _dimsextent(mask, tcrs))
end

function _preparemask(active, tcrs)
    return error("unrecognised `within` argument of type $(typeof(active)); use nothing, a " *
                 "Matrix{Bool}, a LatLong box, or a mask spec (CircleMaskSpec/ShapeSpec/ConstructedSpec).")
end

# A synthetic unitless target `Rasters.Raster` in `crs`, covering the unitful bounds
# `(ylo..yhi, xlo..xhi)` (given in `crs`'s own coordinate unit) in square cells of side `cellside`.
# This is the `size =` override's grid - a uniform step the reference's own grid may not have. Only
# ever built for a *projected* `crs` (`_targetcrs` fails closed otherwise), so `cellside` is a real
# length in the target's own unit and needs no degree conversion: the old degree-only
# `_wgs84template` and its 111.32 km/° approximation are gone.
function _crstemplate(crs, e::Extents.Extent, cellside)
    ylo, yhi, xlo, xhi = _extentvalues(e)
    u = _crsunit(crs)
    # **`Start`, and the lookup values are cell *edges*** - the convention `_syntheticyx`/`_sizedyx`
    # already use, so there is now one grid convention rather than two. It is not cosmetic: `ylo`
    # and `xlo` are the data's own outer **edges** (`_axisbounds`), so labelling them as *centres*
    # made the template straddle the data by half a cell at the origin. Measured on the 9 × 9 /
    # 2.5 km fixture at 10 km: the first cell went from **50 %** covered to **100 %**, and the grid
    # stopped overhanging the data by 5 km at the bottom. Nothing downstream reads the raw values
    # and assumes what they mean any more (`_cellintervals`/`_cellcentres`), which is what makes this
    # a free choice - see the note on locus-blindness above.
    start = DimensionalData.Lookups.Intervals(DimensionalData.Lookups.Start())
    fwd = DimensionalData.Lookups.ForwardOrdered()
    # The step is `cellside` *exactly* - cells are square and of the requested side, which is the
    # whole point of asking for a metric `size`. The cell count therefore comes from covering the
    # span (`ceil`, so the requested extent is always fully enclosed rather than clipped) and the
    # far edge may overshoot slightly; deriving the step from the span instead would silently hand
    # back not-quite-`cellside` cells. **No `+ 1`**: these are cell *starts*, so `n` of them cover
    # `n * cellside` and `ceil` already encloses the span. The `+ 1` belonged to the `Center` reading,
    # where an n-cell run of centres covers only `(n - 1) * cellside`; kept under `Start` it added a
    # spurious wholly-empty row and column that the recut then had to remove again.
    # **`uconvert`, not `ustrip`** - the template keeps its unit, and naming the
    # unit is what makes the *scale* explicit: a `cellside` in km against bounds in m is brought to
    # one unit before the arithmetic, which a bare strip would silently skip.
    # **Measured bit-identical to the old bare-`Float64` arithmetic**, including for a
    # non-representable step (`0.1°`): a unitful `range` and a naive `Quantity` broadcast both agree
    # exactly with `range(bare) .* u`. So the unit costs nothing in precision and there is no reason
    # to strip, compute and re-dress.
    step = uconvert(u, cellside)
    ny = max(2, ceil(Int, uconvert(NoUnits, (yhi - ylo) / cellside)))
    nx = max(2, ceil(Int, uconvert(NoUnits, (xhi - xlo) / cellside)))
    yvals = collect(range(uconvert(u, float(ylo)), step = step, length = ny))
    xvals = collect(range(uconvert(u, float(xlo)), step = step, length = nx))
    yd = Y(Projected(yvals, sampling = start, crs = crs, order = fwd,
                     span = DimensionalData.Lookups.Regular(step)))
    xd = X(Projected(xvals, sampling = start, crs = crs, order = fwd,
                     span = DimensionalData.Lookups.Regular(step)))
    return Rasters.Raster(zeros(length(yvals), length(xvals)), (yd, xd))
end

# A `ClimateRaster`'s own CRS, read from its `Y` dimension's lookup. NB `Rasters.crs` returns
# `nothing` for a bare `AbstractDimArray` (it is a `Rasters.Raster` method) but works on a *dimension*
# - the same access `_striplength` already uses - so the CRS must be taken from the dim, not the array.
_rastercrs(raster::ClimateRaster) = Rasters.crs(dims(raster.array, Y))

# Whether two CRSs are the same, compared by normalised WKT (`ArchGDAL.importCRS`/`toWKT` puts any
# CRS representation - a bare `EPSG` code or a WKT string - into one comparable form, the same
# normalisation `_crsunit` relies on). A missing CRS (`nothing`) only ever equals another missing one.
function _samecrs(a, b)
    (isnothing(a) || isnothing(b)) && return isnothing(a) && isnothing(b)
    return ArchGDAL.toWKT(ArchGDAL.importCRS(a)) ==
           ArchGDAL.toWKT(ArchGDAL.importCRS(b))
end

# The target grid's CRS, by the staged rule: an explicit `crs` wins; else a single projected CRS
# among the inputs is adopted (so a British National Grid layer combined with WGS84 climate keeps
# the *projected* grid, not the degree one); else the reference's own CRS. A physical `size` needs a
# projected target - square metric cells do not exist on a degree grid - so if the resolved CRS is
# geographic we **fail closed** rather than reviving the 111.32 km/° approximation, and name a
# concrete CRS (`_crsadvice`) in the error so the fix is one paste away.
function _targetcrs(regimes::Tuple, crs, size)
    crss = [_rastercrs(r) for r in regimes]
    resolved = if !isnothing(crs)
        crs
    else
        projected = [c for c in crss if _isprojectedcrs(c)]
        unique_projected = isempty(projected) ? projected :
                           [projected[i]
                            for i in eachindex(projected)
                            if all(j -> !_samecrs(projected[i], projected[j]),
                                   1:(i - 1))]
        length(unique_projected) == 1 ? only(unique_projected) : first(crss)
    end
    if !isnothing(size) && !_isprojectedcrs(resolved)
        here = _extentof(_extrema2(_latvals(first(regimes)),
                                   _longvals(first(regimes)))...)
        error("`cellsize = $size` asks for grid cells of a fixed physical side, but the target grid " *
              "is geographic (° coordinates), where a cell's physical size varies with latitude. " *
              "Pass a projected `crs` to build a genuinely metric grid - " *
              "$(_crsadvice(here)) - or omit `cellsize` to keep " *
              "the data's own native resolution.")
    end
    return resolved
end

# The (min, max) of each of two coordinate vectors, flat, for `_extentof` to name - the shape the
# grid rules want bounds in.
function _extrema2(ys, xs)
    return minimum(ys), maximum(ys), minimum(xs), maximum(xs)
end

# A raster's own north-south cell step, in its own coordinate unit - the "cell size" the resolution
# rules below compare, and what `_uniformcellside` reports.
function _rastercellstep(raster::ClimateRaster)
    # The dimension's **declared** step, not one differenced out of the coordinates. Differencing
    # gives a subtly different `Float64` depending on where in the grid you do it - the same
    # WorldClim layer read globally and read cut to Scotland differ in the 13th digit - so two reads
    # of one file would disagree on their own resolution, and the grids built from them would differ
    # by a cell here and there. `_spatialdim` guarantees the span is `Regular`, so this is exact.
    sp = DimensionalData.Lookups.span(dims(raster.array, Y))
    sp isa DimensionalData.Lookups.Regular && return abs(sp.step)
    lat = _latvals(raster)
    return length(lat) > 1 ? abs(lat[2] - lat[1]) : nothing
end

# The layers' own agreed resolution - the fallback when neither an explicit `cellsize` nor an
# alignment layer settles it. Only the layers actually in the target CRS may vote (a layer in another
# CRS states its step in different units, so it is not comparable), and any disagreement **fails
# closed** with an explicit `cellsize` requirement rather than silently adopting whichever layer
# happens to come first.
function _targetcellsize(regimes::Tuple, tcrs)
    steps = [s
             for s in (_rastercellstep(r)
                       for r in regimes if _samecrs(_rastercrs(r), tcrs))
             if !isnothing(s)]
    isempty(steps) && return nothing        # nothing to vote; the caller measures across instead
    all(==(first(steps)), steps) ||
        error("the input layers in the target CRS disagree on cell size ($(join(unique(steps), ", "))), " *
              "so there is no single native resolution to adopt - pass an explicit `cellsize` to " *
              "choose one.")
    return first(steps)
end

# One of `raster`'s own cells, re-expressed as a side length in `tcrs`'s units.
#
# Needed because a step is only meaningful in the CRS it was measured in: asking to put WGS84 climate
# onto British National Grid is entirely reasonable, but the layer states 0.0833° and the target wants
# metres, and nothing in the unanimity rule can bridge that. Rather than refuse, the cell is
# *measured* across the projection - a representative cell near the layer's middle is transformed and
# its area-preserving side taken (`sqrt(ns × ew)`, the same convention `_cellsize` uses throughout),
# so the answer reflects the projection's real local distortion instead of a nominal
# degrees-to-metres constant of the kind this engine deliberately removed.
#
# Necessarily approximate: a projection's scale varies across the grid, so a cell measured at the
# centre is not every cell. It is a *starting* resolution, announced as such, and any user who needs
# an exact one passes `cellsize`.
function _stepacross(raster::ClimateRaster, tcrs)
    lat, long = _latvals(raster), _longvals(raster)
    (length(lat) > 1 && length(long) > 1) || return nothing
    dlat, dlong = lat[2] - lat[1], long[2] - long[1]
    i, j = cld(length(lat), 2), cld(length(long), 2)
    y0, x0 = lat[i], long[j]
    y1, x1 = y0 + dlat, x0 + dlong
    cell = _bboxin(_rastercrs(raster), tcrs,
                   _extentof(min(y0, y1), max(y0, y1), min(x0, x1),
                             max(x0, x1)))
    return sqrt(abs(cell.Y[2] - cell.Y[1]) * abs(cell.X[2] - cell.X[1]))
end

# The resolution to adopt when none was given: the layers' agreed native step where any of them is in
# the target CRS, else the finest layer's own cell measured across the projection. Returns the size
# and where it came from, so the study area can announce which of the two happened.
function _inferredcellsize(regimes::Tuple, tcrs)
    agreed = _targetcellsize(regimes, tcrs)
    isnothing(agreed) || return agreed, AgreedByAllLayers()
    # No layer shares the target CRS. Measure the *finest* available cell, matching `_choosealign`'s
    # preference for the layer carrying the most detail.
    steps = [s
             for s in (_stepacross(r, tcrs) for r in regimes) if !isnothing(s)]
    isempty(steps) &&
        error("cannot infer a cell size: no input layer states a resolution that can be measured in " *
              "the target grid's CRS. Pass an explicit `cellsize`.")
    return minimum(steps), MeasuredAcrossProjection()
end

# The overlap of two extents, erroring when they miss each other or meet only along an edge.
#
# The intersection itself is `Extents.intersection` -- the arithmetic is not ours to write, and
# `Extents` is already the vocabulary here. What stays is the domain rule: a *touching* overlap has
# zero width and is no more use than a miss, which the library has no reason to know.
function _intersectbounds(a::Extents.Extent, b::Extents.Extent, what)
    both = Extents.intersection(a, b)
    (isnothing(both) || both.Y[1] ≥ both.Y[2] || both.X[1] ≥ both.X[2]) &&
        error("$what do not overlap, so there is no grid to build: got $a and $b.")
    return both
end

# Whether a built habitat's grid can be simulated on - the enforcement half of the study area's
# `:geographic` warning.
# A *synthetic* grid carries no CRS at all (plain cell indices) but does carry a genuine metric
# `size`, so it is fine; a **geographic** CRS is not, because dispersal (`genlookups`) uses the single
# scalar `regime.size` as though every cell were that size, while a degree grid's cells genuinely
# change physical size with latitude. Checked at simulation-assembly time (`build_ecosystem` and
# `Ecosystem`'s own constructor) rather than at `GridHabitat`, so a geographic habitat can still
# be built and inspected/plotted - it just cannot be run.
function _checksimulatable(habitat::AbstractHabitat)
    # **Asks the requirement, not a proxy for it.** Testing *is the CRS geographic?* instead
    # which happens to coincide with *is there a uniform metric cell size?* on every area kind we
    # build - but only happens to. Reading the grid makes the refusal self-maintaining: anything that
    # cannot state one metric cell size is refused, whatever its CRS claims.
    native = getcellsizes(habitat.regime)
    isnothing(native) &&
        return error("This environment's regime has no coordinates, so there is no cell size to " *
                     "disperse across. Build it on a `StudyArea`, or give the layer real `(Y, X)` " *
                     "dims - a bare `NoLookup` grid cannot say how far apart its cells are.")
    # The grid must be metric in its OWN frame. Asking for a length unit would convert a degree grid
    # rather than refusing it, so the check is on what the coordinates already are - one call, and
    # the two failures stay distinguishable, where a dimension argument collapsed both to `nothing`.
    s = (y = first(native.y), x = first(native.x))
    s.y isa Unitful.Length ||
        return error("This environment is on a geographic (° coordinate) grid, which cannot be simulated: " *
                     "dispersal assumes one uniform cell size, but a degree grid's cells change physical " *
                     "size with latitude. Rebuild it on a projected grid by passing `crs` (and `cellsize`) " *
                     "to the `StudyArea` - the warning issued when the area was decided names a suitable " *
                     "CRS for this extent.")
    # **Non-square cells are refused here for the same reason a geographic grid is**: dispersal has
    # one lookup table for the whole grid, so it needs one cell size, and `sqrt(dy*dx)` would hand it
    # a number that is right for neither axis. `≈`, because the tolerance absorbs the last-digit
    # drift a *differenced* step carries (see `_rastercellstep`) - not genuinely rectangular cells.
    s.y ≈ s.x ||
        return error("This environment's cells are not square ($(s.y) × $(s.x)), which cannot be " *
                     "simulated: dispersal builds one lookup table for the whole grid and so needs " *
                     "a single cell size. Pass an explicit square `cellsize` to the `StudyArea` to " *
                     "re-grid onto one. (Inspecting and plotting such a grid is fine - only " *
                     "simulating on it is not.)")
    return nothing
end

# Wrap already-sampled, already-canonicalised values as a regime layer.
#
# **Reached from one place** - `_applyrole`, which `materialise` calls and which the builder therefore
# calls too. A second builder-side path would be a function whose whole body is "decide `categorical`,
# sample, canonicalise, wrap", which is exactly what `_materialisefield` and
# `_applyrole(::Type{Condition})` do between them: two spellings of one composition, agreeing by
# inspection rather than by construction.
#
# A 3-D read - a monthly stack, e.g. `SourceSpec(WorldClim{Climate}, :wind, month = 1:12)` -
# becomes a layer holding one slice at a time and carrying the stack as its own change. The layer
# itself is 2-D either way: which slice is current is a property of *elapsed time*, which the change
# knows and the layer does not.
# **No `csize` argument.** `values` carries its own coordinates, so the layer derives its cell size
# from them with `_derivecellsize` - one source of truth, rather than a number threaded down beside
# the grid it is supposed to describe.
function _asregime(values, categorical::Bool, axis::Type{<:NicheAxis})
    categorical &&
        return _reaxis(CategoricalRegime(values, NoLayerChange()), axis)
    ndims(values) == 2 &&
        return _reaxis(ContinuousRegime(values, NoLayerChange()), axis)
    return _reaxis(_setseries!(ContinuousRegime(_firstslice(values),
                                                NoLayerChange()), values), axis)
end

# The concrete supply type for an axis, refusing an axis that is not a resource.
#
# **This is the wind-speed fix, and it works by making the mistake unrepresentable.** The supply
# path picking a type from the *dimension* of the value would mean `UniformSpec(3.0m/s, axis =
# WindSpeed)` silently built a `WaterSupply` - `m/s` and `mm/day` are both `𝐋𝐓^-1`, and nothing ever
# looked at the axis that had been declared right there. Asking the axis cannot get that wrong: an
# axis that declares no `supplytype` is not a resource, and saying so is the whole of the fix.
function _supplytype(axis)
    T = supplytype(axis)
    isnothing(T) && _notaresource(axis)
    return T
end

# A bare value cannot be a supply, because nothing says what it measures. Refused with the remedy
# rather than left to fail as a `MethodError` on `_specaxis`: this is the form that would work by
# guessing from the unit, so someone will still write it.
function _specaxis(x::Union{Number, AbstractArray})
    return error("a supply given as a bare value ($(typeof(x))) has no niche axis, and its meaning " *
                 "cannot be taken from its unit - `m/s` and `mm/day` are the same dimension. " *
                 "Wrap it in a spec that names one, e.g. `UniformSpec(value, axis = SolarRadiation)`.")
end

# The declared niche axis of a data source element: a `SourceSpec` already carries its own
# (resolved from the shipped layer table at construction time, unless overridden - see its
# constructor); a bare `ClimateRaster` has no code left to look up, so `NicheAxis`.
_specaxis(::SourceSpec{A}) where {A} = A

_specaxis(spec::Tuple) = _sourcepairnotaspec(spec)

# A raster read from a source knows its layer, so its axis is a catalogue lookup like any other -
# only a raster built by hand (no `code`) genuinely has none, and that is what `NicheAxis` records.
function _specaxis(raster::ClimateRaster{S}) where {S}
    (isnothing(raster.code) || raster.code isa AbstractVector) &&
        return NicheAxis
    return something(layeraxis(S, raster.code), NicheAxis)
end

# A `ConstructedSpec` carries the niche axis declared at construction (`NicheAxis` by default).
_specaxis(spec::ConstructedSpec) = spec.axis

# Wrap a sampled supply layer as a supply: `cancel` converts the raw per-area rate (at any native
# time unit) to an absolute per-cell one against `cellarea`, stated in the axis's canonical unit, and
# the **axis** picks the supply type - never the value's dimension, on either count. A
# monthly (3-D) stack becomes a supply holding one slice at a time, carrying the stack as its
# change - exactly as `_asregime` does for a regime.
#
# This is also why there is no separate time-varying supply **type**: such a family would need one
# member per resource, and would refuse any resource that had not been given one. A supply that
# varies in time is the same type as one that does not, so every resource can have one.
function _wrapsupply(out, cellareas, axis)
    T = _supplytype(axis)
    abs = cancel.(out, cellareas, axis)
    ndims(abs) == 2 && return T(abs)
    return _setseries!(T(_firstslice(abs)), abs)
end

# Resolve a `ShapeSpec.path` to a local filesystem path: an already-local path passes
# through; a `CachedAsset` (a URL) is downloaded into its cache (if not already there) here, at
# materialisation time, not at `ShapeSpec` construction.
_resolvepath(path::AbstractString) = path

_resolvepath(asset::CachedAsset) = assetpath(asset)

# --- Selecting a named region ---------------------------------------------------------------------
#
# Turning a name into geometry, shared by everything that answers a question about a named region.
# Whatever asks - a bounding box read from the shipped table, or a mask built onto a real grid -
# comes through here, so a box and the shape it claims to describe cannot disagree.

# The projection component areas are measured on: NSIDC EASE-Grid 2.0 Global, an ellipsoidal
# cylindrical equal-area system covering the whole globe.
#
# It must be reached through `_gdalcrs`, as every transform here is. GDAL 3 gives EPSG:4326 its
# authority axis order, latitude first, so a transform built with `importEPSG` reads each coordinate
# pair the wrong way round - which is not an error, merely a different piece of ground: it makes
# Madagascar 27% too small and France 33% too large.
const _EQUALAREACRS = Rasters.EPSG(6933)

# One connected piece of a region, with what is needed to order and bound it.
#
# Unlike `_ShapePart` this leaves its geometry field abstract. These are built once when a region is
# resolved rather than per cell, so the dynamic access costs nothing worth the concrete type.
const _RegionPart = @NamedTuple{geometry::ArchGDAL.IGeometry,
                                envelope::ArchGDAL.GDAL.OGREnvelope,
                                area::typeof(1.0km^2)}

# Every geometry in `level`'s file whose naming attribute equals `value`.
#
# Matching is case-insensitive because the physical file's own names are inconsistently cased -
# `ALLEGHENY PLATEAU` sits beside `Adelie Coast` - so requiring the source's spelling would make
# some names unreachable in practice. Geometries are cloned: an uncloned one borrows its feature's
# storage and does not outlive the dataset.
function _selectfeatures(level::NaturalEarthLevel, value::AbstractString)
    dataset = ArchGDAL.read("/vsizip/" * assetpath(_nesource(level)))
    lyr = ArchGDAL.getlayer(dataset, 0)
    wanted = lowercase(value)
    geoms = ArchGDAL.IGeometry[]
    for feature in lyr
        _fieldmatches(feature, level.field, wanted) || continue
        isnothing(level.within) ||
            _fieldmatches(feature, level.within.first,
                          lowercase(level.within.second)) || continue
        push!(geoms, ArchGDAL.clone(ArchGDAL.getgeom(feature)))
    end
    return geoms
end

# Whether `feature`'s `field` holds `wanted`, which must already be lowercased. A field absent from
# the layer gives an index of -1 rather than raising, and a null cell gives `nothing`; neither
# matches anything.
function _fieldmatches(feature, field::AbstractString, wanted::AbstractString)
    i = ArchGDAL.findfieldindex(feature, field)
    (isnothing(i) || i < 0) && return false
    value = ArchGDAL.getfield(feature, i)
    return !isnothing(value) && lowercase(string(value)) == wanted
end

# Repair a geometry GEOS would refuse to operate on.
#
# GDAL exposes no `makevalid`, and one of Natural Earth's 258 country outlines (Egypt) is invalid.
# Buffering by zero is the standard substitute and on that geometry restores validity with the area
# unchanged to seven figures.
_repairgeom(g) = ArchGDAL.isvalid(g) ? g : ArchGDAL.buffer(g, 0)

# Merge `geoms` into their connected components, largest first, each with its envelope and its area.
#
# Dissolving before splitting is what makes this answer about *ground* rather than about features:
# neighbouring countries that share a border merge into one landmass, so the largest component of a
# continent is its mainland and not merely its largest country.
function _dissolve(geoms)
    isempty(geoms) && return _RegionPart[]
    merged = reduce((a, b) -> ArchGDAL.union(_repairgeom(a), _repairgeom(b)),
                    geoms)
    parts = map(_components(merged)) do g
        return _RegionPart((g, ArchGDAL.envelope(g), _equalarea(g)))
    end
    return sort!(parts, by = p -> p.area, rev = true)
end

# The connected pieces of a dissolved geometry. Merging polygons that touch gives a single polygon;
# merging scattered ones gives a multipolygon whose members are the pieces.
function _components(g)
    ArchGDAL.getgeomtype(g) == ArchGDAL.wkbMultiPolygon || return [g]
    return [ArchGDAL.clone(ArchGDAL.getgeom(g, i))
            for i in 0:(ArchGDAL.ngeom(g) - 1)]
end

# The area of a lat/long geometry, measured on an equal-area projection.
function _equalarea(g)
    projected = ArchGDAL.clone(g)
    ArchGDAL.createcoordtrans(_gdalcrs(Rasters.EPSG(4326)),
                              _gdalcrs(_EQUALAREACRS)) do ct
        return ArchGDAL.transform!(projected, ct)
    end
    return uconvert(km^2, ArchGDAL.geomarea(projected) * m^2)
end

# The components a coverage asks for. `parts` must already be ordered largest first, as `_dissolve`
# returns them.
_coverageof(parts::AbstractVector{_RegionPart}, ::AllTerritories) = parts

function _coverageof(parts::AbstractVector{_RegionPart}, c::LargestLandmass)
    return parts[1:min(c.count, length(parts))]
end

# One feature of a vector file as `_shape` uses it: the prepared geometry to test cells against, and
# the envelope that says which cells those are. The element type is written out because a layer mixes
# `wkbPolygon` and `wkbMultiPolygon` features, so an inferred one keeps only the field names and
# every access through it becomes a dynamic lookup.
const _ShapePart = @NamedTuple{prepared::ArchGDAL.IPreparedGeometry,
                               envelope::ArchGDAL.GDAL.OGREnvelope}

# Read `spec`'s vector file and prepare every feature's geometry in the target grid's own CRS - the
# work `ShapeSpec` defers from construction to materialise time, mirroring `_read(::SourceSpec)`.
function _shapegeoms(spec::ShapeSpec, tcrs)
    path = _resolvepath(spec.path)
    vpath = endswith(path, ".zip") ? "/vsizip/" * path : path
    dataset = ArchGDAL.read(vpath)
    lyr = ArchGDAL.getlayer(dataset, spec.layer)
    sr = ArchGDAL.getspatialref(lyr)
    dest = _gdalcrs(tcrs)
    # A missing `.prj` (no CRS metadata) gives a null spatial ref; assume already WGS84 lat/long.
    src = sr.ptr != C_NULL ? sr : _gdalcrs(Rasters.EPSG(4326))
    u = _crsunit(tcrs)
    ylo, yhi, xlo, xhi = Inf, -Inf, Inf, -Inf
    # The envelope is taken *before* `preparegeom` - `ArchGDAL.envelope` on a prepared geometry
    # segfaults - and kept, both for the overall extent here and for `_shape` to window each
    # geometry onto the grid.
    geoms = map(lyr) do feature
        geom = ArchGDAL.getgeom(feature)
        ArchGDAL.createcoordtrans(src, dest) do ct
            return ArchGDAL.transform!(geom, ct)
        end
        env = ArchGDAL.envelope(geom)
        ylo, yhi = min(ylo, env.MinY), max(yhi, env.MaxY)
        xlo, xhi = min(xlo, env.MinX), max(xhi, env.MaxX)
        return _ShapePart((ArchGDAL.preparegeom(geom), env))
    end
    extent = isempty(geoms) ? nothing :
             _extentof(ylo * u, yhi * u, xlo * u, xhi * u)
    return geoms, extent
end

# Which indices of `axis` fall within `[lo, hi]`. A raster's Y commonly descends, and `searchsorted`
# given the wrong direction returns an empty range rather than failing, so the direction is read off
# the axis; one monotonic in neither direction falls back to all of it.
function _axiswindow(axis, lo, hi)
    issorted(axis) &&
        return searchsortedfirst(axis, lo):searchsortedlast(axis, hi)
    issorted(axis, rev = true) &&
        return searchsortedfirst(axis, hi, rev = true):searchsortedlast(axis,
                                                                        lo,
                                                                        rev = true)
    return firstindex(axis):lastindex(axis)
end

# Mirrors `_circle`: a cell is active if its centre falls inside any of the shapefile's features.
#
# Each geometry is walked over the cells inside its own envelope. A per-cell envelope test on top of
# that costs more than it saves, a prepared geometry already doing one in C; the win is in not
# visiting the cell at all.
function _shape(geoms, tlat, tlong)
    mask = falses(length(tlat), length(tlong))
    # `geoms` are already in the target's own CRS (`_shapegeoms`), so the coordinates just need
    # their unit stripped - whatever it is (° for a geographic grid, m for a projected one) -
    # rather than being forced to degrees.
    lats, longs = ustrip.(tlat), ustrip.(tlong)
    for g in geoms
        prepared, env = g.prepared, g.envelope
        is = _axiswindow(lats, env.MinY, env.MaxY)
        js = _axiswindow(longs, env.MinX, env.MaxX)
        for i in is, j in js
            mask[i, j] && continue
            mask[i, j] = ArchGDAL.contains(prepared,
                                           ArchGDAL.createpoint(longs[j],
                                                                lats[i]))
        end
    end
    return Matrix{Bool}(mask)
end

# Reproject a bare Bool `DimArray` (e.g. from a `ConstructedSpec` mask) onto `target` via
# `_reproject` - nearest-neighbour (`:near`), since a mask must never blend across a class
# boundary. `_reproject` works in Float64 (GDAL has no Bool dtype), so round-trip through
# 0.0/1.0.
function _samplemask(A::DimensionalData.AbstractDimArray, target)
    r = Rasters.Raster(Float64.(A))
    out = _reproject(r, target, method = :near)
    return Matrix{Bool}(Array(out) .> 0.5)
end

# `cm`'s centre as a `SpatialLocation` in the target's own coordinates: an explicit `centre` is
# always a WGS84 `LatLong`, so it is transformed into `crs` (`_pointin` - a no-op for a geographic
# target); the default is the grid's own midpoint, already in target coordinates whatever the CRS.
#
# Both branches return a place rather than a bare pair, so a caller cannot get the two components
# the wrong way round - which a `(y, x)` tuple silently permits and a non-square grid would then
# expose, three files away from here.
function _centrein(cm::CircleMaskSpec, tlat, tlong, crs)
    isnothing(cm.centre) &&
        return SpatialLocation((first(tlat) + last(tlat)) / 2,
                               (first(tlong) + last(tlong)) / 2)
    return _pointin(crs, cm.centre)
end

# Great-circle-ish planar distance from each grid cell to the centre, on a **geographic** (degree)
# target: the degree offsets become physical distances via `_side`, the east-west one shrunk by
# cos(latitude) as the meridians converge. Dispatch mirrors `_cellsize`'s own geographic/projected
# split (`GridHabitat.jl`): this generic method is the degree case, the `Unitful.Length` method below
# the projected one.
function _circle(cm::CircleMaskSpec, tlat, tlong, crs)
    centre = _centrein(cm, tlat, tlong, crs)
    mask = falses(length(tlat), length(tlong))
    for (i, y) in enumerate(tlat), (j, x) in enumerate(tlong)
        dlat = _degreelength(y - centre.y)
        dlong = _degreelength(x - centre.x) * cos(centre.y)
        mask[i, j] = hypot(dlat, dlong) <= cm.radius
    end
    return Matrix{Bool}(mask)
end

# The **projected** (length-coordinate) target: the grid is already a metric plane, so the distance is
# a plain Euclidean one - no degree->length conversion and no cos(latitude) convergence term, making
# this both simpler and more accurate than the geographic approximation above.
function _circle(cm::CircleMaskSpec, tlat::AbstractVector{<:Unitful.Length},
                 tlong::AbstractVector{<:Unitful.Length}, crs)
    centre = _centrein(cm, tlat, tlong, crs)
    mask = falses(length(tlat), length(tlong))
    for (i, y) in enumerate(tlat), (j, x) in enumerate(tlong)
        mask[i, j] = hypot(y - centre.y, x - centre.x) <= cm.radius
    end
    return Matrix{Bool}(mask)
end

# The cells a layer genuinely has data for - its non-`NaN` cells, and for a collection the cells live
# in *every* sub-layer. Deliberately role-generic (`AbstractLayer`, not `AbstractRegime`): a supply's
# gaps mark cells unusable exactly as a regime's do, and the layer types differ only in their `Role`
# parameter, so one rule serves both rather than two that could drift apart.
function _coverage(layer::AbstractLayer)
    # A layer's `matrix` is always `(Y, X)`, holding the values current now - for a time-varying
    # layer, its first stored slice, since coverage is decided before any time has elapsed. That is
    # the same first-slice convention `_nanactive` (this file) uses, but it now falls out of
    # the layer's own shape rather than needing a branch.
    #
    # `isnan` is `false` for the integer niche codes a `CategoricalLayer` holds, so a categorical
    # layer correctly reports full coverage without needing its own branch.
    return Matrix{Bool}(.!isnan.(layer.matrix))
end

function _coverage(layer::LayerCollection)
    covered = _fold(_coverage, values(layer)) do a, b
        return a .& b
    end
    return Matrix{Bool}(covered)
end

# Rasterise a *prepared* mask (`_preparemask`'s `payload`) onto the target grid as a `Matrix{Bool}` -
# one method per recognised payload kind, dispatched on its type rather than branched on at runtime.
# This is only the mask itself; the caller ANDs it with the layers' real coverage (`_coverage`).
#
# **There is deliberately no `::Nothing` method.** "No mask at all" is answered by
# `_rastermaskonly(::Nothing, ...)`, which returns `trues` - the study area's own question. The method
# Returning the regime's *coverage* here instead would have a caller AND it in a second time, which
# is a double count that nothing downstream can undo.
function _rastermask(payload::DimensionalData.AbstractDimArray, regime, target)
    return _samplemask(payload, target)
end

function _rastermask(payload::CircleMaskSpec, regime, target)
    yx = _cellcentres(target)
    return _circle(payload, yx.lat, yx.long, Rasters.crs(target))
end

# `_preparemask(::ShapeSpec, ...)` already read and reprojected the geometries, so the payload is the
# vector of prepared geometries and their envelopes - nothing is re-read here.
function _rastermask(payload::AbstractVector, regime, target)
    yx = _cellcentres(target)
    return _shape(payload, yx.lat, yx.long)
end

function _rastermask(payload::AbstractMatrix{Bool}, regime, target)
    dims_ = Base.size(target)
    Base.size(payload) == dims_ ||
        error("`active` is $(Base.size(payload)) but the grid is $dims_")
    return Matrix{Bool}(payload)
end

# Materialise a `ConstructedSpec` to a raster/array: read each source spec (`_asraster`) and apply
# `combine` - nullary when there are no sources (the thunk produces the layer directly).
function _materialiseconstructed(spec::ConstructedSpec)
    return _combined(spec.combine(map(_asraster, spec.layers)...), spec)
end

# **One contract for every combine: rasters in, a raster out** - whether the spec ends up used as a
# layer or as a mask, which is a decision made *elsewhere* and cannot be known here. A mask is simply
# a raster whose element type is `Bool`.
#
# Making it depend on the later use would mean a layer combine returning a `ClimateRaster` while a
# mask combine returned a bare array, with neither stated where the combine is written. That would put
# an array type into **user** code, which is exactly what the package intends to stay free to change.
#
# Since a raster now broadcasts, satisfying this takes no extra work: `lc .!= code` and
# `sum(bands)` are already rasters, so the natural way to write a combine is the correct one.
function _combined(out, spec::ConstructedSpec)
    out isa ClimateRaster ||
        error("a `ConstructedSpec` combine must return a `ClimateRaster`, but this one returned a " *
              "$(typeof(out)). A raster broadcasts, so operating on the layers directly gives one " *
              "back - write `lc .!= code` rather than `lc.array .!= code`, and `sum(bands)` rather " *
              "than `ClimateRaster(T, sum(b -> b.array, bands))`.")
    return out
end

# --- Mask-led extent: recut to the true active footprint -------------------------------------------
# What a mask asks for but cannot get - cells it marks active where no layer has data, and ground it
# reaches beyond the layers altogether - is classified by `_analyse` and reported by the `StudyArea`
# that decided the grid (the `:mask_lost` and `:mask_clamped` problems).

# The `(Y, X)` index ranges spanning every active cell - the recut. Applied to the *template*, so the
# habitat is then rebuilt on the cropped grid: that keeps regime/supply/active mutually consistent by
# construction (they are all resampled onto one grid) rather than by cropping three arrays in step, and
# leaves each one a single resample straight from its source data.
# Widen a recut range to at least two cells (within the `n` available), because a single-cell axis has
# no derivable step: `_cellsize` reads the spacing from the coordinates themselves, so a 1×N grid would
# fail there. Keeping one spare cell is also more useful than a degenerate grid - dispersal needs
# neighbours. If only one cell exists in total the grid was already degenerate before any recut.
function _atleast2(r, n)
    (length(r) >= 2 || n < 2) && return r
    first(r) > 1 && return (first(r) - 1):last(r)
    return first(r):(last(r) + 1)
end

# The smallest row and column range containing every active cell - what the grid is cropped to once
# the mask is known, so that a study area is not carrying rows and columns nothing lives in.
function _activerange(active, simulate_safely::Bool = false)
    rows = findall(any(active, dims = 2)[:, 1])
    cols = findall(any(active, dims = 1)[1, :])
    (isempty(rows) || isempty(cols)) &&
        error(_noactivemessage(simulate_safely))
    return first(rows):last(rows), first(cols):last(cols)
end

# **Under `simulate_safely` the cause is usually the flag, not the mask.** The general message
# blames the `within` extent and the CRS, which is right when the data and the area genuinely miss
# each other - but with data smaller than a single cell every cell is *partly* covered, all of them
# are refused, and reading about the CRS sends the user looking in the wrong place. So the flag's own
# case names the flag, and names the remedy.
function _noactivemessage(simulate_safely::Bool)
    simulate_safely ||
        return "no cell is active: the mask and the layers' real coverage do not overlap anywhere. " *
               "Check the study area's `within` extent against the data, and its `crs` if given."
    return "no cell is active: with `simulate_safely = true` (the default) a cell is only simulated " *
           "when every layer covers the whole of it, and no cell of this grid is wholly covered. " *
           "That usually means the cells are larger than the ground the data describes - ask for a " *
           "smaller `cellsize`, widen the data, or pass `simulate_safely = false` to simulate " *
           "partly covered cells as well."
end

# Re-wrap a raster's data with a unit, preserving the axes and the source type (so
# `iscategorical` still dispatches correctly). `NoUnits` is a genuine no-op that leaves
# the data as bare numbers, which is what a categorical (class-code) regime wants.
function _attachunit(raster::ClimateRaster{S}, u) where {S}
    A = raster.array
    return ClimateRaster(S, DimArray(A.data .* u, dims(A)), raster.code)
end

# Read a `SourceSpec` into a unit-attached `ClimateRaster` (the eager step deferred by the
# lazy descriptor). A `nothing` code is a whole-dataset spec - read all layers as one multi-band
# raster (`read`'s own default), dimensionless. The spec's own `readkw` (e.g. `month = 1:12`, given
# when it was constructed) forward to `read`, with any keyword passed here overriding them, so a
# caller can still refine a spec's read without rebuilding it.
function _read(sl::SourceSpec; kw...)
    readkw = merge(sl.readkw, NamedTuple(kw))
    # Can these layers honestly share one array? Refused on two counts, for the same underlying
    # reason: an array has **one** eltype and gets **one** resample method, so its layers must agree
    # on both.
    #
    # * **Unit** - a stack of °C and mm and m s^-1 comes back as bare magnitudes, quietly inviting a
    # combine to compare them. Four of the seven shipped datasets are like this.
    # * **Categorical-ness** - one array is resampled one way. A stack mixing class codes with
    # measurements can only be interpolated (meaningless for the codes) or taken by nearest class
    # (lossy for the measurements); there is no right answer, so it must not be built.
    #
    # The second check is *not* implied by the first. Measured across `CHELSA{BioClimPlus}`: every
    # unit group is valuetype-uniform **except the dimensionless one**, which holds all three - so a
    # stack of `kg0` (class codes) and a dimensionless continuous layer passes the unit test and
    # would otherwise reach the resampler with no correct answer available.
    #
    # Both are checked here, at the point a single array is actually demanded and before anything is
    # downloaded, rather than at construction, where they would rule out `WorldClim{BioClim}` and
    # three other shipped datasets outright.
    if sl.code isa AbstractVector
        onelayer = "Name the layers you want instead - `SourceSpec($(sl.source), :code)` for " *
                   "one, or pass a codes vector to `ConstructedSpec`, which reads each on its " *
                   "own terms."
        us = unique(layerunit(sl.source, c) for c in sl.code)
        length(us) == 1 ||
            error("`SourceSpec($(sl.source))` covers $(length(sl.code)) layers with " *
                  "$(length(us)) different units ($(join(us, ", "))), so they cannot be read " *
                  "into one array without losing them. " * onelayer)
        # Called for its refusal: the vector method throws on a stack that mixes class codes with
        # measurements, which is the rule stated once rather than restated here.
        iscategorical(sl.source, sl.code)
    end
    raw = read(sl.source, sl.code; readkw...)
    # Tag the materialised raster with the layer it came from. This is the single point at which a
    # spec becomes data, and the only place that still knows the code - after this the raster is
    # passed around on its own and the shipped table can no longer be consulted for it, which is
    # exactly why `iscategorical` answers `false` for a raster carrying no code.
    data = _applyperiod(raw.array, sl, readkw)
    return _attachunit(ClimateRaster(sl.source, data, sl.code), sl.unit)
end

# Divide a freshly-read array by the interval each of its slices accumulated over, so the values are
# rates before any unit is attached.
#
# **This is the step the whole layer-units subproject exists for.** A monthly total divided by a
# fixed 30.4375-day month is wrong for every real month - 7.7% low for February. Dividing by the
# slice's *own* month makes the rate honest, and because every divisor is a time, all twelve slices
# still share one unit (`sl.unit`, already the rate under `layerrate`); only their values differ.
#
# Which months these slices are comes from `readkw`, never from the array's `Ti` axis - a partial
# read (`month = 2:4`) builds an axis labelled 1:3, so the axis cannot say. Defaults to all twelve,
# which is what an unqualified monthly read returns.
#
# A layer with no period, a per-cell period, or one whose canonical reading is the accumulated amount
# itself (a heat sum) is returned untouched.
function _applyperiod(array, sl::SourceSpec, readkw::NamedTuple)
    sl.code isa AbstractVector && return array          # multi-band: no single period applies
    rec = layerinfo(sl.source, sl.code)
    months = get(readkw, :month, axes(array, ndims(array)))
    divisors = _readdivisors(rec, months)
    isnothing(divisors) && return array
    divisors isa Number && return array ./ divisors
    length(divisors) == size(array, 3) ||
        error("`$(sl.code)` accumulates over one calendar month per slice, but the read returned " *
              "$(size(array, 3)) slices for $(length(divisors)) month(s) - the two must agree, or " *
              "the wrong month's length would be divided into the wrong slice's values.")
    return array ./ reshape(collect(divisors), 1, 1, :)
end

# **A raster is not a spec, and is refused as one.** A `ClimateRaster` holds values and a layer
# *code*; a niche axis is resolved from that code through the shipped catalogue, so a raster carrying
# no code - every derived one - silently becomes `NicheAxis`, with no keyword anywhere to correct
# it. That is the one remaining way a layer's meaning could go undeclared, which is what this whole
# subproject exists to prevent. Wrapping it in a spec is not ceremony: the spec is *where the axis
# is stated*, which is exactly the thing a raster cannot say.
#
# Shared by every entry point (`_asraster` with and without a cache, and `materialise` itself)
# so the message cannot depend on which one a caller happened to reach.
function _rasternotaspec(raster)
    return error("a `$(nameof(typeof(raster)))` is not accepted as a regime, supply or study-area " *
                 "layer: it carries values and a layer code, but no niche axis, so what it means " *
                 "could only be guessed. Name the data instead - `SourceSpec(source, code)` - or, " *
                 "for a raster you already hold, wrap it in a spec that declares the axis: " *
                 "`ConstructedSpec(() -> raster; axis = SomeAxis)`.")
end

# **A bare `(source, code)` pair is not a spec either, and is refused with the spelling that is.**
# **A refusal rather than no method at all**, for the same reason `_rasternotaspec` is one: the form
# reads plausibly, so someone will write it, and the `MethodError` it would otherwise produce names an
# internal function rather than the remedy. The pair says *where* the data is but
# nothing about what it means; a `SourceSpec` is where the axis, the unit and the read options live,
# and it is a single extra word.
#
# A tuple in a `regime`/`supply` position means a **multi-layer environment** - one member per
# layer - which is why the pair form could only ever be a member of one, and why the two forms could
# be told apart by nothing but nesting depth.
function _sourcepairnotaspec(spec)
    return error("a `(source, code)` tuple is no longer accepted as a layer: a tuple `regime` or " *
                 "`supply` describes a **multi-layer** environment, one spec per layer, so a pair " *
                 "inside one could only be told from a layer by how deeply it was nested. Name the " *
                 "data with a spec instead - `SourceSpec($(join(spec, ", ")))` - which is also " *
                 "where its axis, unit and read options live.")
end

# Normalise a regime spec to a `ClimateRaster`: a `SourceSpec` is read and unit-attached as in the
# single-layer source builder. (`ConstructedSpec` layers are always `SourceSpec`s - see
# `_parselayers`; a bare-dataset layer is a whole-dataset `SourceSpec`, read via `_read`.)
_asraster(raster::ClimateRaster) = _rasternotaspec(raster)

_asraster(spec::SourceSpec) = _read(spec)

_asraster(spec::Tuple) = _sourcepairnotaspec(spec)

# ---------------------------------------------------------------------------
# Single-signature GridHabitat (the public spec-based API)
# ---------------------------------------------------------------------------

# A `ConstructedSpec` as a regime/supply layer: materialise it to a raster (its `combine` result).
_asraster(spec::ConstructedSpec) = _materialiseconstructed(spec)

# --- Raster-geometry primitives ----------------------------------------------
# Geometry rather than climate data: what unit a CRS measures in, whether a value is an
# angle, and how to stack layers into one array. Nothing in the submodule used them - the main
# module was reaching in for all three.

# Is this coordinate or cell size an angle (a geographic grid) rather than a length (a projected
# one)? True for `°`, `arcminute`, `arcsecond` and radians alike, so a caller may write
# `cellsize = 30arcsecond` and still get the exact integer-arcsecond grid arithmetic instead of
# silently dropping onto the floating-point path.
#
# Neither half of this test is sufficient alone, which is why it reads oddly. Unitful follows SI
# in modelling an angle as **dimensionless**, so `dimension` only rules out a length's `𝐋` - it
# cannot tell `1°` from a bare `0.5`. Nor does convertibility: `uconvert(arcsecond, 0.5)` does *not*
# fail, it returns 103132.4″, because Unitful reads an unadorned number as radians. An angle is
# precisely the thing that is dimensionless *and yet carries a unit*.
#
# The one theoretical false positive is `u"percent"`, also dimensionless-with-a-unit. It cannot reach
# here: percentages appear in this package as land-cover *values*, never as coordinates or cell sizes.
_isangle(x) = dimension(x) === NoDims && unit(x) !== NoUnits

# The physical unit of a CRS's own coordinate axes: `°` for a geographic CRS (e.g. WGS84), or
# whatever linear unit a projected CRS declares (e.g. metres for British National Grid) - never
# assumed. No high-level ArchGDAL/GDAL accessor exposes "is this geographic/projected" or "what's
# the linear unit" (checked directly against the installed ArchGDAL - genuinely absent at both
# the high-level API and the low-level `GDAL.OSR*` bindings), so this is read directly off the
# CRS's own WKT text: the *last* `UNIT[...]` token, since a projected CRS's `WellKnownText`
# reports its nested geographic CRS's (usually degree) unit first, and its own actual coordinate
# unit last, right before `AXIS`/`AUTHORITY`. Confirmed against a default WGS84 read, a real
# British National Grid (EPSG:27700) file, and a real ERA5-via-`retrieve_era5` netCDF file (whose
# CRS comes back as a bare `GeoFormatTypes.EPSG` code, not `WellKnownText`, since CDS-produced
# netCDF carries no WKT `grid_mapping` - `ArchGDAL.importCRS`/`toWKT` normalises *any* CRS
# representation Rasters/ArchGDAL might hand back to WKT text first, so the same regex covers
# both). A `nothing` CRS (no metadata at all - confirmed this occurs for a netCDF file with no
# `grid_mapping`/CRS attribute) is assumed already WGS84 lat/long, the same fallback `ShapeSpec`
# reading (this file) already uses for a shapefile with no `.prj`.
_crsunit(::Nothing) = °

function _crsunit(crs)
    # A *blank* CRS means the same as no CRS, and must be caught before `importCRS` sees it.
    _isblankcrs(crs) && return _crsunit(nothing)
    wkt = ArchGDAL.toWKT(ArchGDAL.importCRS(crs))
    matches = collect(eachmatch(r"UNIT\[\"([^\"]+)\"", wkt))
    isempty(matches) &&
        error("could not determine the coordinate unit from CRS WKT: $wkt")
    name = lowercase(last(matches).captures[1])
    name in ("degree", "degrees") && return °
    name in ("metre", "metres", "meter", "meters") && return Unitful.m
    return error("unrecognised CRS coordinate unit \"$name\" in WKT - add a mapping in " *
                 "`_crsunit` if this is expected: $wkt")
end

# Combine `aas` (arrays sharing the same leading dims) into one array by stacking them along a new
# trailing dimension `newdim` - a single array is returned unchanged. Used both for per-file reads
# (inputs 2-D `(Y,X)`, output 3-D with `newdim` = `Ti`/`Dim{:layer}`) and for combining several
# already-built layers into one multi-layer array (inputs already `(Y,X,...)`, output with one
# more trailing dimension) - the same operation, generalised to whatever shape the inputs already
# have rather than assuming 2-D.
function _stacklayers(aas::AbstractVector, newdim)
    length(aas) == 1 && return first(aas)
    first_a = first(aas)
    stackdim = ndims(first_a) + 1
    return DimArray(cat(parent.(aas)..., dims = stackdim),
                    (dims(first_a)..., newdim))
end

# A CRS that says nothing - `_crsunit` needs it, so it travelled with it.

# Whether a CRS is present in name only. A file whose projection tag exists but is *empty* - a plain
# GeoTIFF written without a CRS, which `data/Africa.tif` was until it was repaired - comes back from
# Rasters as `WellKnownText("")` rather than `nothing`, so the `::Nothing` fallback above misses it
# and `ArchGDAL.importCRS("")` rejects it with the opaque "Failed to initialize SRS based on WKT
# string (Corrupt data.)". Semantically the two are the same thing (no georeferencing declared), so
# they get the same WGS84 assumption. Only `WellKnownText` can be blank in this way: an `EPSG` code
# holds integers, and a `ProjString` that is empty is malformed rather than absent.
_isblankcrs(crs) = false

function _isblankcrs(crs::Rasters.GeoFormatTypes.WellKnownText)
    return isempty(strip(Rasters.GeoFormatTypes.val(crs)))
end
