# SPDX-License-Identifier: LGPL-3.0-or-later
#
# GETTING RASTER DATA ONTO A GRID. Everything between a file on disk and a value in the right cell:
#
# reading        `_read`, `_asraster`, `_attachunit`, `_applyperiod` - a spec to a `ClimateRaster`
# regridding     `_regrid`, `_blockaggregate`, `_reducer`, `_reproject`, `_targetcrs`, `_cellintervals`
#                - a raster onto a grid, by aggregation of the cells covering each of its cells; the bulk
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
# Data-driven active-area masks are composed with `ConstructedRasterSpec` from a data source plus a
# combine rule; `CircleMaskSpec`/`ShapeSpec` are the synthetic/vector mask specs. The two public
# rules below are reusable building blocks for writing your own combine (`_circle`/`_shapegeoms`+
# `_shape` remain the private geometry helpers for the synthetic/vector masks).

"""
    hasdata(layer)

A [`ConstructedRasterSpec`](@ref) combine returning a `Bool` mask of the cells of `layer` that hold data
(are not missing/`NaN`) - the canonical combine-rule example. Pass it a data source to mask that
source's coverage: `ConstructedRasterSpec(hasdata, WorldClim{BioClim}, 1)`.

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
# reprojecting if `r`'s CRS differs. **Nearest-neighbour only**: each target cell takes the value of
# the source cell its centre falls in, so no value is invented - `_regridsampled`'s point sampling,
# whose aggregation over the samples is what gives the regridding its meaning. A 3-D+ array is
# resampled one 2-D slice at a time - `Rasters.resample` errors directly
# on 3-D input (`GDALError: Too few arguments for '-te'`, confirmed) - and restacked along the
# original axis values, so a real `Ti{DateTime}`/`Dim{:layer}` axis survives unchanged. Returns a
# plain `DimArray` (not a `Rasters.Raster`), matching `ClimateRaster.array`'s existing convention.
function _reproject(r::Rasters.AbstractRaster, target; method = :near)
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
    out = _ontarget(resampled, target)
    return dataunit === NoUnits ? out : out .* dataunit
end

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

# The reducer a coarsening uses: an explicit `fn` as given; otherwise the most frequent class
# for an axis holding class codes, whose mean would name a class nobody observed, and the mean for
# any other. Resolved once, before the aggregate cache key is built, so the key names the reducer
# actually applied.
_reducer(fn, ::Type{<:NicheAxis}) = fn

function _reducer(::Nothing, axis::Type{<:NicheAxis})
    return iscategorical(axis) ? _majorityclass : _meanpresent
end

# **One rule for a block with missing cells, shared by every reducer**: it reduces over the cells
# that carry data, and is absent only where none do. Missing and NaN both count as absent. Whether a
# grid cell is *covered* by a layer is a separate question, answered by the layer having data at the
# cell's centre (see `_coveredgrid`), so a value here never decides activity on its own.
_absent(::Missing) = true

_absent(v::Number) = isnan(v)

_absent(_) = false

_present(block) = (v for v in block if !_absent(v))

_anypresent(block) = any(!_absent, block)

# What an empty reduction is: `missing` where the block can hold one, NaN in the block's own type
# (units included) otherwise, so the result stays in the array it goes back into.
function _absentvalue(block)
    T = eltype(block)
    Missing <: T && return missing
    return convert(T, NaN * oneunit(T))
end

# The mean of the cells that carry data.
function _meanpresent(block)
    _anypresent(block) || return _absentvalue(block)
    return mean(_present(block))
end

# The most frequent value in a block of class codes, ties broken by the smallest code so the answer
# does not depend on iteration order. Read-time and build-time only, so the allocation per block is
# of no consequence. Applied in two stages by `_regridsampled` it is an approximation: the majority
# of block majorities can differ from the majority over every covering cell. An exact class on a much
# coarser grid comes from expanding the codes to one fraction band per class, whose means compose,
# and taking the most frequent class once at the end.
function _majorityclass(block)
    _anypresent(block) || return _absentvalue(block)
    counts = Dict{nonmissingtype(eltype(block)), Int}()
    for v in _present(block)
        counts[v] = get(counts, v, 0) + 1
    end
    best, bestn = first(counts)
    for (v, n) in counts
        (n > bestn || (n == bestn && v < best)) && ((best, bestn) = (v, n))
    end
    return best
end

# Coarsen `A` by a whole factor `f` per spatial axis, each `f × f` block becoming one cell holding
# `fn` of it. The one block-aggregation primitive: the read-time `scale` and `_regrid` both call it,
# so a layer coarsened on read and one coarsened onto the grid are the same computation. Further
# dimensions are left alone, `f == 1` is the identity, and a partial block at the far edge is
# dropped, which is why `_regrid` crops to whole blocks first.
function _blockaggregate(A, f::Integer, fn)
    f == 1 && return A
    return Rasters.aggregate(fn, _anchorsouthwest(_rasterof(A), f),
                             (Y(f), X(f)))
end

# `Rasters.aggregate` wants an `AbstractRaster`; a bare `DimArray` is wrapped, a raster passes.
_rasterof(A::Rasters.AbstractRaster) = A

_rasterof(A) = Rasters.Raster(A)

# Blocks are anchored at the south-west corner, where `_crstemplate` starts the study grid, so the
# partial block a coarsening drops is the northern or eastern one. `Rasters.aggregate` blocks from
# index 1, which is that corner on an ascending axis; a north-first file has its leading remainder
# trimmed here. `_blockrange` widens a windowed read to the same lattice.
function _anchorsouthwest(r, f::Integer)
    for D in (Y, X)
        n = size(r, D)
        offset = _trimoffset(DimensionalData.Lookups.order(DimensionalData.lookup(r,
                                                                                  D)),
                             n, f)
        offset == 0 || (r = r[D((offset + 1):n)])
    end
    return r
end

# How many leading cells a block lattice anchored at the south-west corner skips on an axis of `n`
# cells in blocks of `f`: none on an ascending axis, whose first cell is the south or west edge; the
# remainder on any other, whose first cell is the north or east edge.
function _trimoffset(::DimensionalData.Lookups.ForwardOrdered, n::Integer,
                     f::Integer)
    return 0
end

_trimoffset(::Any, n::Integer, f::Integer) = n % f

# The run of `f · n` source cells (starts, ascending) tiling the `n` target cells from
# `first(targetvals)`, or `nothing` where the target does not start on a source cell boundary or
# reaches past the source. Alignment is `_originaligned`'s test; the run's extent follows from both
# grids being regular. (`_blockrange` in `datasetread.jl` does the same for a windowed read, in the
# file's indices.)
function _tilerange(sourcevals, targetvals, f::Integer)
    (length(sourcevals) >= 2 && !isempty(targetvals)) || return nothing
    sstep = sourcevals[2] - sourcevals[1]
    _originaligned(first(sourcevals), first(targetvals), sstep) ||
        return nothing
    offset = round(Int,
                   uconvert(NoUnits,
                            (first(targetvals) - first(sourcevals)) / sstep))
    i = offset + 1
    j = i + f * length(targetvals) - 1
    (i >= 1 && j <= length(sourcevals)) || return nothing
    return i:j
end

# `raster` on `target`'s grid, every cell an aggregation by `fn` of the source cells covering it, by
# one of two routes. Where `target` is in `raster`'s CRS, its step the same whole multiple `f` of
# `raster`'s on both axes and its cells on `raster`'s cell boundaries - the case `_resamplecost`
# reports as `LayerAggregated(f)`, or `LayerKeptExactly` at `f == 1` - the source is cropped to whole
# blocks and block-aggregated exactly, `f == 1` simply selecting cells. Any other grid goes through
# `_regridsampled`. Nothing is interpolated on either route.
#
# Selecting cells rather than warping them matters even at `f == 1`: a GDAL warp onto an identical
# grid poisons any output cell whose input stencil touches a `NaN`, and which cells that catches
# depends on the extent that was read, so a global and a windowed read of one layer can differ by
# over a hundred coastal cells. Indexing makes the answer invariant by construction.
#
# Source and target are compared with their coordinate units on, so a `km` source matches an `m`
# target of the same ground; stripping one side only would turn an exact selection into a resample,
# which shows in the grid's shape and not in its values.
#
# The reducer runs on bare magnitudes: the unit is taken off before either route and put back on the
# result, since a mean of `°C` is a sum Unitful refuses. `centres = true` asks for the value at each
# target cell's centre alone, which is how `_coveredgrid` decides coverage.
function _regrid(raster::ClimateRaster, target, fn; centres::Bool = false)
    u = _dataunit(raster.array)
    out = _regridbare(raster, _bare(raster.array), target, fn, centres)
    return u === NoUnits ? out : out .* u
end

# The unit an array's values carry, `NoUnits` for plain numbers, and the array without it.
_dataunit(::AbstractArray{<:Unitful.Quantity{T, D, U}}) where {T, D, U} = U()

_dataunit(::AbstractArray) = NoUnits

_bare(A::AbstractArray{<:Unitful.Quantity}) = ustrip.(A)

_bare(A::AbstractArray) = A

function _regridbare(raster::ClimateRaster, A, target, fn, centres::Bool)
    centres && return _regridsampled(raster, A, target, fn, centres = true)
    yd, xd = dims(A, Y), dims(A, X)
    _samecrs(Rasters.crs(yd), Rasters.crs(target)) ||
        return _regridsampled(raster, A, target, fn)
    sy,
    sx = parent(DimensionalData.lookup(yd)),
         parent(DimensionalData.lookup(xd))
    ty,
    tx = parent(DimensionalData.lookup(target, Y)),
         parent(DimensionalData.lookup(target, X))
    (length(sy) >= 2 && length(sx) >= 2 && length(ty) >= 2 && length(tx) >= 2) ||
        return _regridsampled(raster, A, target, fn)
    ry = _stepratio(_lookupstep(yd), _lookupstep(dims(target, Y)))
    rx = _stepratio(_lookupstep(xd), _lookupstep(dims(target, X)))
    (!isnothing(ry.factor) && ry.factor == rx.factor && !ry.finer && !rx.finer) ||
        return _regridsampled(raster, A, target, fn)
    f = ry.factor
    rows = _tilerange(sy, ty, f)
    cols = _tilerange(sx, tx, f)
    (isnothing(rows) || isnothing(cols)) &&
        return _regridsampled(raster, A, target, fn)
    return _ontarget(_blockaggregate(A[Y(rows), X(cols)], f, fn), target)
end

# The route for a grid that is not an aligned whole multiple of the source's cells, in two stages.
# First the source is block-aggregated on its own lattice by the whole part of the step ratio, which
# is exact and leaves it at most twice as fine as the target. Then it is sampled nearest-neighbour
# onto a lattice `k` times finer than `target` and block-aggregated by `k` with `fn`: every fine
# cell carries a real source value, so this is an area-weighted aggregation of the covering cells
# to `k²` samples per target cell, for any reducer, and a finer target comes out as repetition. The
# first stage is what keeps `k` small however coarse the grid; without it a grid eighty times
# coarser than its source would either need a fine lattice of 25 000 cells per target cell or, capped,
# see a fraction of the covering cells. `fn` is applied at both stages, so it must compose: a mean,
# a maximum or a minimum do exactly, the most frequent class only approximately. `A` is `raster`'s
# array with its unit taken off; `raster` is needed for the step across a change of CRS. With
# `centres` the target itself is the sampling lattice and nothing is aggregated: each cell takes the
# source value at its centre.
function _regridsampled(raster::ClimateRaster, A, target, fn;
                        centres::Bool = false)
    centres && return _ontarget(_reproject(Rasters.Raster(A), target), target)
    ratio = _stepratioacross(raster, target)
    p = isnothing(ratio) ? 1 : max(1, floor(Int, ratio))
    coarse = _blockaggregate(A, p, fn)
    k = _oversampling(isnothing(ratio) ? 1.0 : ratio / p)
    sampled = _reproject(Rasters.Raster(coarse), _finetemplate(target, k))
    return _ontarget(_blockaggregate(sampled, k, fn), target)
end

# `A` on `target`'s own `(Y, X)` dims, its other dims kept: both routes end here, so a layer sampled
# either way carries the identical, `Regular`-spanned coordinates.
function _ontarget(A, target)
    yd, xd = _targetyx(target)
    return DimArray(parent(A), map(d -> _replaceyx(d, yd, xd), dims(A)))
end

# The target's own `Y` or `X` in place of a sampled array's; any other dimension is kept.
_replaceyx(::Y, yd, xd) = yd

_replaceyx(::X, yd, xd) = xd

_replaceyx(d, yd, xd) = d

# The target's step over the source's, as a plain number - measured across a change of CRS by
# `_stepacross` - or `nothing` where the source is too small to measure.
function _stepratioacross(raster::ClimateRaster, target)
    tcrs = Rasters.crs(target)
    tstep = _lookupstep(dims(target, Y))
    sstep = _samecrs(_rastercrs(raster), tcrs) ?
            _lookupstep(dims(raster.array, Y)) : _stepacross(raster, tcrs)
    isnothing(sstep) && return nothing
    return uconvert(NoUnits, tstep / sstep)
end

# Fine samples per target cell side for the residual step ratio left after pre-aggregation: twice
# the ratio rounded up, so every contributing source cell is sampled at least twice per axis. The
# residual is below two, so this is 2 or 4; anything larger means the pre-aggregation was skipped.
function _oversampling(residual::Real)
    k = 2 * ceil(Int, residual)
    2 <= k <= 4 ||
        error("an oversampling factor of $k means the source was not pre-aggregated first")
    return k
end

# A template `k` times finer than `target` that tiles it exactly: `k` fine cells per target cell on
# each axis, starting where the target starts, so `_blockaggregate` by `k` lands every block on one
# target cell. Built directly from the target's lookups rather than through `_crstemplate`, whose
# `ceil` over a span could gain or lose a fine row to rounding.
function _finetemplate(target, k::Integer)
    yd, xd = _targetyx(target)
    ys,
    xs = parent(DimensionalData.lookup(yd)),
         parent(DimensionalData.lookup(xd))
    sy, sx = _lookupstep(yd) / k, _lookupstep(xd) / k
    crs = Rasters.crs(target)
    start = DimensionalData.Lookups.Intervals(DimensionalData.Lookups.Start())
    fwd = DimensionalData.Lookups.ForwardOrdered()
    yvals = collect(range(first(ys), step = sy, length = k * length(ys)))
    xvals = collect(range(first(xs), step = sx, length = k * length(xs)))
    fy = Y(Projected(yvals, sampling = start, crs = crs, order = fwd,
                     span = DimensionalData.Lookups.Regular(sy)))
    fx = X(Projected(xvals, sampling = start, crs = crs, order = fwd,
                     span = DimensionalData.Lookups.Regular(sx)))
    return Rasters.Raster(zeros(length(yvals), length(xvals)), (fy, fx))
end

# A dimension's cell step as its lookup declares it - the `Regular` span, exact from the geotransform
# or the grid's construction - and only where no span is declared the difference of its first two
# coordinates. The step test needs the exact value: the rounding in a difference of two coordinates
# near 60° exceeds the arcsecond tolerance, and a ratio of exactly 1 must read as 1.
function _lookupstep(d)
    return _spanstep(DimensionalData.Lookups.span(DimensionalData.lookup(d)), d)
end

function _spanstep(span::DimensionalData.Lookups.Regular, d)
    return abs(DimensionalData.Lookups.val(span))
end

function _spanstep(::Any, d)
    v = parent(DimensionalData.lookup(d))
    return abs(v[2] - v[1])
end

# Put `raster` on `target` through `_regrid`, erroring clearly if it has no coverage of `target` at
# all rather than silently building an all-inactive environment.
#
# Deliberately says nothing about a change of resolution: what the target grid costs each layer -
# kept exactly, aggregated by a whole factor, or regridded from the covering cells, and why - is
# classified per layer by `_analyse` and reported by the `StudyArea` that decided the grid.
function _sampledata(raster::ClimateRaster, target; name = "raster",
                     categorical::Bool, fn = nothing, centres::Bool = false)
    # The reducer: the caller's, else the axis's - the majority class for class codes, the mean of the
    # cells carrying data for values - the same rule `_reducer` applies on the read path. A **mask**
    # is class codes: `true` and `false` must never be blended.
    reducer = something(fn, categorical ? _majorityclass : _meanpresent)
    # `_regrid` returns a real `(Y, X[, extra])` `DimArray` on the target's own dims - kept as-is (not
    # stripped to a bare `Array`) so the regime/supply built from it carries real CRS provenance, per
    # the `(y, x)` order used throughout.
    out = _regrid(raster, target, reducer, centres = centres)
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
    ys,
    xs = parent(DimensionalData.lookup(yx[1])),
         parent(DimensionalData.lookup(yx[2]))
    dy, dx = _axisstep(yx[1]), _axisstep(yx[2])
    ny, nx = length(ys), length(xs)
    ey,
    ex = Matrix{typeof(1.0°)}(undef, ny, nx),
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
# or a `ConstructedRasterSpec`'s data read happening twice - once to learn the extent and again to rasterise.
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

# A shape - a vector file, a named region, a combination - either as its outline or as the box
# around it.
#
# `outline = false` returns the extent with **no** payload, which is how `_preparemask` already says
# "restrict the grid to this box, but leave every cell in it active" - the same answer an
# `Extents.Extent` gives. The geometries are still prepared to get there, which is a little wasted
# work once per build against a read or download that dominates it.
function _preparemask(active::AbstractShapeSpec, tcrs)
    geoms, extent = _shapegeoms(active, tcrs)
    # No geometry means no ground, which is never a usable mask and is usually a coverage that
    # filtered everything out, a combination whose members do not meet, or a file with no polygon
    # in it. Saying so beats handing back a grid with no active cells, or an extent at the origin.
    isempty(geoms) &&
        error("`$active` selects no ground. A `LandmassesAbove` threshold may have excluded every " *
              "component, a file may hold no polygon, or the members of a combination may not " *
              "overlap - Natural Earth's physical outlines are drawn per landmass, so a " *
              "continent's polygon does not contain its offshore islands.")
    return (payload = active.outline ? geoms : nothing, extent = extent)
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

# Materialising a `ConstructedRasterSpec` mask yields a `Bool` array on its own real grid, so its extent comes
# free from that array's dims - and materialising here rather than later avoids reading the (possibly
# global) source data a second time.
function _preparemask(active::ConstructedRasterSpec, tcrs)
    # Unwrapped **here and only here**: the combine's contract is that it hands back a raster, and
    # everything downstream of this point works in plain arrays. The package owns the wrapper, so
    # reaching through it internally is free; what the contract buys is that *user* code never has to.
    mask = _materialiseconstructed(active).array
    return (payload = mask, extent = _dimsextent(mask, tcrs))
end

function _preparemask(active, tcrs)
    return error("unrecognised `within` argument of type $(typeof(active)); use nothing, a " *
                 "Matrix{Bool}, a LatLong box, or a mask spec (CircleMaskSpec/ShapeSpec/ConstructedRasterSpec).")
end

# A synthetic unitless target `Rasters.Raster` in `crs`, covering the unitful bounds
# `(ylo..yhi, xlo..xhi)` (given in `crs`'s own coordinate unit) in square cells of side `cellside`.
# This is the `size =` override's grid - a uniform step the reference's own grid may not have.
# `cellside` is already the kind of quantity the target is laid out in - a length on a projected
# `crs`, an angle on a geographic one, `_targetcrs` having refused the other pairing - so the only
# conversion is within that kind (`km` to `m`, `arcminute` to `°`). Nothing here converts degrees
# to kilometres.
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
# the *projected* grid, not the degree one); else the reference's own CRS. A `size` must match the
# kind of grid: a length needs a projected target - square metric cells do not exist on a degree
# grid - so if the resolved CRS is geographic we **fail closed** rather than reviving the
# 111.32 km/° approximation, and name a concrete CRS (`_crsadvice`) in the error so the fix is one
# paste away; an angle (`30arcminute`) needs a geographic target, and is refused on a projected one
# for the mirror-image reason.
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
    if !isnothing(size) && _isangle(size) && _isprojectedcrs(resolved)
        error("`cellsize = $size` is an angle, but the target grid is projected, where a cell's " *
              "side is a length. Pass a length (`cellsize = 1km`), or a geographic `crs` such " *
              "as `EPSG(4326)` if a degree grid is what is wanted.")
    elseif !isnothing(size) && !_isangle(size) && !_isprojectedcrs(resolved)
        here = _extentof(_extrema2(_latvals(first(regimes)),
                                   _longvals(first(regimes)))...)
        error("`cellsize = $size` asks for grid cells of a fixed physical side, but the target grid " *
              "is geographic (° coordinates), where a cell's physical size varies with latitude. " *
              "Pass a projected `crs` to build a genuinely metric grid - " *
              "$(_crsadvice(here)) - or give the cell size as an angle " *
              "(`cellsize = 30arcminute`), or omit `cellsize` to keep the data's own native " *
              "resolution.")
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
    return _stepacross(_rastercrs(raster), _latvals(raster), _longvals(raster),
                       tcrs)
end

# The same from a CRS and its two coordinate vectors, so a file opened lazily can be measured
# before it is read.
function _stepacross(crs, lat, long, tcrs)
    (length(lat) > 1 && length(long) > 1) || return nothing
    dlat, dlong = lat[2] - lat[1], long[2] - long[1]
    i, j = cld(length(lat), 2), cld(length(long), 2)
    y0, x0 = lat[i], long[j]
    y1, x1 = y0 + dlat, x0 + dlong
    cell = _bboxin(crs, tcrs,
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
function _asregime(values, categorical::Bool, axis::Type{<:NicheAxis},
                   series::NamedTuple)
    categorical &&
        return _reaxis(CategoricalRegime(values, NoLayerChange()), axis)
    ndims(values) == 2 &&
        return _reaxis(ContinuousRegime(values, NoLayerChange()), axis)
    return _reaxis(_setseries!(ContinuousRegime(_firstslice(values),
                                                NoLayerChange()), values;
                               series...), axis)
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
_specaxis(::RasterSpec{A}) where {A} = A

# The reducer a spec asks for, or `nothing` where its axis decides; no other spec has one. What the
# read-time `scale` uses is what the study grid uses, so a caller's choice holds at both coarsening
# sites.
_specfn(spec::RasterSpec) = spec.fn

_specfn(::Any) = nothing

_specaxis(spec::Tuple) = _sourcepairnotaspec(spec)

# A raster read from a source knows its layer, so its axis is a catalogue lookup like any other -
# only a raster built by hand (no `code`) genuinely has none, and that is what `NicheAxis` records.
function _specaxis(raster::ClimateRaster{S}) where {S}
    (isnothing(raster.code) || raster.code isa AbstractVector) &&
        return NicheAxis
    return something(layeraxis(S, raster.code), NicheAxis)
end

# A `ConstructedRasterSpec` carries the niche axis declared at construction (`NicheAxis` by default).
_specaxis(spec::ConstructedRasterSpec) = spec.axis

# Wrap a sampled supply layer as a supply: `cancel` converts the raw per-area rate (at any native
# time unit) to an absolute per-cell one against `cellarea`, stated in the axis's canonical unit, and
# the **axis** picks the supply type - never the value's dimension, on either count. A
# monthly (3-D) stack becomes a supply holding one slice at a time, carrying the stack as its
# change - exactly as `_asregime` does for a regime.
#
# This is also why there is no separate time-varying supply **type**: such a family would need one
# member per resource, and would refuse any resource that had not been given one. A supply that
# varies in time is the same type as one that does not, so every resource can have one.
function _wrapsupply(out, cellareas, axis, series::NamedTuple)
    T = _supplytype(axis)
    abs = cancel.(out, cellareas, axis)
    ndims(abs) == 2 && return T(abs)
    return _setseries!(T(_firstslice(abs)), abs; series...)
end

# Resolve a `ShapeSpec.path` to a local filesystem path: an already-local path passes
# through; a `CachedAsset` (a URL) is downloaded into its cache (if not already there) here, at
# materialisation time, not at `ShapeSpec` construction.
_resolvepath(path::AbstractString) = path

_resolvepath(asset::CachedAsset) = assetpath(asset)

_resolvepath(request::CDSRequest) = assetpath(request)

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
const _ShapeComponent = @NamedTuple{geometry::ArchGDAL.IGeometry,
                                    envelope::ArchGDAL.GDAL.OGREnvelope,
                                    area::typeof(1.0km^2)}

# Every feature of `level`'s file, grouped by the value of its naming attribute and keyed by the
# lowercased form of it. Each group keeps the source's own spelling of the name alongside.
#
# This is the ONE place a name becomes geometry: `_selectfeatures` is a lookup into it and
# `_levelvalues` reads its keys, so the shipped table's generator and a built shape cannot disagree
# about which features a name covers.
#
# Grouped rather than filtered per name because selection is a linear scan - a shapefile carries no
# attribute index - so answering every name separately costs the scan once per name. Measured on the
# countries file: the scan is 33 ms of a 58 ms single selection, against 5.6 ms to open the zip, so
# grouping is where the generator's time is, not caching the open.
const _REGION_GROUPS = Dict{String,
                            Dict{String,
                                 @NamedTuple{name::String,
                                             geometries::Vector{ArchGDAL.IGeometry}}}}()

# Memoised per level. One level's geometries are held for the process's life, which is a few tens of
# megabytes for the largest file; a session touches one or two levels, and the shipped region table
# answers a bounding box without coming here at all.
function _groupfeatures(level::NaturalEarthLevel)
    return get!(_REGION_GROUPS, level.name) do
        dataset = ArchGDAL.read("/vsizip/" * assetpath(_nesource(level)))
        lyr = ArchGDAL.getlayer(dataset, 0)
        groups = Dict{String,
                      @NamedTuple{name::String,
                                  geometries::Vector{ArchGDAL.IGeometry}}}()
        for feature in lyr
            isnothing(level.within) ||
                _fieldmatches(feature, level.within.first,
                              lowercase(level.within.second)) || continue
            name = _fieldvalue(feature, level.field)
            # "-99" is Natural Earth's unset marker. Left in, it would become a region of its own
            # holding every unassigned feature - a name that looks real and spans the globe.
            (isnothing(name) || isempty(name) || name == "-99") && continue
            group = get!(groups, lowercase(name)) do
                return (name = name, geometries = ArchGDAL.IGeometry[])
            end
            push!(group.geometries, ArchGDAL.clone(ArchGDAL.getgeom(feature)))
        end
        return groups
    end
end

# The geometries `value` names at `level`, empty where the name does not exist there - a name absent
# from a level is a legitimate answer, not an error: the United Kingdom has no map unit of its own.
#
# Matching is case-insensitive because the physical file's own names are inconsistently cased -
# `ALLEGHENY PLATEAU` sits beside `Adelie Coast` - so requiring the source's spelling would make some
# names unreachable in practice.
function _selectfeatures(level::NaturalEarthLevel, value::AbstractString)
    group = get(_groupfeatures(level), lowercase(value), nothing)
    return isnothing(group) ? ArchGDAL.IGeometry[] : group.geometries
end

# Every name defined at `level`, in the source's own spelling, sorted for a stable table.
function _levelvalues(level::NaturalEarthLevel)
    return sort!([g.name for g in values(_groupfeatures(level))])
end

# The string in `feature`'s `field`, or `nothing` where the field is absent from the layer or the
# cell is null. A missing field gives an index of -1 rather than raising.
function _fieldvalue(feature, field::AbstractString)
    i = ArchGDAL.findfieldindex(feature, field)
    (isnothing(i) || i < 0) && return nothing
    value = ArchGDAL.getfield(feature, i)
    return isnothing(value) ? nothing : string(value)
end

# Whether `feature`'s `field` holds `wanted`, which must already be lowercased.
function _fieldmatches(feature, field::AbstractString, wanted::AbstractString)
    value = _fieldvalue(feature, field)
    return !isnothing(value) && lowercase(value) == wanted
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
    # An empty geometry is not a component: a set operation that found no common ground returns
    # one, and left in it would become a part of zero area whose envelope is the origin - a mask of
    # nothing, reported as if it were somewhere.
    present = filter(g -> !ArchGDAL.isempty(g), geoms)
    isempty(present) && return _ShapeComponent[]
    merged = reduce((a, b) -> ArchGDAL.union(_repairgeom(a), _repairgeom(b)),
                    present)
    parts = map(filter(g -> !ArchGDAL.isempty(g), _components(merged))) do g
        return _ShapeComponent((g, ArchGDAL.envelope(g), _equalarea(g)))
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

# The bounding box of a set of components, and whether it crosses the antimeridian.
#
# Natural Earth splits its polygons at the date line, so a region reaching across it arrives as
# components either side and a naive smallest-to-largest longitude spans the globe: the United
# States would read -179.14 to 179.78, which is true and useless.
#
# The east-west extent is therefore taken as the complement of the *widest* longitude gap. Where that
# gap is the one running the long way round outside the data, the box is ordinary; where the widest
# gap lies inside it, the box wraps and `west > east`, which is how RFC 7946 writes one. Latitude
# needs none of this, having no seam.
function _regionbox(parts)
    isempty(parts) &&
        return (west = nothing, south = nothing, east = nothing,
                north = nothing,
                wraps = false)
    south = minimum(p -> p.envelope.MinY, parts)
    north = maximum(p -> p.envelope.MaxY, parts)
    spans = sort!([(p.envelope.MinX, p.envelope.MaxX) for p in parts],
                  by = first)
    merged = Tuple{Float64, Float64}[]
    for (lo, hi) in spans
        if !isempty(merged) && lo <= merged[end][2]
            merged[end] = (merged[end][1], max(merged[end][2], hi))
        else
            push!(merged, (lo, hi))
        end
    end
    lo, hi = merged[1][1], merged[end][2]
    # The gap outside the data, running the long way round through the date line.
    widest, atgap = 360.0 - (hi - lo), 0
    for i in 1:(length(merged) - 1)
        gap = merged[i + 1][1] - merged[i][2]
        gap > widest && ((widest, atgap) = (gap, i))
    end
    atgap == 0 &&
        return (west = lo, south = south, east = hi, north = north,
                wraps = false)
    return (west = merged[atgap + 1][1], south = south, east = merged[atgap][2],
            north = north, wraps = true)
end

# The components a coverage asks for. `parts` must already be ordered largest first, as `_dissolve`
# returns them.
_coverageof(parts::AbstractVector{_ShapeComponent}, ::AllTerritories) = parts

function _coverageof(parts::AbstractVector{_ShapeComponent}, c::LargestLandmass)
    return parts[1:min(c.count, length(parts))]
end

# Components are ordered largest first, so this is a prefix too - but expressed as a threshold, which
# is what "everything except the specks" needs when the count is not known in advance.
function _coverageof(parts::AbstractVector{_ShapeComponent}, c::LandmassesAbove)
    isempty(parts) && return parts
    # A share is of the region's own total, so the threshold is only an area once the parts are in.
    threshold = _thresholdarea(c, sum(p -> p.area, parts))
    return filter(p -> p.area >= threshold, parts)
end

# One feature of a vector file as `_shape` uses it: the prepared geometry to test cells against, and
# the envelope that says which cells those are. The element type is written out because a layer mixes
# `wkbPolygon` and `wkbMultiPolygon` features, so an inferred one keeps only the field names and
# every access through it becomes a dynamic lookup.
const _ShapePart = @NamedTuple{prepared::ArchGDAL.IPreparedGeometry,
                               envelope::ArchGDAL.GDAL.OGREnvelope}

# Reproject `geoms` from `src` into the target grid's own CRS and prepare each for the per-cell
# containment test, with the extent they jointly cover.
#
# Shared by every route that turns geometry into a mask - a vector file through `ShapeSpec`, a named
# region through `NaturalEarthSpec` - so the two cannot come to differ in how a geometry reaches a
# grid.
#
# Each geometry is CLONED before being transformed, and that is load-bearing. `transform!` rewrites in place, and
# `_groupfeatures` hands out geometries it is memoising: transforming those would silently leave the
# cache holding coordinates in whatever CRS was asked for last, so a second build on a different grid
# would reproject already-reprojected ground.
function _preparegeoms(geoms, src, tcrs)
    dest = _gdalcrs(tcrs)
    u = _crsunit(tcrs)
    ylo, yhi, xlo, xhi = Inf, -Inf, Inf, -Inf
    # The envelope is taken *before* `preparegeom` - `ArchGDAL.envelope` on a prepared geometry
    # segfaults - and kept, both for the overall extent here and for `_shape` to window each
    # geometry onto the grid.
    parts = map(geoms) do geom
        g = ArchGDAL.clone(geom)
        ArchGDAL.createcoordtrans(src, dest) do ct
            return ArchGDAL.transform!(g, ct)
        end
        env = ArchGDAL.envelope(g)
        ylo, yhi = min(ylo, env.MinY), max(yhi, env.MaxY)
        xlo, xhi = min(xlo, env.MinX), max(xhi, env.MaxX)
        return _ShapePart((ArchGDAL.preparegeom(g), env))
    end
    extent = isempty(parts) ? nothing :
             _extentof(ylo * u, yhi * u, xlo * u, xhi * u)
    return parts, extent
end

# Every geometry in `spec`'s vector file, with the CRS they are in - the read `ShapeSpec` defers
# from construction.
function _readshapefile(spec::ShapeSpec)
    path = _resolvepath(spec.path)
    vpath = endswith(path, ".zip") ? "/vsizip/" * path : path
    dataset = ArchGDAL.read(vpath)
    lyr = ArchGDAL.getlayer(dataset, spec.layer)
    sr = ArchGDAL.getspatialref(lyr)
    # A missing `.prj` (no CRS metadata) gives a null spatial ref; assume already WGS84 lat/long.
    src = sr.ptr != C_NULL ? sr : _gdalcrs(Rasters.EPSG(4326))
    return [ArchGDAL.clone(ArchGDAL.getgeom(f)) for f in lyr], src
end

"""
    read(spec::AbstractShapeSpec)

Read the ground a shape spec names - a vector file's polygons, a named region's, or what a
combination of shapes builds - as the connected pieces it is, in WGS84 and before any grid
exists: the eager step the spec defers. Each piece is a named tuple of its `geometry` (an
`ArchGDAL` geometry), its `envelope` (the bounding box GDAL reports for it) and its `area` in
square kilometres, ordered largest first, after the spec's `coverage` has said which to keep.

Building a study area resolves a shape through this, so calling it directly is for inspection -
to see how many pieces a name is, or how large the one a `coverage` would drop is.

# Arguments

  - `spec`: what to read.
"""
Base.read(spec::AbstractShapeSpec) = _shapecomponents(spec)

# The components a shape spec resolves to, in WGS84 and before any grid exists: one method per
# leaf, and what `read` returns.
#
# A vector file's own geometry, dissolved into components in WGS84 so that it can take part in a
# combination on equal terms with a named region, then filtered by the spec's coverage as a named
# region's is.
function _shapecomponents(spec::ShapeSpec)
    geoms, src = _readshapefile(spec)
    wgs = _gdalcrs(Rasters.EPSG(4326))
    ArchGDAL.createcoordtrans(src, wgs) do ct
        return foreach(g -> ArchGDAL.transform!(g, ct), geoms)
    end
    return _coverageof(_dissolve(geoms), spec.coverage)
end

# For a single name these are the same calls the shipped table's generator makes, which is what
# keeps a built shape agreeing with the box `boundingbox` reports for it.

function _shapecomponents(spec::NaturalEarthSpec)
    return _coverageof(_dissolve(_selectfeatures(_findlevel(spec.level),
                                                 spec.name)),
                       spec.coverage)
end

# A combination resolves its members first, applies the set operation to their geometry, and only
# then splits the result into components - so the coverage acts on what was built rather than on any
# member. That order is what makes "the British Isles, dropping anything under a square kilometre"
# one expression instead of a filter that could only be applied per member.
function _shapecomponents(spec::ConstructedShapeSpec)
    combined = _applyshapeop(spec.operation,
                             map(_shapegeometry, spec.members))
    return _coverageof(_dissolve([combined]), spec.coverage)
end

# One geometry for a whole spec, its components merged back together, for use as a member of a
# combination.
function _shapegeometry(spec::AbstractShapeSpec)
    parts = _shapecomponents(spec)
    isempty(parts) &&
        error("`$spec` selects no geometry, so it cannot take part in a combination.")
    return _mergegeoms([p.geometry for p in parts])
end

function _mergegeoms(gs)
    return reduce((a, b) -> ArchGDAL.union(_repairgeom(a), _repairgeom(b)), gs)
end

# One method per operation rather than a branch, so an unsupported one is a `MethodError` naming it.
_applyshapeop(::ShapeUnion, gs) = _mergegeoms(gs)

function _applyshapeop(::ShapeIntersection, gs)
    return reduce((a,
                   b) -> ArchGDAL.intersection(_repairgeom(a),
                                               _repairgeom(b)),
                  gs)
end

function _applyshapeop(::ShapeDifference, gs)
    return reduce((a, b) -> ArchGDAL.difference(_repairgeom(a), _repairgeom(b)),
                  gs)
end

# The transforms take one shape. `buffer` and `simplify` want a plain number in the geometry's own
# coordinates, which for these is degrees; a length is converted at the equator, where a degree of
# longitude is widest, so the buffer is never narrower than asked for.
function _applyshapeop(o::ShapeBuffer, gs)
    return ArchGDAL.buffer(only(gs), _indegrees(o.distance))
end

function _applyshapeop(o::ShapeSimplify, gs)
    return ArchGDAL.simplify(only(gs), _indegrees(o.tolerance))
end

_applyshapeop(::ShapeConvexHull, gs) = ArchGDAL.convexhull(only(gs))

# A bare function is the escape hatch that mirrors `ConstructedRasterSpec`'s `combine`: it is handed
# one geometry per member and returns one geometry.
_applyshapeop(f::Function, gs) = f(gs...)

# A distance as a number of degrees, which is what GDAL's `buffer` and `simplify` want here: these
# geometries are in WGS84 lat/long, so their coordinates are angles.
#
# A length is divided by the length of a degree at the equator - the widest a degree of longitude
# gets - so a converted distance never under-reaches the one asked for.
_indegrees(d::Real) = float(d)
_indegrees(d::Unitful.DimensionlessQuantity) = ustrip(°, uconvert(°, d))
_indegrees(d::Unitful.Length) = ustrip(NoUnits, d / _degreelength(1°))

# The geometry a shape spec resolves to, prepared in the target grid's own CRS, with the extent it
# covers. Every shape's components are in WGS84 - a file's are reprojected there on read - so that
# is the source; a file already in the target CRS makes the round trip once per build.
function _shapegeoms(spec::AbstractShapeSpec, tcrs)
    return _preparegeoms([p.geometry for p in _shapecomponents(spec)],
                         _gdalcrs(Rasters.EPSG(4326)), tcrs)
end

# Which indices of `axis` fall within `[lo, hi]`. A raster's Y commonly descends, and `searchsorted`
# given the wrong direction returns an empty range rather than failing, so the direction is read off
# the axis; one monotonic in neither direction falls back to all of it.
function _axiswindow(axis, lo, hi)
    issorted(axis) &&
        return searchsortedfirst(axis, lo):searchsortedlast(axis, hi)
    issorted(axis, rev = true) &&
        return searchsortedfirst(axis, hi,
                                 rev = true):searchsortedlast(axis,
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

# Put a bare Bool `DimArray` (a `ConstructedRasterSpec` mask, say) onto `target` through the same
# pipeline as a layer, a target cell being active when the majority of the source cells it covers
# are. Round-tripped through 0.0/1.0 because GDAL has no Bool dtype; a `NaN`, no covering cell
# present, compares as inactive.
function _samplemask(A::DimensionalData.AbstractDimArray, target)
    out = _regrid(ClimateRaster(SyntheticData, Float64.(A)), target,
                  _majorityclass)
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

# `_preparemask(::AbstractShapeSpec, ...)` already read and reprojected the geometries, so the
# payload is the vector of prepared geometries and their envelopes - nothing is re-read here.
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

# Materialise a `ConstructedRasterSpec` to a raster/array: read each source spec (`_asraster`) and apply
# `combine` - nullary when there are no sources (the thunk produces the layer directly).
function _materialiseconstructed(spec::ConstructedRasterSpec)
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
function _combined(out, spec::ConstructedRasterSpec)
    out isa ClimateRaster ||
        error("a `ConstructedRasterSpec` combine must return a `ClimateRaster`, but this one returned a " *
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

"""
    read(spec::RasterSpec; cut = spec.cut, scale = spec.scale, fn = spec.fn)

Read the data a [`RasterSpec`](@ref) names into a [`ClimateRaster`](@ref) on the source's own grid,
values in the spec's unit and the layer's `code` attached - the eager step the spec defers. A
source that fetches its own files does so here (a catalogued dataset's download), and a spec
naming its files reads them: each opened by the backend its source's catalogue row names, stacked
along time with the coordinates the files carry, or `times` where given, or monthly ordinals where
they have none; a file in the 0 to 360 longitude convention is rolled onto -180 to 180. What the
files say about themselves is checked against the row on the way.

Building a layer reads a spec through this, windowed and coarsened to suit the study grid, so
calling it directly is for inspection - a whole global layer at its own resolution can be large.

# Arguments

  - `spec`: what to read.
  - `cut`: an `Extents.Extent` of `°` intervals to window the read to; the spec's own by default.
  - `scale`: an integer factor to coarsen by on read, each block of `scale × scale` cells becoming
    one; the spec's own by default, and `1` where it states none.
  - `fn`: how a block is reduced to one cell; the spec's own by default, and the axis's choice
    where it states none.
"""
function Base.read(spec::RasterSpec; cut = spec.cut, scale = spec.scale,
                   fn = spec.fn)
    return _readspec(spec, spec.files, cut = cut, scale = something(scale, 1),
                     fn = fn)
end

"""
    fetchfiles(spec::RasterSpec; dryrun = false)
    fetchfiles(specs; dryrun = false)

Fetch every file a spec reads that is not there yet, reading nothing, and return their local paths
in the spec's order - the step to run on a node with a network before a run on nodes without one.
A vector of specs is fetched spec by spec.

# Arguments

  - `spec`, `specs`: what to fetch for.
  - `dryrun`: `true` fetches nothing and instead returns, per file, a named tuple of the `entry`,
    whether it is `present`, and its `bytes` where the server states them (a plain download's
    `Content-Length`; a Climate Data Store request cannot say), so a user on a slow link can decide.
    A source that resolves its own files through RasterDataSources is reported as one entry with
    no size.
"""
function fetchfiles(spec::RasterSpec; dryrun::Bool = false)
    dryrun || return _fetchedpaths(spec, spec.files)
    isnothing(spec.files) &&
        return [(entry = spec.source, present = nothing, bytes = nothing)]
    return [_dryrunentry(entry) for entry in spec.files]
end

function fetchfiles(specs::AbstractVector; dryrun::Bool = false)
    return reduce(vcat, (fetchfiles(s, dryrun = dryrun) for s in specs),
                  init = dryrun ? Any[] : String[])
end

# The local paths of a spec's files, fetched where missing: named entries through `assetpath`, a
# source's own through its fetch hook.
function _fetchedpaths(spec::RasterSpec, files::AbstractVector)
    return String.(_resolvepath.(files))
end

function _fetchedpaths(spec::RasterSpec, ::Nothing)
    readkw = (; _getrasterkw(spec.source)..., spec.readkw...)
    return _allfiles(_fetchfiles(spec.source, spec.code; readkw...))
end

# Every path in whatever shape a fetch hook returned: one, a vector, or named layers per time.
function _allfiles(raw::Vector{<:NamedTuple})
    return String[String(p) for nt in raw for p in values(nt)]
end

_allfiles(raw) = _filelist(raw)

# One file's dry-run report: where it would land, whether it is there, and what the server says
# it weighs.
function _dryrunentry(entry)
    path = _localpath(entry)
    return (entry = entry, present = isfile(path), bytes = _remotesize(entry))
end

# Where an entry lives, or would, without fetching it.
_localpath(path::AbstractString) = String(path)

function _localpath(asset::CachedAsset)
    return something(asset.path,
                     joinpath(assetdir(owner = asset.owner),
                              basename(asset.url)))
end

_localpath(request::CDSRequest) = request.path

# The size a server states for a download, from a `HEAD`, or `nothing` where it states none or
# the entry is not a plain download.
function _remotesize(asset::CachedAsset)
    response = try
        Downloads.request(asset.url, method = "HEAD",
                          downloader = _http11downloader(), throw = false)
    catch
        return nothing
    end
    i = findfirst(h -> lowercase(first(h)) == "content-length",
                  response.headers)
    return isnothing(i) ? nothing :
           tryparse(Int, String(last(response.headers[i])))
end

_remotesize(::Any) = nothing

function provenance(spec::RasterSpec)
    entries = isnothing(spec.files) ? _presentfiles(spec) : spec.files
    return Union{Nothing, InputRecord}[provenance(_localpath(e))
                                       for e in entries]
end

"""
    verifyassets(spec::RasterSpec)

Check every file a spec reads that is present against the checksum its provenance record holds,
erroring on the first that differs - a truncated or replaced copy - and return how many were
checked. A file with no record, or whose record carries no checksum, has nothing to check
against and is not counted; nothing is fetched.

# Arguments

  - `spec`: the spec whose files to check.
"""
function verifyassets(spec::RasterSpec)
    entries = isnothing(spec.files) ? _presentfiles(spec) : spec.files
    checked = 0
    for path in _localpath.(entries)
        isfile(path) || continue
        checked += _verifyfile(path)
    end
    return checked
end

# The files a source has on disk for a spec that resolves its own, fetching nothing.
function _presentfiles(spec::RasterSpec)
    readkw = (; _getrasterkw(spec.source)..., spec.readkw...)
    return _localfiles(spec.source, spec.code; readkw...)
end

# A spec naming its files reads them itself: each through `_cachedlayer` - the same step a dataset
# layer takes, so a coarsened read of a whole file is memoised on disk exactly as a dataset's is -
# with the backend and the variable name the source's row and code decide, its magnitudes expressed
# in the layer table's unit where the file states its own, stacked in time, rolled and checked
# against the row. A bare file - no code - consults no row: its `source` is provenance it records,
# not a dataset whose layers it claims to be.
function _readspec(spec::RasterSpec{A}, files::AbstractVector; cut, scale,
                   fn) where {A}
    spec.code isa AbstractVector &&
        error("a spec naming its files reads one layer: name the layer you want, or pass a codes " *
              "vector to `ConstructedRasterSpec`, which reads each on its own terms.")
    rec = _specrecord(spec)
    open = _openkw(spec)
    paths = _resolvepath.(files)
    isnothing(rec) || _checkcatalogue(rec, _lazyopen(first(paths); open...))
    layers = map(paths) do p
        # The `;` is load-bearing: `open` splats as keywords only after it.
        return _cachedlayer(p, scale, fn, NoUnits; cut = cut, axis = A,
                            expressedin = _tableunit(spec), open...)
    end
    world = _stackfiles(layers, spec.times)
    if !isnothing(rec) && _needswrap(rec.longituderange, world)
        world = _wraplong180(world)
    end
    world = _applycut(world, cut)
    data = isnothing(spec.code) ? world :
           _applyfold(_applyperiod(world, spec, spec.readkw), spec)
    return _attachunit(ClimateRaster(spec.source, data, spec.code), spec.unit)
end

# How a spec's files are opened: the backend its source's catalogue row names (`nothing` sniffs
# the filename, for a file that belongs to no dataset), and, for a netCDF file, the variable to
# take - the spec's code - and the level to select on a file holding several, the top of the
# layer the code's row spans.
function _openkw(spec::RasterSpec)
    rec = _specrecord(spec)
    isnothing(rec) && return (source = nothing, name = nothing, level = nothing)
    netcdf = rec.format === :netCDF && spec.code isa CODE_TYPE
    name = netcdf ? Symbol(spec.code) : nothing
    level = netcdf ? _layertop(layerinfo(spec.source, spec.code)) : nothing
    return (source = _rasterssource(rec.format), name = name, level = level)
end

# The depth of the top of the layer a row's `VerticalExtent` range spans, as a file states a soil
# level - positive downward - or `nothing` for a row naming a height or no extent.
function _layertop(rec::LayerRecord)
    rec.verticalextent isa Tuple || return nothing
    return -max(rec.verticalextent...)
end

# The `datasets.csv` row a spec reads under: its source's, for a spec naming a layer of it, and
# none for a bare file.
function _specrecord(spec::RasterSpec)
    return isnothing(spec.code) ? nothing :
           _datasetrecord(spec.source)
end

# The unit a catalogued layer's table declares - the amount, before any accumulation period turns
# it into a rate - which is what a file stating its own unit is converted into on read, so that
# `_applyperiod` and the spec's rate unit then apply exactly as they do to a fetched layer. A bare
# file has no table, and its magnitudes are taken in the spec's unit as they are.
function _tableunit(spec::RasterSpec)
    spec.code isa CODE_TYPE || return spec.unit
    return layerunit(spec.source, spec.code)
end

# A catalogued spec whose source resolves its own files: fetched through the source's hook with
# its own keywords, read in whatever shape they arrive, and the published-scale correction, the
# accumulation period and the unit applied.
function _readspec(spec::RasterSpec{A}, ::Nothing; cut, scale, fn) where {A}
    S = spec.source
    # The spec's own axis decides the aggregation reducer unless the read options say otherwise -
    # it is the one statement of what the layer holds, whether the catalogue or the caller made it.
    readkw = (; _getrasterkw(S)..., spec.readkw...)
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
    if spec.code isa AbstractVector
        onelayer = "Name the layers you want instead - `SourceSpec($(spec.source), :code)` for " *
                   "one, or pass a codes vector to `ConstructedRasterSpec`, which reads each on its " *
                   "own terms."
        us = unique(layerunit(spec.source, c) for c in spec.code)
        length(us) == 1 ||
            error("`SourceSpec($(spec.source))` covers $(length(spec.code)) layers with " *
                  "$(length(us)) different units ($(join(us, ", "))), so they cannot be read " *
                  "into one array without losing them. " * onelayer)
        # Called for its refusal: the vector method throws on a stack that mixes class codes with
        # measurements, which is the rule stated once rather than restated here.
        iscategorical(spec.source, spec.code)
    end
    raw = _readraw(S, _fetchfiles(S, spec.code; readkw...), cut = cut,
                   scale = scale, fn = fn,
                   slices = get(readkw, :month, nothing), axis = A)
    corrected = _rescalepublished(S, spec.code, raw)
    # Tag the materialised raster with the layer it came from. This is the single point at which a
    # spec becomes data, and the only place that still knows the code - after this the raster is
    # passed around on its own and the shipped table can no longer be consulted for it, which is
    # exactly why `iscategorical` answers `false` for a raster carrying no code.
    data = _applyfold(_applyperiod(corrected.array, spec, readkw), spec)
    return _attachunit(ClimateRaster(S, data, spec.code), spec.unit)
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
function _applyperiod(array, sl::RasterSpec, readkw::NamedTuple)
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

# Multiply a freshly-read array by the factor its layer's reading on its axis carries (`_fold`): a
# volumetric fraction by the thickness of the layer it was measured over, a mass of water per
# area by the volume that mass of water fills, so the values are a depth over the cell before the
# unit `_foldedunit` gave the spec is attached. Most layers carry a factor of one and are returned
# untouched.
function _applyfold(array, sl::RasterSpec)
    sl.code isa AbstractVector && return array
    factor = _foldfactor(layerinfo(sl.source, sl.code))
    factor == 1 && return array
    return array .* factor
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
                 "`ConstructedRasterSpec(() -> raster; axis = SomeAxis)`.")
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
# single-layer source builder. (`ConstructedRasterSpec` layers are always `SourceSpec`s - see
# `_parselayers`; a bare-dataset layer is a whole-dataset `SourceSpec`, read via `_read`.)
_asraster(raster::ClimateRaster) = _rasternotaspec(raster)

_asraster(spec::RasterSpec) = read(spec)

_asraster(spec::Tuple) = _sourcepairnotaspec(spec)

# ---------------------------------------------------------------------------
# Single-signature GridHabitat (the public spec-based API)
# ---------------------------------------------------------------------------

# A `ConstructedRasterSpec` as a regime/supply layer: materialise it to a raster (its `combine` result).
_asraster(spec::ConstructedRasterSpec) = _materialiseconstructed(spec)

# A shape has no grid, so it cannot become a raster on its own account: a resolution would have to be
# invented before the study grid had been decided. This is the wrong question rather than a missing
# method, so it says which question was meant instead of failing as a `MethodError` on a private
# function.
function _asraster(spec::AbstractShapeSpec)
    return error("`$(nameof(typeof(spec)))` is geometry, not a raster, so it cannot be a member " *
                 "of a `ConstructedRasterSpec`: it has no grid of its own to be read onto, and " *
                 "choosing one here would fix a resolution before the study grid was decided. " *
                 "To combine it with other geometry use `ConstructedShapeSpec`, which composes " *
                 "shapes exactly and at no resolution; to use it as a mask on a grid pass it as " *
                 "`within` to `StudyArea`.")
end

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
