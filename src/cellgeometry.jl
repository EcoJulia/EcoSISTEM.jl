# SPDX-License-Identifier: LGPL-3.0-or-later

using FillArrays: Fill

# --- Cell geometry: one question, asked of anything that knows the grid -----
#
# **Why these exist.** A cell's size was previously stored on each layer as a scalar `size` field,
# supplied at construction and therefore able to disagree with the coordinates beside it — which is
# exactly what happened: every supply reported `1.0 m` whatever its grid. These read the grid itself,
# so there is one answer and one place it comes from.
#
# **Setup-only, deliberately.** `dynamics.jl` never calls these: dispersal lookups are built once
# in the `Ecosystem` constructor and the hot loop reads the cached table. So a `Union{T, Nothing}`
# return costs nothing here, and these must not be called per-cell.

# --- The public grid questions ------------------------------------------------------------

"""
    getgridshape(x)

Return the size of `x`'s grid in cells, as `(ny, nx)`, or `nothing` where `x` has no grid.

This is how big the grid *is* — how many rows and how many columns — and is what the dispersal and
intervention code needs to turn a flat location index into a `(y, x)` position.

  - `x` — anything that knows the grid, as [`getcellareas`](@ref): a [`StudyArea`](@ref), its report,
    a [`ClimateRaster`](@ref), a raster, a layer or collection, a habitat, or an ecosystem.
  - `force` — as [`getcellareas`](@ref): read a data-backed spec's data rather than answering
    `missing`.

`shape` rather than `size` is doing real work in the name. `getcellareas` and [`getgridarea`](@ref)
are the same quantity at two scales and read as a pair, but a cell's *size* is a length while a
grid's is a count, so those two should not. `getgridsize` keeps its released meaning — one cell's
side length — and is deprecated onto it.

Deliberately a plain `Tuple` rather than a [`SpatialSize`](@ref): these are counts, with no units and
no frame, so nothing about them is spatial in the way a distance is.

One method covers everything, because [`getcellareas`](@ref) and the rest route through the same
`_gridyx` dispatch table — so a habitat answers directly and no caller needs to reach inside it for a
regime.

See also [`getcellcount`](@ref), which is `prod` of this and can exclude inactive cells.
"""
function getgridshape(x; force = false)
    yx = _resolvegrid(x, force)
    (isnothing(yx) || ismissing(yx)) && return yx
    return (length(yx[1]), length(yx[2]))
end

# --- The field-returning accessors -------------------------------------------------------------
#
# Every one of these answers "one value per cell", grid-shaped, so a caller can index it, take its
# `size`, broadcast over it or reduce it without ever asking whether the quantity happens to be
# uniform. Where it is, the answer is a `Fill` and costs 24 bytes however large the grid; where it is
# not, it is materialised. The alternative — an `(ny, 1)` column that broadcasts but cannot be indexed
# `[i, j]` — would make the shape of the answer depend on the data, which is what these exist to stop.
#
# Reductions are Base's. There is deliberately no `mean` keyword: `mean(getcellareas(x))` says it
# better, and `sum` over a mask is how `getgridarea` is built.
"""
    getcellareas([units,] x)

Return the area of every cell of `x`'s grid, one value per cell.

The result is always grid-shaped, so it can be indexed, reduced or broadcast without asking whether
the grid happens to be uniform. Where every cell has the same area the answer is a `Fill`, which
costs the same however large the grid; where areas vary -- as they do with latitude on a geographic
grid -- it is materialised.

  - `units` — the unit to answer in, and **optional**: a length-squared unit (`km^2`) asks for a
    physical area, `°^2` or `sr` for a true solid angle, and omitting it gives the grid's own.
  - `x` — anything that knows the grid: a [`StudyArea`](@ref), its report, a
    [`ClimateRaster`](@ref), a raster, a layer or collection, a habitat, or an ecosystem. Something
    with no grid answers `nothing`.
  - `force` — whether to **read the data** to find out. A data-backed spec has a grid, but only
    reading it would say what, so it answers `missing` rather than `nothing`: the value exists and
    is unknown here. `force = true` performs that read, which may involve a download. Only a
    [`SourceSpec`](@ref) names data that can be read on its own; any other lazy spec describes a
    computation whose grid comes from the area it is built on, and refuses.

Reduce with Base: `mean(getcellareas(x))`, `extrema(...)`, `sum(...)`.

See also [`getcellsizes`](@ref) and [`getgridarea`](@ref).
"""
function getcellareas(units, x; force = false)
    _checkareaunit(units)
    yx = _resolvegrid(x, force)
    (isnothing(yx) || ismissing(yx)) && return yx
    shape = (length(yx[1]), length(yx[2]))
    native = _nativecellsizeyx(yx)
    isnothing(native) && return nothing
    if native.y isa Unitful.Length
        # A projected grid: every cell is the same rectangle, so one value covers the whole grid.
        a = native.y * native.x
        ismissing(units) && return Fill(a, shape...)
        units isa Unitful.AreaUnits || return _noangularfromprojected()
        return Fill(uconvert(units, a), shape...)
    end
    # A geographic grid: the area falls towards the poles, so there is a value per row. Angular
    # units (and no unit at all) want the true solid angle; a length unit wants physical ground.
    if ismissing(units) || units isa Unitful.FreeUnits{<:Any, NoDims}
        col = _angularareacolumn(_gridyx(x))
        out = ismissing(units) ? uconvert.(°^2, col) : uconvert.(units, col)
        return _spreadrows(out, shape)
    end
    # A physical area from a geographic grid: each row's true area on the sphere, spread across the
    # grid's columns. Routed through `_cellareasyx`, which is what the supply constructor uses, so
    # the two cannot disagree about how much ground a cell covers.
    return _spreadrows(uconvert.(units, vec(_cellareasyx(_gridyx(x)))), shape)
end
getcellareas(x; kw...) = getcellareas(missing, x; kw...)

"""
    getcellsizes([units,] x)

Return the extent of every cell of `x`'s grid as a [`SpatialSize`](@ref), one value per cell.

Each component is grid-shaped on the same terms as [`getcellareas`](@ref): a `Fill` where the extent
is uniform, materialised where it is not.

  - `units` — the unit to answer in, and **optional**: a length (`km`, `m`) asks for a metric
    answer, an angle (`°`) for an angular one, and omitting it gives the grid's own. The unit
    chooses the *frame*, not merely a conversion.
  - `x`, `force` — as [`getcellareas`](@ref).

This is not the square root of [`getcellareas`](@ref). On a geographic grid a cell's extent is
angular and constant while its area is physical and shrinks towards the poles.
"""
function getcellsizes(units, x; force = false)
    yx = _resolvegrid(x, force)
    (isnothing(yx) || ismissing(yx)) && return yx
    shape = (length(yx[1]), length(yx[2]))
    native = _nativecellsizeyx(yx)
    isnothing(native) && return nothing
    if ismissing(units)
        return SpatialSize(Fill(native.y, shape...), Fill(native.x, shape...))
    end
    if _samekind(native.y, units)
        return SpatialSize(Fill(uconvert(units, native.y), shape...),
                           Fill(uconvert(units, native.x), shape...))
    end
    # A metric size asked of a geographic grid. The north-south extent of a cell is the same
    # everywhere -- a degree of latitude is a degree of latitude -- while the east-west extent
    # shrinks as the cosine of latitude. So `y` is constant and `x` varies by row.
    #
    # Both are materialised at grid shape even though `y` is uniform, because the two components of a
    # `SpatialSize` share one type. Compressing `y` to a `Fill` and leaving `x` a matrix would not
    # typecheck, and compressing both to their average would claim a uniformity the grid does not
    # have -- which is what the old stored `size` field was removed for doing.
    # Angular asked of a projected grid: the inverse transform, one answer per cell. It varies, and
    # on an equal-area projection it is exactly what varies while the area does not.
    if !(units isa Unitful.LengthFreeUnits) && native.y isa Unitful.Length
        crs = _gridcrs(x)
        isnothing(crs) && return nothing
        e = _angularextents(_gridyx(x), crs)
        return SpatialSize(uconvert.(units, e.y), uconvert.(units, e.x))
    end
    (units isa Unitful.LengthFreeUnits && !(native.y isa Unitful.Length)) ||
        return nothing
    yx = _gridyx(x)
    ns = uconvert(units, _axisstep(yx[1]) * LONGITUDE_DEGREE_LENGTH)
    ew = uconvert.(units,
                   _axisstep(yx[2]) * LONGITUDE_DEGREE_LENGTH .*
                   cos.(_dimcentres(yx[1])))
    return SpatialSize(fill(ns, shape...), _spreadrows(ew, shape))
end
getcellsizes(x; kw...) = getcellsizes(missing, x; kw...)

"""
    getcellat(x, place)

Return the `CartesianIndex` of the cell of `x`'s grid that contains `place`, or `nothing` where `x`
has no grid.

  - `x` — anything that knows the grid, as [`getcellareas`](@ref).
  - `place` — a [`SpatialLocation`](@ref) (or a [`LatLong`](@ref)) in the grid's own frame.

This is the one function that turns a position into a grid position; everything else is then
indexing, which is why no other accessor takes a location as a keyword.

Selection is by the lookup's own `Contains`, so it respects where in each cell the coordinate sits
rather than re-deriving it. A place outside the grid throws: it is a caller mistake rather than a
missing value.
"""
function getcellat(x, place)
    yx = _gridyx(x)
    isnothing(yx) && return nothing
    _checkplaceframe(yx, place)
    return CartesianIndex(_axisindexat(yx[1], _placey(place)),
                          _axisindexat(yx[2], _placex(place)))
end

"""
    getcellarea([units,] x, place)

Return the area of the single cell of `x`'s grid containing `place`.

The plural [`getcellareas`](@ref) answers for every cell; this answers for one, and is exactly
`getcellareas(units, x)[getcellat(x, place)]`.

  - `units` — as [`getcellareas`](@ref).
  - `x`, `place` — as [`getcellat`](@ref).

There is no `force` here: a place has to be stated in the grid's own frame, so a caller who holds one
already knows the grid. Force the read through [`getcellareas`](@ref) if that is what is wanted.
"""
function getcellarea(units, x, place)
    areas = getcellareas(units, x)
    isnothing(areas) && return nothing
    return areas[getcellat(x, place)]
end
getcellarea(x, place) = getcellarea(missing, x, place)

"""
    getcellsize([units,] x, place)

Return the extent of the single cell of `x`'s grid containing `place`, as a
[`SpatialSize`](@ref) of scalars.

The plural [`getcellsizes`](@ref) answers for every cell; this answers for one.

  - `units` — as [`getcellsizes`](@ref).
  - `x`, `place` — as [`getcellat`](@ref).

There is no `force` here, for the reason given on [`getcellarea`](@ref).
"""
function getcellsize(units, x, place)
    sizes = getcellsizes(units, x)
    isnothing(sizes) && return nothing
    i = getcellat(x, place)
    return SpatialSize(sizes.y[i], sizes.x[i])
end
getcellsize(x, place) = getcellsize(missing, x, place)

"""
    getgridarea([units,] x; active = false)

Return the total area of `x`'s grid, or of its active cells only when `active = true`.

This is `sum(getcellareas(x))`, and is correct on every projection because it adds each cell's own
area rather than multiplying one cell by a count — which would be wrong wherever cells differ.

  - `units`, `x`, `force` — as [`getcellareas`](@ref).
  - `active` — count only the cells the study area decided to simulate.
"""
function getgridarea(units, x; active = false, force = false)
    areas = getcellareas(units, x, force = force)
    (isnothing(areas) || ismissing(areas)) && return areas
    active || return sum(areas)
    mask = _activemask(x)
    return isnothing(mask) ? sum(areas) : sum(areas[mask])
end
getgridarea(x; kw...) = getgridarea(missing, x; kw...)

"""
    getcellcount(x; active = false)

Return how many cells `x`'s grid has, or how many are active when `active = true`.

  - `x`, `force` — as [`getcellareas`](@ref).
  - `active` — count only the cells the study area decided to simulate.
"""
function getcellcount(x; active = false, force = false)
    shape = getgridshape(x, force = force)
    (isnothing(shape) || ismissing(shape)) && return shape
    active || return prod(shape)
    mask = _activemask(x)
    return isnothing(mask) ? prod(shape) : count(mask)
end

# --- Private helpers ----------------------------------------------------------------------

# The `(Y, X)` dims of whatever `x` is, or `nothing` when it has no grid at all. One method per kind
# of thing that knows, so "who can answer" is a dispatch table rather than a chain of checks.
_gridyx(::Any) = nothing
_gridyx(area::StudyArea) = _gridyx(area.report)
_gridyx(report::StudyAreaReport) = dims(report.active, (Y, X))
_gridyx(A::DimensionalData.AbstractDimArray) = dims(A, (Y, X))
_gridyx(raster::ClimateRaster) = dims(raster.array, (Y, X))
_gridyx(layer::AbstractLayer) = _yx(layer)
_gridyx(habitat::GridHabitat) = _yx(habitat.regime)
_gridyx(eco::AbstractEcosystem) = _gridyx(eco.habitat)

# A **synthetic** spec has no grid at all — it is a rule, not data, and adopts whatever grid it is
# placed on. `nothing` is the honest answer, and it is the same signal `_rastercellstep` already uses
# for "cannot say", which `_targetcellsize` filters on when deciding a study area's resolution.
_gridyx(::AbstractSyntheticLayerSpec) = nothing
_gridyx(::AbstractSyntheticMaskSpec) = nothing

# A **data-backed** spec has a grid, but only reading it -- potentially a download -- would say what.
# That is `missing`, not `nothing`: the value exists and is unknown here, where `nothing` means there
# is none to know. A caller who wants the read asks for it with `force = true`.
_gridyx(::AbstractLazySpec) = missing

# One axis's step, from the lookup's **declared** `Regular` span where it has one, and only otherwise
# differenced out of the coordinates.
#
# **The declared step is the reliable one**, for the reason `_rastercellstep` records: differencing
# gives a subtly different `Float64` depending on where in the grid you do it — the same WorldClim
# layer read globally and read cut to Scotland differ in the 13th digit. Every grid this package
# builds carries a `Regular` span (and an index-*range* crop preserves it), so the fallback is for
# grids that arrive from elsewhere.
function _axisstep(d)
    sp = DimensionalData.Lookups.span(d)
    sp isa DimensionalData.Lookups.Regular && return abs(sp.step)
    l = parent(DimensionalData.lookup(d))
    return length(l) > 1 ? abs(l[2] - l[1]) : oneunit(eltype(l))
end

# The same, from dims already in hand, so an accessor resolves the grid once rather than twice.
_nativecellsizeyx(::Nothing) = nothing
_nativecellsizeyx(::Missing) = missing
function _nativecellsizeyx(yx)
    all(_isreallookup, yx) || return nothing
    return (y = _axisstep(yx[1]), x = _axisstep(yx[2]))
end

# Is this the grid's own kind of unit? A length answer from a metric grid, or an angular one from a
# geographic grid, needs no conversion beyond `uconvert`; the cross-kind cases need a location.
function _samekind(::Unitful.Quantity{<:Any, D},
                   ::Unitful.FreeUnits{<:Any, D}) where {D}
    return true
end
_samekind(_, _) = false

# The lower and upper edges of cell `i` on axis `d`.
#
# **The locus is load-bearing and must be read, not assumed.** A coordinate labels its cell's
# `Start`, `Center` or `End`, so the same number means three different pieces of ground; a solid angle
# computed from the wrong pair is wrong by up to a whole cell. Every grid this package builds declares
# `Intervals(Start)`, but a raster from elsewhere may not, and assuming would be silent.
function _celledges(d, i)
    v = parent(DimensionalData.lookup(d))[i]
    step = _axisstep(d)
    loc = DimensionalData.Lookups.locus(DimensionalData.sampling(d))
    lo = loc isa DimensionalData.Lookups.Start ? v :
         loc isa DimensionalData.Lookups.Center ? v - step / 2 : v - step
    return (lo, lo + step)
end

# The solid angle of a cell spanning latitudes `φ1…φ2` and `dlong` of longitude, in steradians.
#
# `Δλ · (sin φ₂ − sin φ₁)` — the exact expression, not a small-angle approximation, and the reason
# a cell's angular area **shrinks towards the poles** while its angular *extent* does not.
function _solidangle(φ1, φ2, dlong)
    return ustrip(u"rad", dlong) * (sin(φ2) - sin(φ1)) * u"sr"
end

# Refuse a unit that is not an area, naming the fix.
#
# `getcellareas` insists on a real area unit — `km²`, `m²`, `°²`, `sr` — rather than squaring a
# base unit for you. `getcellareas(km, x)` is a *caller mistake*, not an unanswerable question, so it
# errors where `nothing` would swallow a typo; `nothing` stays reserved for "cannot answer".
# The point is that the first argument means **the same thing in both functions**: the unit the
# answer comes back in. A `km` that meant kilometres to `getcellsizes` and kilometres-squared to
# `getcellareas` is the kind of thing that reads fine and confuses later.
#
# **The angular set still has to be named, and Unitful gives no way round it.** `°`, `°²` and `sr`
# are all `NoDims`, so a solid-angle unit cannot be told from an angle by dimension — measured:
# `uconvert(°, 1.0sr)` cheerfully returns `57.29577951308232°` rather than refusing. `Unitful.AreaUnits`
# catches `km²` because that is `𝐋²`; there is no `SolidAngleUnits` to catch `sr`. So this is the same
# small hardcoded set as before, repurposed from *squaring* to *validating*.
function _isareaunit(u)
    return u isa Unitful.AreaUnits || u == u"sr" || u == °^2 ||
           u == u"rad"^2
end

# Refuse anything that is not an area unit before it reaches the conversion, so the error names the
# argument the caller wrote. A bare `Dimensions` passes: asking for a dimension keeps native units
# rather than converting, which is a different request.
function _checkareaunit(units)
    (ismissing(units) || units isa Unitful.Dimensions || _isareaunit(units)) &&
        return nothing
    return error("`getcellareas` answers in an area unit, and `$units` is not one — pass " *
                 "`$(units)^2` if that is what you meant (`km^2`, `°^2`), or `sr` for a solid " *
                 "angle. `getcellsizes` is the one that takes a length or an angle.")
end

# Resolve `x`'s grid, honouring `force`. Three answers, kept distinct on purpose:
#
# dims     there is the grid
# nothing  there is no grid -- `x` is a rule, not data, and adopts whatever it is placed on
# missing  there IS one, but only reading would say what, and reading was not asked for
#
# Collapsing the last two would make "we did not look" indistinguishable from "there is nothing to
# look at", which is the difference `force` exists to let a caller act on.
_resolvegrid(x, force::Bool) = _gridyx(x)
function _resolvegrid(spec::AbstractLazySpec, force::Bool)
    force || return missing
    return _forcedgrid(spec)
end

# Read a data-backed spec to find its grid. Only a `SourceSpec` names data that can be read on its
# own; the others describe a computation whose grid comes from the area they are built on, so there
# is nothing to read in isolation and the study area is the only thing that can answer.
_forcedgrid(spec::SourceSpec) = _gridyx(_read(spec).array)
function _forcedgrid(spec::AbstractLazySpec)
    return error("`force = true` cannot read a $(nameof(typeof(spec))) on its own: it describes a " *
                 "computation over inputs, and the grid it lands on is decided by the `StudyArea` " *
                 "it is built with. Ask that area instead.")
end

# The one frame conversion still not written: a projected cell's angular AREA -- the solid angle it
# subtends. `getcellsizes` now answers the matching size question by transforming each cell back
# through its CRS (`_angularextents`); the area needs those transformed corners fed through
# `_solidangle` rather than differenced, which is not done.
#
# It WARNS and returns `nothing`, where a genuine "there is no answer" stays silent. The distinction
# is worth the noise: a `nothing` that reflects unwritten code should be distinguishable from one
# that is the answer to the question asked -- as a grid with no CRS at all now gives, silently,
# because for that grid there is genuinely nothing to transform through.
function _noangularfromprojected()
    @warn "`getcellareas` cannot yet give an angular area for a projected grid — that needs each " *
          "cell's transformed corners fed through a solid angle, which is not implemented. " *
          "`getcellsizes` does answer the matching size question. Ask in an area unit, or omit " *
          "the unit for the grid's own. Returning `nothing`."
    return nothing
end

# The CRS a grid's own `Y` dim declares, or `nothing` where it declares none -- a synthetic grid is
# plain cell indices with no real-world position, so there is nothing to transform through.
function _gridcrs(x)
    yx = _resolvegrid(x, false)
    (isnothing(yx) || ismissing(yx)) && return nothing
    return Rasters.crs(yx[1])
end

# Spread a per-row column across the grid's columns. The varying case is materialised rather than
# left as an `(ny, 1)` that only broadcasting can use.
_spreadrows(col, shape) = repeat(reshape(col, :, 1), 1, shape[2])

# One solid angle per row of a geographic grid, from each row's own latitude interval.
function _angularareacolumn(yx)
    dlong = _axisstep(yx[2])
    lats = parent(DimensionalData.lookup(yx[1]))
    return [_solidangle(_celledges(yx[1], i)..., dlong)
            for i in eachindex(lats)]
end

# The two components of a place. One method each: a `LatLong` is a `SpatialLocation` whose type
# parameter is pinned to degrees, so it needs no method of its own.
_placey(p::SpatialLocation) = p.y
_placex(p::SpatialLocation) = p.x

# The index containing `v` on one axis, by the lookup's own `Contains` -- as an `Int` for one value
# and a `Vector{Int}` for several, so a caller's result shape falls out of what it asked for and
# needs no separate flag.
#
# `Contains`, never `Near`: the two agree at cell boundaries and differ off the grid, where `Near`
# silently invents an answer and `Contains` refuses. Refusing is what makes "the value at this place"
# a claim rather than a guess.
function _axisindexat(d, v)
    return DimensionalData.Lookups.selectindices(DimensionalData.lookup(d),
                                                 DimensionalData.Lookups.Contains(v))
end
_axisindexat(d, vs::AbstractVector) = [_axisindexat(d, v) for v in vs]

# Refuse a place stated in a frame the grid is not in, rather than letting `Contains` fail on a unit
# mismatch with a message about lookups. A degree place on a projected grid is the common slip.
function _checkplaceframe(yx, place)
    step = _axisstep(yx[1])
    _samekind(step, Unitful.unit(_placey(place))) ||
        error("this grid's coordinates are in $(Unitful.unit(step)) and the place given is in " *
              "$(Unitful.unit(_placey(place))), so it cannot say which cell is meant. Give the " *
              "place in the grid's own frame.")
    return nothing
end

# The active-cell mask of whatever knows one, or `nothing` when the thing has no notion of inactive
# cells. Only the study area and what is built on it decide activity; a bare raster does not.
# `parent` is used on every branch, and is a no-op on a plain `Matrix`, so a mask arrives as
# something a `Fill` or a materialised field can be indexed by whichever way it was stored.
_activemask(::Any) = nothing
_activemask(area::StudyArea) = _activemask(area.report)
_activemask(report::StudyAreaReport) = parent(report.active)
_activemask(habitat::GridHabitat) = parent(habitat.active)
_activemask(eco::AbstractEcosystem) = _activemask(eco.habitat)
