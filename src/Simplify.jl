# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Distributions
using Random
using AxisArrays
using RasterDataSources
import ArchGDAL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# The niche axis of a threaded trait / layer (`Unclassified` when none was declared).
axisof(::NicheTolerance{A}) where {A} = A
axisof(::ContinuousLayer{R, A}) where {R, A} = A
axisof(::DiscreteLayer{A}) where {A} = A
axisof(_) = Unclassified

# A trait may only be matched to a layer on the same niche axis (`Unclassified` matches
# anything, for the un-threaded/back-compat paths).
function _checkaxis(t, h)
    at, ah = axisof(t), axisof(h)
    (at === Unclassified || ah === Unclassified || at === ah) ||
        error("A trait on the `$at` axis cannot be matched to a layer on the `$ah` axis — check the order (or axes) of the species niches against the environment layers.")
    return nothing
end

# Infer the trait–environment nichefit from the trait type, checking the niche axes.
# A `NicheTolerance`'s response distribution is built per-species in `_suitability`; its nichefit's `NF`
# is the tolerance's own frame (what the distribution's parameters are actually expressed in) — not the
# regime's incidental unit, so a genuine tolerance/regime unit disagreement is caught by `tematch`
# (Ecosystem construction) instead of nichefit silently mirroring whatever the regime happens to be.
function _default_suitability(t::NicheTolerance, regime)
    _checkaxis(t, regime)
    return NicheSuitability{eltype(t)}()
end
# Same reasoning as the `NicheTolerance` method above: `NF` comes from the tolerance, not the regime.
function _default_suitability(t::DiscreteTolerance, regime)
    return MatchSuitability{eltype(t)}()
end
function _default_suitability(t::LandCoverTolerance, regime)
    return LandCoverSuitability{eltype(t)}()
end
# A collection of tolerances over a matching collection of regimes infers each sub-tolerance's
# nichefit and combines them multiplicatively.
function _default_suitability(t::ToleranceCollection2, h::RegimeCollection2)
    return multiplicativeFit2(_default_suitability(t.one, h.one),
                              _default_suitability(t.two, h.two))
end
function _default_suitability(t::ToleranceCollection3, h::RegimeCollection3)
    return multiplicativeFit3(_default_suitability(t.one, h.one),
                              _default_suitability(t.two, h.two),
                              _default_suitability(t.three, h.three))
end
function _default_suitability(traits::AbstractTolerance, regime)
    return error("Cannot infer a trait nichefit for traits of type $(typeof(traits)); pass one explicitly with `nichefit = …`.")
end

# ---------------------------------------------------------------------------
# Public builders
# ---------------------------------------------------------------------------

# Turn a physical grid extent (x, y lengths) and a cell `size` into an integer
# `(nx, ny)` cell count, warning if the extent is not close to a whole number of
# cells. Returns both the cell counts and the total area for the `*AE` builders.
function _extent_dims(extent::Tuple{<:Unitful.Length, <:Unitful.Length},
                      size::Unitful.Length)
    rx = uconvert(NoUnits, extent[1] / size)
    ry = uconvert(NoUnits, extent[2] / size)
    nx = round(Int, rx)
    ny = round(Int, ry)
    (nx ≥ 1 && ny ≥ 1 && isapprox(rx, nx; rtol = 1.0e-2) &&
     isapprox(ry, ny; rtol = 1.0e-2)) ||
        @warn "Grid extent $(extent[1]) × $(extent[2]) is not close to a whole number of $size cells; using a $nx × $ny grid."
    return (Int64(nx), Int64(ny)), float(extent[1] * extent[2])
end

# `_canonical(value, axis)` (axis-driven, in `src/Layer.jl`) canonicalises a regime value to its axis's
# unit; `_reaxis(layer, axis)` (also `src/Layer.jl`) tags the layer with that axis.

# Generic linear change of a continuous regime: every cell shifts by `rate × timestep`,
# for any unit. Subsumes `TempChange`/`RainfallChange`, which differ only in this
# canonicalising `uconvert`, and neither of which is ever dispatched on.
function LinearChange(eco::AbstractEcosystem, regime::ContinuousRegime,
                      timestep::Unitful.Time)
    v = uconvert(unit(eltype(regime.matrix)) / unit(timestep),
                 regime.dynamics.rate)
    return regime.matrix .+= v * timestep
end

# Build a top→bottom linear gradient regime for any unit — one function (cf. the
# duplicated `tempgrad`/`raingrad`), since the value is already canonicalised here.
function _gradient_regime(low, high, size, dim, rate, axis::Type{<:NicheAxis})
    dim[1] > 1 || error("First grid dimension must be > 1 for a gradient")
    lo = _canonical(low, axis)
    hi = _canonical(high, axis)
    M = Array{typeof(lo)}(undef, dim)
    grad = collect(range(lo, stop = hi, length = dim[1]))
    for i in 1:dim[1]
        M[i, :] .= grad[i]
    end
    return ContinuousRegime(M, size,
                            LayerUpdate(LinearChange, float(rate),
                                        typeof(dimension(lo))))
end

# --- Materialise a synthetic spec into a regime on a grid -------------------
# `materialise(spec, dim, size)` produces today's `ContinuousRegime`/`DiscreteRegime` (stage 1
# of the layer redesign: no downstream type changes yet). Gradients go through one
# generic path for any unit, not the unit-specific `tempgrad`/`raingrad`. `_reaxis` (tag a layer with its
# axis) and `_canonical` (canonicalise a value to its axis's unit) live in `src/Layer.jl`.

function materialise(s::UniformSpec{A}, dim::Tuple,
                     size::Unitful.Length) where {A}
    return _reaxis(simpleregime(_canonical(s.value, A), size, dim, A), A)
end
function materialise(s::GradientSpec{A}, dim::Tuple,
                     size::Unitful.Length) where {A}
    return _reaxis(_gradient_regime(s.low, s.high, size, dim, s.rate, A), A)
end
function materialise(s::PeakedSpec{A}, dim::Tuple,
                     size::Unitful.Length) where {A}
    lo = _canonical(s.low, A)
    hi = _canonical(s.high, A)
    regime = _gradient_regime(lo, hi + (hi - lo), size, dim, s.rate, A)
    regime.matrix[(ceil(Int, dim[1] / 2) + 1):end, :] = regime.matrix[floor(Int,
                                                                            dim[1] / 2):-1:1,
                                                                      :]
    return _reaxis(regime, A)
end
function materialise(s::NicheSpec{A}, dim::Tuple,
                     size::Unitful.Length) where {A}
    n = s.numniches
    return _reaxis(randomniches(dim, collect(1:n), 0.5, fill(1.0 / n, n), size),
                   A)
end

# Cell counts, total area (km²) and cell side (km) of an (extent, cell-size) grid.
function _grid_geometry(extent, size)
    dim, area = _extent_dims(extent, size)
    a = uconvert(km^2, float(area))
    return dim, a, sqrt(a / (dim[1] * dim[2]))
end

# Assemble a regime + a per-area `maxsupply` supply + active mask into a GridHabitat.
function _synthetic_env(regime, maxsupply, area, active, dim)
    B = cancel(float(maxsupply), area)
    T = supplytype(typeof(B))
    supply = fill(B / (dim[1] * dim[2]), dim)
    act = active === nothing ? fill(true, dim) : active
    return GridHabitat{typeof(regime), T}(regime, act, T(supply))
end

"""
    build_environment(extent, size, regime::Unitful.Temperature; supply = 1.0kJ/(km^2*day), active = nothing)
    build_environment(extent, size, regime::Tuple; rate = 0.0K/month, peaked = false, supply = 1.0kJ/(km^2*day), active = nothing)
    build_environment(extent, size, regime::Integer; supply = 1.0kJ/(km^2*day), active = nothing)

Build a `GridHabitat` over a grid of physical `extent` (an `(x, y)` tuple of
lengths, e.g. `(10km, 10km)`) divided into square cells of side `size`. The cell
counts are derived as `extent ./ size` (a warning is issued if the extent is not
close to a whole number of cells).

The third argument chooses the kind of environment by its type (multiple
dispatch, no pattern flag):

  - a `Unitful.Temperature` (e.g. `298.0K`) — a flat temperature regime;
  - a pair of temperatures `(low, high)` — a temperature gradient from `low` at
    the bottom to `high` at the top, changing over time at `rate`
    (default `0.0K/month`); pass `peaked = true` for a gradient that instead
    peaks in the middle of the grid;
  - a pair of lengths `(low, high)` in mm — a rainfall gradient; if `supply` is
    omitted the rainfall values themselves become a [`WaterSupply`](@ref), interpreted
    as delivered per day;
  - an `Integer` — that many discrete random niches.

`supply` is a per-area maximum resource rate whose dimension selects the supply type
(`kJ/km^2/day` → solar, `L/km^2/day` → water, …), or a pre-built `AbstractSupply`.
`active` is an optional `Matrix{Bool}` of live cells (default: all active).
"""
function build_environment(extent::Tuple{<:Unitful.Length, <:Unitful.Length},
                           size::Unitful.Length,
                           regime::Union{Unitful.Temperature, AbstractFloat};
                           axis::Type{<:NicheAxis} = Unclassified,
                           supply = 1.0 *
                                    canonicalunit(Resource, SolarRadiation()) /
                                    km^2,
                           active::Union{Nothing, Matrix{Bool}} = nothing)
    dim, area, cs = _grid_geometry(extent, size)
    regime = materialise(UniformSpec(axis, regime), dim, cs)
    return _synthetic_env(regime, supply, area, active, dim)
end

function build_environment(extent::Tuple{<:Unitful.Length, <:Unitful.Length},
                           size::Unitful.Length,
                           regime::Tuple{<:Unitful.Temperature,
                                         <:Unitful.Temperature};
                           rate = 0.0K / month,
                           peaked::Bool = false,
                           axis::Type{<:NicheAxis} = Unclassified,
                           supply = 1.0 *
                                    canonicalunit(Resource, SolarRadiation()) /
                                    km^2,
                           active::Union{Nothing, Matrix{Bool}} = nothing)
    dim, area, cs = _grid_geometry(extent, size)
    spec = peaked ? PeakedSpec(axis, regime[1], regime[2], rate) :
           GradientSpec(axis, regime[1], regime[2], rate)
    regime = materialise(spec, dim, cs)
    return _synthetic_env(regime, supply, area, active, dim)
end

function build_environment(extent::Tuple{<:Unitful.Length, <:Unitful.Length},
                           size::Unitful.Length,
                           regime::Tuple{<:Unitful.Length, <:Unitful.Length};
                           rate = 0.0mm / month,
                           axis::Type{<:NicheAxis} = Unclassified,
                           supply = nothing,
                           active::Union{Nothing, Matrix{Bool}} = nothing)
    dim, area, cs = _grid_geometry(extent, size)
    regime = materialise(GradientSpec(axis, regime[1], regime[2], rate), dim,
                         cs)
    if supply === nothing
        # the rainfall values themselves become the water supply, interpreted as "this
        # many mm delivered per day" and converted to an absolute L/day via `cancel`
        # against the cell area — the same areal-rate-to-absolute path every other
        # supply construction in this file goes through.
        supply = WaterSupply(cancel.(regime.matrix ./ day, Ref(cs^2)))
        act = active === nothing ? fill(true, dim) : active
        return GridHabitat{typeof(regime), typeof(supply)}(regime, act, supply)
    end
    return _synthetic_env(regime, supply, area, active, dim)
end

function build_environment(extent::Tuple{<:Unitful.Length, <:Unitful.Length},
                           size::Unitful.Length,
                           regime::Integer;
                           axis::Type{<:NicheAxis} = Unclassified,
                           supply = 1.0 *
                                    canonicalunit(Resource, SolarRadiation()) /
                                    km^2,
                           active::Union{Nothing, Matrix{Bool}} = nothing)
    dim, area, cs = _grid_geometry(extent, size)
    regime = materialise(NicheSpec(axis, regime), dim, cs)
    return _synthetic_env(regime, supply, area, active, dim)
end

"""
    build_environment(spec::AbstractSyntheticSpec, extent, size; supply = 1.0kJ/(km^2*day),
        active = nothing)

Build a `GridHabitat` from a synthetic layer `spec` ([`UniformSpec`](@ref),
[`GradientSpec`](@ref), [`PeakedSpec`](@ref) or [`NicheSpec`](@ref)) over a grid of
physical `extent` in square cells of side `size`, with a per-area `supply`.
"""
function build_environment(spec::AbstractSyntheticSpec,
                           extent::Tuple{<:Unitful.Length, <:Unitful.Length},
                           size::Unitful.Length;
                           supply = 1.0 *
                                    canonicalunit(Resource, SolarRadiation()) /
                                    km^2,
                           active::Union{Nothing, Matrix{Bool}} = nothing)
    dim, area, cs = _grid_geometry(extent, size)
    regime = materialise(spec, dim, cs)
    return _synthetic_env(regime, supply, area, active, dim)
end

# ---------------------------------------------------------------------------
# Data-driven environments: reconcile one or more rasters onto a common grid
# and build regime + supply + active mask from them.
# ---------------------------------------------------------------------------

# Degrees north (latitude) and east (longitude) of a raster's cell centres.
_latvals(raster::ClimateRaster) = raster.array.axes[1].val
_longvals(raster::ClimateRaster) = raster.array.axes[2].val

# Land cover is the one discrete data regime; every other source is continuous.
_isdiscrete(::ClimateRaster{<:EarthEnv{<:LandCover}}) = true
_isdiscrete(::ClimateRaster) = false

# Convert a physical cell side to a step in degrees (≈111.32 km per degree, as
# the existing `*AE` constructors assume), and back.
_deg(side::Unitful.Length) = uconvert(NoUnits, side / 111.32km) * °
_side(step_deg) = abs(uconvert(NoUnits, step_deg / °)) * 111.32km

# Intersection (lo, hi) of several axis-value ranges; errors if they are disjoint.
function _overlap(ranges)
    lo = maximum(minimum.(ranges))
    hi = minimum(maximum.(ranges))
    lo < hi ||
        error("the rasters do not overlap in space, so there is no common region to simulate — pass an explicit `region = LatLong(lat, long)`.")
    return lo, hi
end

# Nearest source index for each target value (both in degrees).
_nearest(src, targets) = [argmin(abs.(src .- t)) for t in targets]

# The common target grid (target latitude range, target longitude range) shared
# by all `regimes`: region defaults to their intersection, cell size to the
# first regime's native resolution unless `size` is given.
function _target_grid(regimes::Tuple; size = nothing, region = nothing)
    if region === nothing
        latlo, lathi = _overlap([_latvals(h) for h in regimes])
        longlo, longhi = _overlap([_longvals(h) for h in regimes])
    else
        latlo, lathi = minimum(region.lat), maximum(region.lat)
        longlo, longhi = minimum(region.long), maximum(region.long)
    end
    if size === nothing
        # No target size given: use the primary raster's own cell centres inside
        # the region, so the common case resamples nothing.
        lat0 = _latvals(first(regimes))
        long0 = _longvals(first(regimes))
        return filter(x -> latlo <= x <= lathi, lat0),
               filter(x -> longlo <= x <= longhi, long0)
    else
        step_deg = _deg(size)
        nlat = max(2, round(Int, (lathi - latlo) / step_deg))
        nlong = max(2, round(Int, (longhi - longlo) / step_deg))
        return range(latlo, lathi; length = nlat),
               range(longlo, longhi; length = nlong)
    end
end

# Sample a raster's data onto the target grid by nearest neighbour, warning if
# this actually changes its resolution.
function _sample_onto(raster::ClimateRaster, tlat, tlong; name = "raster")
    A = Array(raster.array)
    li = _nearest(_latvals(raster), tlat)
    lj = _nearest(_longvals(raster), tlong)
    out = ndims(A) == 2 ? A[li, lj] : A[li, lj, :]
    srcstep = abs(uconvert(NoUnits,
                           (_latvals(raster)[2] - _latvals(raster)[1]) / °))
    tgtstep = length(tlat) > 1 ?
              abs(uconvert(NoUnits, (tlat[2] - tlat[1]) / °)) :
              srcstep
    isapprox(srcstep, tgtstep; rtol = 1.0e-3) ||
        @warn "Resampling $name from ~$(round(srcstep, digits = 4))° to ~$(round(tgtstep, digits = 4))° cells to match the simulation grid — check this is intended."
    return out
end

# Build the regime (continuous or discrete) from a raster on the target grid, tagged with
# niche axis `axis` (from the shipped layer table via the source `code`; `Unclassified` for a
# bare raster whose code is unknown). The low-level `*Hab` constructors default to
# `Unclassified`, so `_reaxis` re-tags the result with the real axis (sharing the array).
function _dataregime(raster::ClimateRaster, tlat, tlong,
                     axis::Type{<:NicheAxis} = Unclassified)
    out = _sample_onto(raster, tlat, tlong; name = "regime")
    csize = _cellsize(tlat, tlong)
    if _isdiscrete(raster)
        # discrete land-cover classes: keep the categorical values as read
        regime = DiscreteRegime(out, csize, LayerUpdate(NoChange, 0.0 / s))
    else
        # canonicalise the read values (in their actual raster unit, e.g. °C) to the axis's
        # canonical unit (e.g. K) — a proper affine conversion for an absolute Temperature,
        # a no-op for an interval/dimensionless axis whose data is already canonical.
        regime = ContinuousRegime(_canonical.(out, Ref(axis)), csize,
                                  LayerUpdate(NoChange, 0.0 / s,
                                              Unitful.Dimensions{()}))
    end
    return _reaxis(regime, axis)
end

# Wrap a sampled supply layer as a supply: `cancel` (dispatched on dimension, any native
# time unit) converts the raw per-area rate to an absolute per-cell one against `cellarea`,
# and `supplytype` (also dispatched, on the result's dimension) picks Solar/Water/Simple — a
# monthly (3-D) stack becomes the time-varying variant.
function _wrap_supply(out, cellarea)
    abs = cancel.(out, Ref(cellarea))
    T = supplytype(eltype(abs))
    return ndims(abs) == 2 ? T(abs) : T(abs, 1)
end

# Resolve the `supply` keyword into an AbstractSupply on the target grid — one method per
# recognised kind of `supply`, dispatched on its type rather than branched on at runtime.
_resolve_supply(supply::AbstractSupply, tlat, tlong) = supply
function _resolve_supply(supply::ClimateRaster, tlat, tlong)
    cellarea = _cellsize(tlat, tlong)^2
    return _wrap_supply(_sample_onto(supply, tlat, tlong; name = "supply"),
                        cellarea)
end
function _resolve_supply(supply, tlat, tlong)
    dims = (length(tlat), length(tlong))
    maxsupply = float(supply === nothing ?
                      1.0 * canonicalunit(Resource, SolarRadiation()) /
                      km^2 : supply)
    percell = cancel(maxsupply, _side(tlat[2] - tlat[1])^2)
    T = supplytype(typeof(percell))
    return T(fill(percell, dims))
end

# --- Active-area masks -----------------------------------------------------
# Raster-derived masks (`landmask`, `datamask`) return an `AxisArray{Bool}` that
# carries its own lat/long grid, so it can be resampled onto the environment's
# grid at build time. Geometric masks are lazy specs resolved against the grid.

"""
    datamask(raster::ClimateRaster)

Mask of the cells of `raster` that hold data (are not missing/`NaN`), for use as the
`active =` argument of [`build_environment`](@ref).
"""
function datamask(raster::ClimateRaster)
    A = Array(raster.array)
    layer = ndims(A) == 2 ? A : A[:, :, 1]
    return AxisArray(Matrix{Bool}(.!isnan.(layer)),
                     raster.array.axes[1], raster.array.axes[2])
end

"""
    landmask(landcover::ClimateRaster{<:EarthEnv{<:LandCover}}; water = [4])

Mask of the non-water cells of a land-cover raster `landcover` (excluding the classes in
`water`), for use as `active =` in [`build_environment`](@ref).
"""
function landmask(landcover::ClimateRaster{<:EarthEnv{<:LandCover}};
                  water = [4])
    A = Array(landcover.array)
    layer = ndims(A) == 2 ? A : A[:, :, 1]
    return AxisArray(Matrix{Bool}(map(v -> !(v in water), layer)),
                     landcover.array.axes[1], landcover.array.axes[2])
end

struct CircleMask{C, R <: Unitful.Length}
    centre::C          # (lat, long) in degrees, or `nothing` for the grid centre
    radius::R
end

"""
    circlemask(; radius, centre = nothing)

A circular active-area mask of the given `radius` (a length) about `centre`
(`(lat, long)` degrees, default the grid centre), resolved onto the grid by
[`build_environment`](@ref).
"""
function circlemask(; radius::Unitful.Length, centre = nothing)
    return CircleMask(centre, radius)
end

struct ShapeMask
    geoms::Vector{ArchGDAL.IPreparedGeometry}   # reprojected to WGS84 lat/long (traditional x/y order)
end

"""
    shapemask(path::AbstractString; layer = 0)

An active-area mask from every polygon feature in `layer` (0-indexed) of the vector file `path` (a
shapefile, GeoJSON, GeoPackage, or any format GDAL/OGR reads — a `path` ending in `.zip` is read directly,
no need to unzip first) — a cell is active if its centre falls inside *any* feature. Reprojected to WGS84
lat/long if the file's CRS differs (e.g. a British National Grid shapefile), resolved onto the grid by
[`build_environment`](@ref), like [`circlemask`](@ref).
"""
function shapemask(path::AbstractString; layer::Integer = 0)
    vpath = endswith(path, ".zip") ? "/vsizip/" * path : path
    dataset = ArchGDAL.read(vpath)
    lyr = ArchGDAL.getlayer(dataset, layer)
    sr = ArchGDAL.getspatialref(lyr)
    # `order = :trad` forces traditional (long, lat) GIS axis order on the target CRS — EPSG:4326's
    # authority-defined order is (lat, long), which would silently mismatch `createpoint(long, lat)` below.
    wgs84 = ArchGDAL.importEPSG(4326; order = :trad)
    geoms = map(lyr) do feature
        geom = ArchGDAL.getgeom(feature)
        # A missing `.prj` (no CRS metadata) gives a null spatial ref; assume already WGS84 lat/long.
        if sr.ptr != C_NULL
            ArchGDAL.createcoordtrans(sr, wgs84) do ct
                return ArchGDAL.transform!(geom, ct)
            end
        end
        return ArchGDAL.preparegeom(geom)
    end
    return ShapeMask(geoms)
end

# Mirrors `_circle`: a cell is active if its centre falls inside any of the shapefile's features.
function _shape(sm::ShapeMask, tlat, tlong)
    mask = falses(length(tlat), length(tlong))
    for (i, lat) in enumerate(tlat), (j, long) in enumerate(tlong)
        point = ArchGDAL.createpoint(ustrip(°, long), ustrip(°, lat))
        mask[i, j] = any(g -> ArchGDAL.contains(g, point), sm.geoms)
    end
    return Matrix{Bool}(mask)
end

# Nearest-sample a bare AxisArray (e.g. a Bool mask) onto the target grid.
function _sample_axes(A::AxisArray, tlat, tlong)
    li = _nearest(A.axes[1].val, tlat)
    lj = _nearest(A.axes[2].val, tlong)
    return Matrix{Bool}(Array(A)[li, lj])
end

# Great-circle-ish planar distance (km) from each grid cell to a centre.
function _circle(cm::CircleMask, tlat, tlong)
    clat = cm.centre === nothing ? (first(tlat) + last(tlat)) / 2 : cm.centre[1]
    clong = cm.centre === nothing ? (first(tlong) + last(tlong)) / 2 :
            cm.centre[2]
    mask = falses(length(tlat), length(tlong))
    for (i, lat) in enumerate(tlat), (j, long) in enumerate(tlong)
        dlat = _side(lat - clat)
        dlong = _side(long - clong) * cos(clat)
        mask[i, j] = hypot(dlat, dlong) <= cm.radius
    end
    return Matrix{Bool}(mask)
end

# Default active mask for a regime: a continuous regime's non-missing cells, all cells
# for a discrete one, and (for a collection) the cells that are live in every sub-regime.
function _defaultactive(regime::AbstractRegime)
    return iscontinuous(regime) ? Matrix{Bool}(.!isnan.(regime.matrix)) :
           fill(true, Base.size(regime.matrix))
end
function _defaultactive(regime::RegimeCollection2)
    return Matrix{Bool}(_defaultactive(regime.one) .&
                        _defaultactive(regime.two))
end
function _defaultactive(regime::RegimeCollection3)
    return Matrix{Bool}(_defaultactive(regime.one) .&
                        _defaultactive(regime.two) .&
                        _defaultactive(regime.three))
end

# Resolve the `active` keyword into a Matrix{Bool} on the target grid — one method per
# recognised kind of `active`, dispatched on its type rather than branched on at runtime.
_resolve_active(active::Nothing, regime, tlat, tlong) = _defaultactive(regime)
function _resolve_active(active::AxisArray, regime, tlat, tlong)
    return _sample_axes(active, tlat, tlong)
end
function _resolve_active(active::CircleMask, regime, tlat, tlong)
    return _circle(active, tlat, tlong)
end
function _resolve_active(active::ShapeMask, regime, tlat, tlong)
    return _shape(active, tlat, tlong)
end
function _resolve_active(active::AbstractMatrix{Bool}, regime, tlat, tlong)
    dims = (length(tlat), length(tlong))
    Base.size(active) == dims ||
        error("`active` is $(Base.size(active)) but the grid is $dims")
    return Matrix{Bool}(active)
end
function _resolve_active(active, regime, tlat, tlong)
    return error("unrecognised `active` argument of type $(typeof(active)); use nothing, a Matrix{Bool}, or a mask (landmask/datamask/circlemask/shapemask).")
end

"""
    build_environment(regime_data::ClimateRaster; supply = nothing, size = nothing,
        active = nothing, region = nothing)

Build a `GridHabitat` from a climate/land-cover raster. The grid comes from the
raster: `region` defaults to the raster's extent (the intersection of all rasters, for
the tuple/multi form) and cell `size` defaults to the raster's native resolution
(pass `size` to resample). `supply` accepts a per-area quantity, a prebuilt
`AbstractSupply`, or a data layer (its units pick solar/water, a 12-month
stack giving the time-varying variant). `active` defaults to the non-missing cells of
the regime, or takes a `Matrix{Bool}` / mask (see [`landmask`](@ref)/[`shapemask`](@ref) etc.).

`axis` tags the regime layer with a [`NicheAxis`](@ref) so it is matched to species niches
by name; it defaults to `Unclassified` here (a bare raster has no layer code to look up), but
the `source`/`layer` and [`SourceSpec`](@ref) builders derive it from the shipped layer table.
"""
function build_environment(regime_data::ClimateRaster;
                           axis::Type{<:NicheAxis} = Unclassified,
                           supply = nothing,
                           size = nothing,
                           active = nothing,
                           region = nothing)
    tlat, tlong = _target_grid((regime_data,); size = size, region = region)
    regime = _dataregime(regime_data, tlat, tlong, axis)
    supply = _resolve_supply(supply, tlat, tlong)
    act = _resolve_active(active, regime, tlat, tlong)
    return GridHabitat{typeof(regime), typeof(supply)}(regime, act, supply)
end

# Re-wrap a raster's data with a unit, preserving the axes and the source type (so
# `_isdiscrete` still dispatches correctly). `NoUnits` is a genuine no-op that leaves
# the data as bare numbers, which is what discrete (land cover) regimes want.
function _attach_unit(raster::ClimateRaster{S}, u) where {S}
    A = raster.array
    return ClimateRaster(S, AxisArray(A.data .* u, AxisArrays.axes(A)...))
end

# Default a `SourceSpec`'s unit from the shipped layer table when not given explicitly.
function SourceSpec(source::Type, code)
    return SourceSpec(source, code, layerunit(source, code))
end

# Read a `SourceSpec` into a unit-attached `ClimateRaster` (the eager step deferred by the
# lazy descriptor). Extra keywords (e.g. `month`) forward to `read`.
function _read(sl::SourceSpec; kw...)
    return _attach_unit(read(sl.source, sl.code; kw...), sl.unit)
end

"""
    build_environment(source::Type{<:RasterDataSources.RasterDataSource}, layer;
        unit = layerunit(source, layer), region = nothing, supply = nothing,
        size = nothing, active = nothing, kw...)
    build_environment(spec::SourceSpec; region = nothing, supply = nothing,
        size = nothing, active = nothing, kw...)

Read a single `layer` (an `Int` layer number or a `Symbol`/`String` layer key) from a
RasterDataSources `source` (e.g. `build_environment(WorldClim{BioClim}, 1)`), attach its
physical unit from the shipped layer table (temperature comes back in `K`, precipitation
in `mm`, land cover dimensionless) and build a [`GridHabitat`](@ref) from it. Pass
`unit =` to override the table lookup, or construct a [`SourceSpec`](@ref) directly. The
regime layer is tagged with the [`NicheAxis`](@ref) named for this `source`/`layer` in the
shipped layer table (e.g. WorldClim BioClim 1 → `Temperature`), so it is matched to species
niches by name; pass `axis =` to override. Extra keywords (`kw...`, e.g. `month`) are forwarded
to `read`; `region`, `supply`, `size` and `active` behave as in the [`ClimateRaster`](@ref)
method.
"""
function build_environment(source::Type{<:RasterDataSources.RasterDataSource},
                           layer;
                           unit = layerunit(source, layer),
                           kw...)
    return build_environment(SourceSpec(source, layer, unit); kw...)
end

function build_environment(spec::SourceSpec;
                           axis::Type{<:NicheAxis} = _specaxis(spec),
                           region = nothing,
                           supply = nothing,
                           size = nothing,
                           active = nothing,
                           kw...)
    return build_environment(_read(spec; kw...); axis = axis, region = region,
                             supply = supply, size = size, active = active)
end

# The declared niche axis of a data source element, from the shipped layer table (via the
# source `code`), or `Unclassified` when the source has no table axis — or when the code is
# unknown (a bare `ClimateRaster`, which has dropped it).
function _specaxis(spec::SourceSpec)
    return something(layeraxis(spec.source, spec.code),
                     Unclassified)
end
_specaxis(spec::Tuple) = something(layeraxis(spec[1], spec[2]), Unclassified)
_specaxis(::ClimateRaster) = Unclassified

