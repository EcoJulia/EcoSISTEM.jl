# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The layer recipes — what a caller writes to declare a layer before anything is built. The
# synthetic ones describe a grid outright; `Varying` wraps any spec to say how it changes in time.

using Unitful
using Unitful.DefaultSymbols

"""
    AbstractSpec

Abstract root of every build-time layer/mask recipe. A spec holds no grid array itself;
[`materialise`](@ref EcoSISTEM.materialise) (or the data-driven resolvers) turns it into a real grid
layer only at [`GridHabitat`](@ref) time. Two branches: [`AbstractLazySpec`](@ref) (data-backed,
resolved against the target grid, usable as *either* a layer or a mask) and
[`AbstractSyntheticSpec`](@ref) (self-contained generators, role-fixed).
"""
abstract type AbstractSpec end

"""
    AbstractSyntheticSpec <: AbstractSpec

Abstract supertype of the synthetic (generated, non-data) specs, split by role into
[`AbstractSyntheticLayerSpec`](@ref) (regime/supply layers) and [`AbstractSyntheticMaskSpec`](@ref)
(active masks).
"""
abstract type AbstractSyntheticSpec <: AbstractSpec end

"""
    AbstractSyntheticLayerSpec <: AbstractSyntheticSpec

Synthetic layer recipes — [`UniformSpec`](@ref), [`GradientSpec`](@ref), [`PeakedSpec`](@ref),
[`NicheSpec`](@ref) — used as a regime or supply, never as a mask.
"""
abstract type AbstractSyntheticLayerSpec <: AbstractSyntheticSpec end

"""
    AbstractSyntheticMaskSpec <: AbstractSyntheticSpec

Synthetic active-mask recipes — [`CircleMaskSpec`](@ref) — used only as an `active` mask.
"""
abstract type AbstractSyntheticMaskSpec <: AbstractSyntheticSpec end

"""
    UniformSpec{A <: NicheAxis, V}

A layer with a single constant value in every cell.

# Arguments

  - `value`: what every cell holds — a `Unitful` quantity, or a bare number.
  - `axis`: the niche axis this layer is on. **Required**: it is what makes the layer matchable
    against a species' tolerances, and no default could be right for every model. Pass
    [`NicheAxis`](@ref) itself for data whose meaning is genuinely not being claimed, such as a mask.
"""
struct UniformSpec{A <: NicheAxis, V} <: AbstractSyntheticLayerSpec
    value::V
    # The sole constructor — axis is a keyword (not a leading positional), which is what lets
    # this collapse to one constructor with no separate no-axis overload.
    function UniformSpec(value::V;
                         axis::Type{A}) where {A <: NicheAxis, V}
        return new{A, V}(value)
    end
end

"""
    GradientSpec{A <: NicheAxis, V}

A layer whose value runs linearly across the grid.

The pattern is fixed in space and does not change of itself; wrap the spec in [`Varying`](@ref) to
give it a change in time.

# Arguments

  - `low`, `high`: the values at the two ends. `low` must be strictly less than `high`.
  - `axis`: the niche axis this layer is on, as [`UniformSpec`](@ref). Required.
  - `orientation`: which way the gradient runs, as a compass bearing clockwise from North. `0°`, the
    default, runs south to north; `90°` west to east; `180°` north to south; and so on for any
    bearing between. Must be a real angle `Quantity` — `90°`, or `(pi/2)rad` — since a bare number
    would be read as radians by Unitful's own counter-intuitive convention. Stored in degrees, so the
    type is not parametric on the angle unit.
"""
struct GradientSpec{A <: NicheAxis, V} <: AbstractSyntheticLayerSpec
    low::V
    high::V
    orientation::typeof(1.0°)
    # The sole constructor — `axis` is a keyword (not a leading positional), which is what lets
    # this collapse to one constructor with no separate no-axis overload, since defining it here
    # suppresses Julia's default inner constructor: `low < high` and the canonicalising
    # `orientation` check can't be bypassed by constructing a `GradientSpec{A, V}` directly.
    function GradientSpec(low::V, high::V;
                          axis::Type{A},
                          orientation::Unitful.Quantity{<:Real, Unitful.NoDims} = 0.0°) where {A <:
                                                                                               NicheAxis,
                                                                                               V}
        low < high ||
            error("`low` ($low) must be strictly less than `high` ($high) for a `GradientSpec`.")
        return new{A, V}(low, high, float(uconvert(°, orientation)))
    end
end

"""
    PeakedSpec{A <: NicheAxis, V}

A layer whose value peaks — or dips — in the middle of the grid: `inside` at the centre, falling off
to `outside` at each end.

As [`GradientSpec`](@ref), the pattern is fixed in space; wrap the spec in [`Varying`](@ref) to give
it a change in time.

# Arguments

  - `outside`, `inside`: the values at the ends and at the centre. They must differ, or there is no
    peak.
  - `axis`: the niche axis this layer is on. Required.
  - `orientation`: which way the peak runs, with the same meaning and the same angle-`Quantity`
    requirement as [`GradientSpec`](@ref)'s.
"""
struct PeakedSpec{A <: NicheAxis, V} <: AbstractSyntheticLayerSpec
    outside::V
    inside::V
    orientation::typeof(1.0°)
    # See `GradientSpec`'s inner constructor comment: the sole constructor, so it can't be
    # bypassed with a bare-number `orientation` or with `outside ≈ inside`.
    function PeakedSpec(outside::V, inside::V;
                        axis::Type{A},
                        orientation::Unitful.Quantity{<:Real, Unitful.NoDims} = 0.0°) where {A <:
                                                                                             NicheAxis,
                                                                                             V}
        outside ≈ inside &&
            error("`outside` ($outside) and `inside` ($inside) must differ for a `PeakedSpec` (otherwise there is no peak).")
        return new{A, V}(outside, inside, float(uconvert(°, orientation)))
    end
end

"""
    NicheSpec{A <: NicheAxis}

A **categorical** layer of randomly assigned niche classes — the synthetic counterpart of a land
cover map.

# Arguments

  - `n`: how many distinct classes to draw from.
  - `axis`: the niche axis this layer is on. Required.
"""
struct NicheSpec{A <: NicheAxis} <: AbstractSyntheticLayerSpec
    numniches::Int64
    function NicheSpec(n::Integer;
                       axis::Type{A}) where {A <: NicheAxis}
        return new{A}(Int64(n))
    end
end

"""
    CircleMaskSpec(; radius, centre = nothing)

A circular active-area mask: the cells within a given distance of a point are simulated and the rest
are not. Resolved onto the grid when the habitat is built.

# Arguments

  - `radius`: how far the circle reaches, as a **physical distance** (`100.0km`) rather than an angle.
  - `centre`: where it is centred, as a [`LatLong`](@ref). Defaults to the grid's own centre — and on
    a **synthetic** grid, which has no real-world coordinates, that default is the only thing
    supported, since a geographic centre could not be placed.
"""
struct CircleMaskSpec{V <: Unitful.Length} <: AbstractSyntheticMaskSpec
    centre::Union{LatLong, Nothing}   # a lat/long point, or `nothing` for the grid centre
    radius::V
    # The sole constructor — `V` is inferable from the `radius` keyword's value in the `where`
    # clause regardless of it being a keyword, so this collapses to one constructor with no
    # separate outer method.
    function CircleMaskSpec(; radius::V,
                            centre::Union{LatLong, Nothing} = nothing) where {V <:
                                                                              Unitful.Length}
        return new{V}(centre, radius)
    end
end
# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# **A spec's one-liner is the expression that builds it.** A spec is a recipe the user typed, so
# echoing it back is both the most useful thing to show and self-documenting — it can be pasted
# straight into a session. That is a different rule from the layer displays in `Layer.jl`, which
# describe a built object rather than a recipe, and it is why these carry `axis = …` rather than
# folding the axis into the type name.
#
# Optional arguments appear only when they are not at their default, for the same reason: a printed
# `orientation = 0.0°` is noise in an expression nobody would type.
function Base.show(io::IO, spec::UniformSpec{A}) where {A}
    return print(io, "UniformSpec($(spec.value), axis = $(nameof(A)))")
end

function Base.show(io::IO, spec::GradientSpec{A}) where {A}
    orient = iszero(spec.orientation) ? "" :
             ", orientation = $(spec.orientation)"
    return print(io,
                 "GradientSpec($(spec.low), $(spec.high), axis = $(nameof(A))$(orient))")
end

function Base.show(io::IO, spec::PeakedSpec{A}) where {A}
    orient = iszero(spec.orientation) ? "" :
             ", orientation = $(spec.orientation)"
    return print(io,
                 "PeakedSpec($(spec.outside), $(spec.inside), axis = $(nameof(A))$(orient))")
end

function Base.show(io::IO, spec::NicheSpec{A}) where {A}
    return print(io, "NicheSpec($(spec.numniches), axis = $(nameof(A)))")
end

function Base.show(io::IO, spec::CircleMaskSpec)
    centre = isnothing(spec.centre) ? "" : ", centre = $(spec.centre)"
    return print(io, "CircleMaskSpec(radius = $(spec.radius)$(centre))")
end

# ══ Functions ══════════════════════════════════════════════════════════════════════════════════

# ---------------------------------------------------------------------------
# Layer specs (recipes)
# ---------------------------------------------------------------------------
# A spec describes *how to produce* a gridded layer without holding any grid data;
# `materialise(spec, area)` turns it into a layer on a decided grid. Specs are
# build-time only and never appear in the simulation hot loop.

"""
    SurfaceSpec(fraction = 1.0; axis = SurfaceArea)

A uniform space layer: the `fraction` of each cell that is available, from `0.0`
to all of it (`1.0`, the default).

A thin convenience over [`UniformSpec`](@ref): it fixes the axis to a space one and checks the
fraction is not negative, and nothing more. `SurfaceSpec()` is the whole-cell case, which is what
makes a space supply easy to ask for without its being automatic.

Not capped at `1.0`. A cover fraction cannot exceed a whole cell, but a canopy legitimately can,
since crowns overlap, so the ceiling belongs to the axis's own `bounds` rather than here.

# Arguments

  - `fraction`: the share of each cell that is available, `1.0` by default.
  - `axis`: which stratum, for when more than one is declared. [`SurfaceArea`](@ref) is the only
    [`SpaceAxis`](@ref) leaf today.
"""
function SurfaceSpec(fraction::Real = 1.0;
                     axis::Type{<:SpaceAxis} = SurfaceArea)
    fraction >= 0 ||
        error("a `SurfaceSpec` fraction is the share of a cell that is available, so it cannot be " *
              "negative; got $fraction.")
    return UniformSpec(float(fraction), axis = axis)
end

# --- Reading a spec ----------------------------------------------------------
# What a spec says about itself, asked by the study-area machinery that has to plan around it.

# Does this spec generate its values rather than read them — and so have no grid of its own?
_issyntheticspec(::AbstractSyntheticSpec) = true

_issyntheticspec(_) = false

# An already-built layer carries its axis in its own type, so it needs no resolving at all. That is
# the property that makes passing one in safe.
_specaxis(layer::AbstractLayer) = _layeraxis(layer)

# The niche axis of a synthetic layer spec, which carries it as its first type parameter (unlike a
# `SourceSpec`, whose axis is resolved from the shipped table at construction).
_specaxis(::UniformSpec{A}) where {A} = A

_specaxis(::GradientSpec{A}) where {A} = A

_specaxis(::PeakedSpec{A}) where {A} = A

_specaxis(::NicheSpec{A}) where {A} = A

# A synthetic spec's values over a grid of `dim` cells. `NicheSpec` is categorical, so it gets its
# own random-niche generator rather than a gradient field.
#
# The signature matches `_syntheticsupplyfield`'s exactly, which is the point: between them these are
# the whole synthetic family, and the builder reaches every one of them the same way.
function _specfield(spec::NicheSpec, dim, rowsnorth)
    n = spec.numniches
    return _nichefield(dim, collect(1:n), 0.5, fill(1.0 / n, n))
end

function _specfield(spec::AbstractSyntheticLayerSpec, dim, rowsnorth)
    return _syntheticsupplyfield(spec, dim, rowsnorth)
end

# --- Generating a layer from a synthetic spec ---------------------------------
# What a `UniformSpec`/`GradientSpec`/`PeakedSpec`/`NicheSpec` actually produces on a grid, and the
# mask a `CircleMaskSpec` resolves to. Called from here, from `StudyArea.jl` and from `materialise.jl`.

# Turn a physical grid extent and a cell `size` into an integer `(ny, nx)` cell count, warning if the
# extent is not close to a whole number of cells. Returns both the cell counts and the total area.
#
# The extent is `(y, x)` — north-south first, then east-west — which is the dimension order used
# everywhere else in the package, and the order of the `(ny, nx)` returned here. A square extent
# cannot tell the two apart, so only a non-square grid ever exposes a crossing.
function _extentdims(extent::Tuple{<:Unitful.Length, <:Unitful.Length},
                     size::Unitful.Length)
    ry = uconvert(NoUnits, extent[1] / size)
    rx = uconvert(NoUnits, extent[2] / size)
    ny = round(Int, ry)
    nx = round(Int, rx)
    (nx ≥ 1 && ny ≥ 1 && isapprox(rx, nx, rtol = 1.0e-2) &&
     isapprox(ry, ny, rtol = 1.0e-2)) ||
        @warn "Grid extent $(extent[1]) (y) × $(extent[2]) (x) is not close to a whole number of $size cells; using a $ny × $nx grid."
    return (Int64(ny), Int64(nx)), float(extent[1] * extent[2])
end

# Normalised (0 to 1) linear field over `dim` along compass bearing `orientation` (degrees
# clockwise from North — `0°` runs low row -> high row, `90°` runs low column -> high column
# (west -> east), and so on for any bearing in between). `rowsincreasenorth` says whether
# increasing row index means increasing latitude: `true` (default) for the synthetic grid's
# arbitrary index convention (no real `Y` lookup at all, row 1 = south by documented convention);
# for a *data-driven* grid this must instead reflect the regime's own real `Y` lookup order — a
# GDAL-sourced raster is conventionally north-up (`Y` stored `ReverseOrdered`, row 1 = north),
# the opposite of the synthetic convention — see `_rowsincreasenorth`.
function _gradientfield(dim, orientation, rowsincreasenorth::Bool = true)
    ny, nx = dim
    north, east = cos(orientation), sin(orientation)
    row(i) = rowsincreasenorth ? i : ny + 1 - i
    proj(i, j) = row(i) * north + j * east
    corners = (proj(1, 1), proj(1, nx), proj(ny, 1), proj(ny, nx))
    lo, hi = extrema(corners)
    lo < hi ||
        error("`orientation` ($orientation) gives no variation across a $(ny)×$(nx) grid")
    return [(proj(i, j) - lo) / (hi - lo) for i in 1:ny, j in 1:nx]
end

# Whether increasing row index in a `(Y, X)` dims tuple corresponds to increasing latitude
# (north) — see `_gradientfield`'s docstring comment above for why this matters.
function _rowsincreasenorth(yx)
    return !(DimensionalData.order(yx[1]) isa
             DimensionalData.Lookups.ReverseOrdered)
end

# Cell counts, total area (km²) and cell side (km) of an (extent, cell-size) grid.
function _gridgeometry(extent, size)
    dim, area = _extentdims(extent, size)
    a = uconvert(km^2, float(area))
    return (dim = dim, area = a, cellsize = sqrt(a / (dim[1] * dim[2])))
end

# Resolve the `active` keyword for a SYNTHETIC (`extent`/`size`-based) grid — one method per
# recognised kind, mirroring `_rastermask`'s dispatch-per-kind pattern for the data-driven
# path, but working in cell-index × `cellsize` (Cartesian) rather than lat/long geodesy, since a
# synthetic grid has no real-world coordinates at all.
_resolvesyntheticactive(active::Nothing, dim, cellsize) = fill(true, dim)

_resolvesyntheticactive(active::AbstractMatrix{Bool}, dim, cellsize) = active

function _resolvesyntheticactive(active::CircleMaskSpec, dim, cellsize)
    isnothing(active.centre) ||
        error("An explicit `centre` for `CircleMaskSpec` needs real-world " *
              "coordinates and is only supported for data-driven (geographic) environments; " *
              "a synthetic grid has none, so only the default `centre = nothing` (the grid's " *
              "own physical centre) is supported here.")
    ny, nx = dim
    cy, cx = (ny + 1) / 2, (nx + 1) / 2
    mask = falses(ny, nx)
    for y in 1:ny, x in 1:nx
        mask[y, x] = hypot((y - cy) * cellsize, (x - cx) * cellsize) <=
                     active.radius
    end
    return mask
end

function _resolvesyntheticactive(active, dim, cellsize)
    return error("unrecognised `active` argument of type $(typeof(active)) for a synthetic " *
                 "environment; use nothing, a Matrix{Bool}, or CircleMaskSpec(...).")
end

# A per-cell field of per-area supply-rate values over `dim` — one method per recognised
# `AbstractSyntheticSpec` kind (a bare quantity is uniform; `UniformSpec` unwraps its value;
# `GradientSpec` and `PeakedSpec` reuse `_gradientfield`'s spatial pattern, mirroring how those
# specs materialise a regime but without any axis tagging — a supply is a `Resource`-role rate, not
# a `Condition`-role niche value). Reached through `_specfield`, which adds the one kind this family
# cannot answer for itself — `NicheSpec`, whose values are class codes — so every synthetic layer,
# either role and either kind of study area, is generated here.
function _syntheticsupplyfield(maxsupply, dim, rowsincreasenorth::Bool = true)
    return fill(float(maxsupply),
                dim)
end

function _syntheticsupplyfield(spec::UniformSpec, dim,
                               rowsincreasenorth::Bool = true)
    return fill(float(spec.value), dim)
end

function _syntheticsupplyfield(spec::GradientSpec, dim,
                               rowsincreasenorth::Bool = true)
    lo, hi = float(spec.low), float(spec.high)
    return lo .+
           (hi - lo) .* _gradientfield(dim, spec.orientation, rowsincreasenorth)
end

function _syntheticsupplyfield(spec::PeakedSpec, dim,
                               rowsincreasenorth::Bool = true)
    lo, hi = float(spec.outside), float(spec.inside)
    peaked = 1 .-
             abs.(2 .*
                  _gradientfield(dim, spec.orientation, rowsincreasenorth) .- 1)
    return lo .+ (hi - lo) .* peaked
end
