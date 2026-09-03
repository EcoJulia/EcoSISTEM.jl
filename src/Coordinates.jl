# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The coordinate vocabulary the study area, the climate readers and the mask specs all state
# positions and sizes in.
#
# Two kinds of two-dimensional quantity live here and are deliberately different types. A *position*
# says where something is; a *size* says how far apart two things are. They are the point and the
# vector of an affine space: subtracting two positions gives a size, adding a size to a position
# gives a position, and adding two positions is meaningless and therefore undefined. Keeping them
# apart in the type system is what stops a cell's extent being passed where a location is wanted.

using Unitful

using Unitful.DefaultSymbols

import GeoInterface

import Extents

using CSV

"""
    SpatialKind

Supertype of the tags that say what kind of two-dimensional quantity a [`Spatial2D`](@ref) is.

The two are `AbsolutePosition` (a place) and `RelativeOffset` (a separation between places). Both
are internal: a caller names [`SpatialLocation`](@ref) or [`SpatialSize`](@ref) instead, and the tag
appears only inside the type.
"""
abstract type SpatialKind end

# A place: somewhere a thing can be. Validated where the units make a valid range statable.
struct AbsolutePosition <: SpatialKind end

# A separation between places: an extent, a cell size, a step. Unbounded, and signed, because a
# displacement may legitimately run in either direction.
struct RelativeOffset <: SpatialKind end

"""
    Spatial2D{Kind, T}

A two-dimensional quantity in `(y, x)` order, with `Kind` saying whether it is a place or a
separation and `T` the type of each component.

Write [`SpatialLocation`](@ref) or [`SpatialSize`](@ref) rather than this; both are aliases that pin
`Kind`, so a method annotated with either refuses the other.

The order is `y` then `x` -- north-south before east-west -- matching every other pair in the
package, and the field names are always `y`/`x` whatever the units. A degree quantity makes a
geographic pair and a length quantity a projected one, so the units say which frame a value is in
without a second type having to.

# Fields

  - `y`: the north-south component -- a latitude for a geographic position, a northing otherwise.
  - `x`: the east-west component -- a longitude for a geographic position, an easting otherwise.
"""
struct Spatial2D{Kind <: SpatialKind, T}
    y::T
    x::T
    function Spatial2D{Kind, T}(y::T, x::T) where {Kind <: SpatialKind, T}
        _checkspatial(Kind, y, x)
        return new{Kind, T}(y, x)
    end
end

function Spatial2D{Kind}(y::T, x::T) where {Kind <: SpatialKind, T}
    return Spatial2D{Kind, T}(y, x)
end

"""
    SpatialLocation(y, x)

A place, in `(y, x)` order: `SpatialLocation(50.0°, -3.0°)` geographically, or
`SpatialLocation(5000.0m, 3000.0m)` on a projected grid.

Degrees are checked against the valid geographic range on construction; other units are not, because
what counts as in range is a property of the coordinate reference system and is not known here.
"""
const SpatialLocation{T} = Spatial2D{AbsolutePosition, T}

"""
    SpatialSize(y, x)

A separation in `(y, x)` order -- the extent of a cell, the size of a grid, the step between two
coordinates: `SpatialSize(1.0km, 1.0km)`.

Never range-checked, and signed: a displacement may run either way, and an angular separation may
exceed the range a latitude is confined to.
"""
const SpatialSize{T} = Spatial2D{RelativeOffset, T}

"""
    LatLong(lat, long)
    LatLong(; lat, long)

A geographic **point** in Unitful degrees, validated on construction: `LatLong(50.0°, -3.0°)`.
Latitude must lie in `[-90°, 90°]` and longitude in `(-180°, 180°]`.

A **region** is an `Extents.Extent` rather than one of these - `Extent(Y = (54.6°, 58.7°),
X = (-6.2°, -1.8°))`. That is the bounding-box vocabulary the whole Julia geo ecosystem shares, and
reusing it beats a parallel type of our own.

Reach the components with [`getlat`](@ref) and [`getlong`](@ref). It implements the GeoInterface
`PointTrait`, with the usual `X = longitude`, `Y = latitude` order, stripped to plain degrees.
"""
const LatLong = SpatialLocation{typeof(1.0°)}

LatLong(; lat, long) = LatLong(lat, long)

function Base.show(io::IO, p::SpatialLocation{typeof(1.0°)})
    return print(io, "SpatialLocation(", _compass(p.y, 'N', 'S'), ", ",
                 _compass(p.x, 'E', 'W'), ")")
end

function Base.show(io::IO, p::SpatialLocation)
    return print(io, "SpatialLocation(y = ", p.y, ", x = ", p.x, ")")
end

function Base.show(io::IO, d::SpatialSize)
    return print(io, "SpatialSize(Δy = ", d.y, ", Δx = ", d.x, ")")
end

"""
    getlat(p)
    getlong(p)

Return the latitude or longitude of a geographic point.

# Arguments

  - `p`: the point to read, which must be a [`LatLong`](@ref).

Defined only for a [`LatLong`](@ref), so asking a *projected* location for its latitude is a
`MethodError` rather than a plausible-looking northing: a projected `y` is not a latitude, and the
type system is where that should be caught.
"""
getlat(p::LatLong) = p.y

@doc (@doc getlat) getlong(p::LatLong) = p.x

# --- Named regions -------------------------------------------------------------------------------
# A region's geographic extent, looked up by name from the shipped `data/NaturalEarth/regions.csv`.
# Here rather than with the readers because it is geographic vocabulary, not climate data: it names
# no dataset and reads no raster, and what it returns is the `Extents.Extent` this file deals in.

"""
    boundingbox(region::AbstractString; level = nothing, coverage = LargestLandmass(),
                round = false)

Return a named region's geographic bounding box, as an `Extents.Extent` of `°` intervals ready to
pass as the `cut` keyword to [`read`](@ref) and the other raster readers.

The answer costs no download: it is read from a table generated from Natural Earth's own polygons at
1:10m, so the box agrees with the shape the same name gives. Rounding onto a source's own lattice is
what lets a layer be **aggregated** exactly rather than resampled: WorldClim's cells are 10
arcminutes, EarthEnv's and CHELSA's 30 arcseconds.

# Arguments

  - `region`: the region's name, matched case-insensitively.
  - `level`: which kind of region the name means - `"ADMIN"` for a country, `"SUBUNIT"` for a
    constituent country, `"Physical Island"` for a landmass, and so on;
    `EcoSISTEM.NATURALEARTH_LEVELS` lists them all. Optional, and only needed where a name means
    genuinely different ground at different levels: "Scotland" is the same box whether asked for as
    a map unit or a map subunit, while "Africa" as a continent stops 54 degrees west of "Africa" as
    a UN region, and naming a level is then required rather than guessed at.
  - `coverage`: how much of what the name covers to take - [`LargestLandmass`](@ref), the default,
    for the principal landmass alone, or [`AllTerritories`](@ref) for every scattered territory too.
  - `round`: an angular step to snap the box **outwards** onto, so the result fully contains the
    exact box. Any angular unit will do: `round = 5°`, `round = 10arcminute`,
    `round = 30arcsecond`. The result is in degrees whichever was used. `false`, the default, leaves
    it unrounded.
"""
function boundingbox(region::AbstractString; level = nothing,
                     coverage::AbstractCoverage = LargestLandmass(),
                     round = false)
    lvl = isnothing(level) ? _resolvelevel(region, coverage) : String(level)
    row = _regionrow(lvl, region)
    isnothing(row) &&
        error("No region named \"$region\" at level \"$lvl\". " *
              "`EcoSISTEM.NATURALEARTH_LEVELS` lists the levels.")
    box = _coveragebox(row, coverage, region, lvl)
    # An `Extents.Extent` states an interval, and an interval cannot run the long way round: a
    # wrapping region written into one would read as empty, or as its own complement. Refusing is
    # the only honest answer, and the two halves are what a caller needs to proceed.
    box.wraps &&
        error("\"$region\" ($lvl) crosses the antimeridian, from $(box.west)° east to " *
              "$(box.east)°, so it has no bounding box that is a single interval of longitude. " *
              "Read the two halves separately ($(box.west)° to 180° and -180° to $(box.east)°), " *
              "or take `coverage = LargestLandmass()`, which cannot wrap.")
    south, north = box.south * °, box.north * °
    west, east = box.west * °, box.east * °
    if round !== false
        south, west = _snapout(floor, south, round),
                      _snapout(floor, west, round)
        north, east = _snapout(ceil, north, round), _snapout(ceil, east, round)
    end
    return Extents.Extent(Y = (south, north), X = (west, east))
end

# --- The affine algebra ------------------------------------------------------------------------
#
# Two positions give a separation; a position plus a separation is a position; separations add,
# subtract and scale among themselves. Adding two positions is deliberately absent: there is no
# meaningful answer, and a `MethodError` at the call site says so better than a wrong number would.
function Base.:-(a::SpatialLocation, b::SpatialLocation)
    return SpatialSize(a.y - b.y, a.x - b.x)
end

function Base.:-(a::SpatialLocation, d::SpatialSize)
    return SpatialLocation(a.y - d.y, a.x - d.x)
end

Base.:-(a::SpatialSize, b::SpatialSize) = SpatialSize(a.y - b.y, a.x - b.x)

Base.:-(d::SpatialSize) = SpatialSize(-d.y, -d.x)

function Base.:+(a::SpatialLocation, d::SpatialSize)
    return SpatialLocation(a.y + d.y, a.x + d.x)
end

Base.:+(d::SpatialSize, a::SpatialLocation) = a + d

Base.:+(a::SpatialSize, b::SpatialSize) = SpatialSize(a.y + b.y, a.x + b.x)

Base.:*(d::SpatialSize, k::Number) = SpatialSize(d.y * k, d.x * k)

Base.:*(k::Number, d::SpatialSize) = d * k

Base.:/(d::SpatialSize, k::Number) = SpatialSize(d.y / k, d.x / k)

# --- GeoInterface: a point is a PointTrait geometry. GeoInterface fixes the coordinate order as
# (X, Y); by convention X = longitude, Y = latitude, so coord 1 = long and coord 2 = lat - the
# only place the lat/long-vs-X/Y reversal lives. Coordinates are stripped to plain degrees. ---
GeoInterface.isgeometry(::Type{LatLong}) = true

function GeoInterface.geomtrait(::LatLong)
    return GeoInterface.PointTrait()
end

function GeoInterface.ncoord(::GeoInterface.PointTrait,
                             ::LatLong)
    return 2
end

function GeoInterface.getcoord(::GeoInterface.PointTrait,
                               p::LatLong, i::Integer)
    return i == 1 ? ustrip(°, getlong(p)) : ustrip(°, getlat(p))
end

# The construction guard, dispatched on the tag and the component type together. A geographic
# *position* is the one case with a range that can be stated without knowing the CRS; everything
# else -- projected positions, and every separation -- is left alone.
_checkspatial(::Type{<:SpatialKind}, y, x) = nothing

function _checkspatial(::Type{AbsolutePosition}, lat::typeof(1.0°),
                       long::typeof(1.0°))
    return _checkcoords(lat, long)
end

# == Functions ==================================================================================

# --- Display -----------------------------------------------------------------------------------
#
# There are four things a `Spatial2D` can be -- a place or a separation, geographic or projected --
# and the default `show` distinguishes none of them readably, printing the full parameterised type
# ahead of the values. These make each case identifiable at a glance:
#
#   SpatialLocation(50.0°N, 3.0°W)          a geographic place
#   SpatialLocation(y = 5000.0 m, x = 3000.0 m)   a projected place
#   SpatialSize(Δy = 0.5°, Δx = 0.5°)       an angular separation
#   SpatialSize(Δy = 1.0 km, Δx = 1.0 km)   a projected separation
#
# The alias name gives the kind, the units give the frame, and `Δ` marks a separation. Compass
# letters are used only for a geographic *place*, where a hemisphere is meaningful; a separation has
# none, which is why it keeps signed `Δy`/`Δx`.

# One component of a geographic place, written with its hemisphere letter so a southern latitude
# cannot be misread as a northern one at a glance.
function _compass(v, positive::Char, negative::Char)
    return string(abs(v), v < zero(v) ? negative : positive)
end

# Error unless (`lat`, `long`) fall within valid geographic bounds: latitude in [-90°, 90°],
# longitude in (-180°, 180°].
function _checkcoords(lat::typeof(1.0°), long::typeof(1.0°))
    -90.0° ≤ lat ≤ 90.0° || error("Latitude coordinate is out of bounds")
    -180.0° < long ≤ 180.0° || error("Longitude coordinate is out of bounds")
    return nothing
end

# Refuse a geographic region that could not describe ground: reversed either way, or off the globe.
# An `Extents.Extent` belongs to `Extents`, so there is no construction hook to hang this on and it
# runs at the point of use instead.
function _checkgeographicextent(e::Extents.Extent)
    y, x = e.Y, e.X
    y[1] ≤ y[2] ||
        error("Latitude interval is reversed (south > north): got $(y[1]) ... $(y[2]).")
    x[1] ≤ x[2] ||
        error("Longitude interval is reversed (west > east): got $(x[1]) ... $(x[2]). This is how " *
              "an area crossing the antimeridian (±180°) would be expressed, but EcoSISTEM does " *
              "not support dateline-crossing regions - split it into two, one either side.")
    -90° ≤ y[1] && y[2] ≤ 90° ||
        error("Latitude interval is out of bounds: got $(y[1]) ... $(y[2]).")
    -180° < x[1] && x[2] ≤ 180° ||
        error("Longitude interval is out of bounds: got $(x[1]) ... $(x[2]).")
    return nothing
end

# Snap a `°` coordinate `x` to the nearest multiple of step `r` (any angular quantity) in direction
# `dirn` (`floor` = down, `ceil` = up). Applied as floor to the low edge and ceil to the high edge,
# this always rounds *outward* so a rounded box encloses the exact one.
#
# Converted back to `°` rather than left in `r`'s unit: the step may legitimately be given in
# arcminutes or arcseconds, but `boundingbox` promises degrees, and `3270.0 ′` - correct, but not
# degrees - would break that for every caller downstream.
_snapout(dirn, x, r) = uconvert(°, dirn(uconvert(NoUnits, x / r)) * r)
