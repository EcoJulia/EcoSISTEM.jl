# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using IntervalSets: ClosedInterval, leftendpoint, rightendpoint
import GeoInterface
import Extents

# Error unless (`lat`, `long`) fall within valid geographic bounds: latitude in [-90°, 90°],
# longitude in (-180°, 180°].
function _checkcoords(lat::typeof(1.0°), long::typeof(1.0°))
    -90.0° ≤ lat ≤ 90.0° || error("Latitude coordinate is out of bounds")
    -180.0° < long ≤ 180.0° || error("Longitude coordinate is out of bounds")
    return nothing
end

# Validate a `LatLong`'s arguments: a scalar point delegates to `_checkcoords`; a pair of intervals
# checks each endpoint is in bounds and that neither interval is reversed.
_checklatlong(lat::typeof(1.0°), long::typeof(1.0°)) = _checkcoords(lat, long)
function _checklatlong(lat::ClosedInterval, long::ClosedInterval)
    leftendpoint(lat) ≤ rightendpoint(lat) ||
        error("Latitude interval is reversed (south > north)")
    leftendpoint(long) ≤ rightendpoint(long) ||
        error("Longitude interval is reversed (west > east). This is how an area crossing the " *
              "antimeridian (±180°) would be expressed, but EcoSISTEM does not support " *
              "dateline-crossing regions — split it into two boxes, one either side of ±180°.")
    -90° ≤ leftendpoint(lat) && rightendpoint(lat) ≤ 90° ||
        error("Latitude interval is out of bounds")
    -180° < leftendpoint(long) && rightendpoint(long) ≤ 180° ||
        error("Longitude interval is out of bounds")
    return nothing
end

"""
    LatLong(lat, long)
    LatLong(; lat, long)

A geographic location or region in Unitful degrees, validated on construction. Pass two angles for a
point (`LatLong(50.0°, -3.0°)`), or two `°` intervals for a bounding box
(`LatLong(54.6°..58.7°, -6.2°.. -1.8°)`). Latitude must lie in `[-90°, 90°]` and longitude in
`(-180°, 180°]` — for a box, both endpoints, and neither interval may be reversed.

The fields `.lat` and `.long` hold the latitude and longitude (each a scalar or an interval), so a
`LatLong` is a drop-in replacement for the old `(lat = …, long = …)` coordinate/`cut`/`region`
NamedTuples. A point implements the GeoInterface `PointTrait` (with the usual `X = longitude`,
`Y = latitude` order, stripped to plain degrees); a box implements `Extents.extent`.
"""
struct LatLong{T}
    lat::T
    long::T
    function LatLong(lat::T, long::T) where {T}
        _checklatlong(lat, long)
        return new{T}(lat, long)
    end
end
LatLong(; lat, long) = LatLong(lat, long)

# --- GeoInterface: a point is a PointTrait geometry. GeoInterface fixes the coordinate order as
# (X, Y); by convention X = longitude, Y = latitude, so coord 1 = long and coord 2 = lat — the
# only place the lat/long-vs-X/Y reversal lives. Coordinates are stripped to plain degrees. ---
GeoInterface.isgeometry(::Type{<:LatLong{<:Unitful.Quantity}}) = true
function GeoInterface.geomtrait(::LatLong{<:Unitful.Quantity})
    return GeoInterface.PointTrait()
end
function GeoInterface.ncoord(::GeoInterface.PointTrait,
                             ::LatLong{<:Unitful.Quantity})
    return 2
end
function GeoInterface.getcoord(::GeoInterface.PointTrait,
                               p::LatLong{<:Unitful.Quantity}, i::Integer)
    return i == 1 ? ustrip(°, p.long) : ustrip(°, p.lat)
end

# --- Extents: a box exposes its bounding extent (plain °) — the Rasters/DimensionalData currency. ---
function Extents.extent(r::LatLong{<:ClosedInterval})
    return Extents.Extent(X = (ustrip(°, leftendpoint(r.long)),
                               ustrip(°, rightendpoint(r.long))),
                          Y = (ustrip(°, leftendpoint(r.lat)),
                               ustrip(°, rightendpoint(r.lat))))
end
