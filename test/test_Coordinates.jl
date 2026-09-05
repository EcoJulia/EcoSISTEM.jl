# SPDX-License-Identifier: LGPL-3.0-or-later

module TestCoordinates

using EcoSISTEM
# `[C7-VIS]` C: the unit-general spatial aliases are `public`, not exported.
using EcoSISTEM: SpatialLocation, SpatialSize, getlat, getlong
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
import GeoInterface
import Extents
using Test

@testset "LatLong is a geographic SpatialLocation" begin
    p = LatLong(50.0°, -3.0°)
    # `LatLong` is now an alias, not a type of its own: a `SpatialLocation` whose components are
    # pinned to degrees. So it *is* a `Spatial2D`, and everything written against those works on it.
    @test p isa SpatialLocation
    @test p isa EcoSISTEM.Spatial2D
    @test !(SpatialLocation(5.0km, 3.0km) isa LatLong)

    # Fields are `y`/`x` like every other pair in the package; `lat`/`long` are accessors.
    @test p.y == 50.0° && p.x == -3.0°
    @test getlat(p) == 50.0° && getlong(p) == -3.0°
    @test LatLong(lat = 50.0°, long = -3.0°) == p

    # Asking a *projected* location for its latitude is a `MethodError`, not a plausible-looking
    # northing. That is the whole reason the accessors exist rather than the fields being named
    # `lat`/`long`: a projected `y` is not a latitude.
    @test_throws MethodError getlat(SpatialLocation(5.0km, 3.0km))

    # out-of-bounds points error, as before
    @test_throws ErrorException LatLong(100.0°, 0.0°)     # latitude > 90°
    @test_throws ErrorException LatLong(0.0°, 200.0°)     # longitude > 180°

    # A *separation* is unbounded, and the same numbers are legal there: 200° of longitude is a
    # perfectly good displacement, where 200° of longitude is not a place.
    @test SpatialSize(200.0°, 200.0°) isa SpatialSize
end

@testset "a region is an Extents.Extent, not a LatLong" begin
    # A `LatLong` must be a point and never also a box, told apart by whether its fields hold
    # scalars or intervals. One name for two kinds of thing is what `Spatial2D` exists to stop, and
    # a bounding box already has a shared vocabulary - so a region is an `Extents.Extent`.
    b = Extents.Extent(Y = (54.6°, 58.7°), X = (-6.2°, -1.8°))
    @test b.Y == (54.6°, 58.7°) && b.X == (-6.2°, -1.8°)
    # The interval form is not merely rejected - it is unconstructible: `LatLong` is pinned to
    # `SpatialLocation{typeof(1.0°)}`, so an interval cannot inhabit it. Asserting that needed a
    # whole `IntervalSets` dependency for one marginal check, so it is left to the alias to say.

    # `Extents.Extent` is Extents' type, so it has no construction hook for a guard. The check
    # a box could make on construction is made where a region is *used* instead - the error
    # survives, only its timing moves.
    ok = EcoSISTEM._checkgeographicextent
    @test isnothing(ok(b))
    @test_throws ErrorException ok(Extents.Extent(Y = (0.0°, 100.0°),
                                                  X = (0.0°, 1.0°)))
    @test_throws ErrorException ok(Extents.Extent(Y = (58.7°, 54.6°),
                                                  X = (-6.2°, -1.8°)))
    # a region crossing the antimeridian (west > east) is still not supported
    @test_throws ErrorException ok(Extents.Extent(Y = (0.0°, 1.0°),
                                                  X = (170.0°, -170.0°)))
end

@testset "LatLong point -> GeoInterface PointTrait" begin
    p = LatLong(50.0°, -3.0°)
    @test GeoInterface.isgeometry(p)
    @test GeoInterface.geomtrait(p) isa GeoInterface.PointTrait
    @test GeoInterface.ncoord(p) == 2
    # GeoInterface order is (X, Y) = (longitude, latitude), in plain degrees
    @test GeoInterface.getcoord(p, 1) == -3.0
    @test GeoInterface.getcoord(p, 2) == 50.0
    @test GeoInterface.x(p) == -3.0
    @test GeoInterface.y(p) == 50.0
end

# `boundingbox` only reads the shipped table, so it runs on every platform and needs no network.
# It returns an `Extents.Extent`, so its components are `.Y`/`.X` - the shared geo vocabulary -
# rather than a `LatLong` box's old `.lat`/`.long`.
@testset "Bounding boxes" begin
    # The default coverage is everything the name covers, which is what the source means by it: for
    # Scotland that reaches Rockall in the west and Shetland in the north.
    isl = EcoSISTEM.boundingbox("Scotland")
    @test maximum(isl.Y) ≈ 60.85° atol=0.01°
    @test minimum(isl.X) ≈ -13.69° atol=0.01°
    @test EcoSISTEM.boundingbox("Scotland", coverage = AllTerritories()) == isl

    # `LargestLandmass` narrows it to the principal landmass - the Scottish mainland.
    scot = EcoSISTEM.boundingbox("Scotland", coverage = LargestLandmass())
    @test minimum(scot.Y) ≈ 54.63° atol=0.01°
    @test maximum(scot.Y) ≈ 58.68° atol=0.01°
    @test minimum(scot.X) ≈ -6.24° atol=0.01°
    @test maximum(scot.X) ≈ -1.76° atol=0.01°
    # Every coverage of a name is contained by `AllTerritories`, which is what "all" has to mean.
    @test minimum(isl.X) ≤ minimum(scot.X) && maximum(isl.X) ≥ maximum(scot.X)
    @test minimum(isl.Y) ≤ minimum(scot.Y) && maximum(isl.Y) ≥ maximum(scot.Y)

    # Asking for more landmasses widens the box monotonically, and never past `AllTerritories`.
    two = EcoSISTEM.boundingbox("Scotland",
                                coverage = LargestLandmass(count = 2))
    @test minimum(two.X) ≤ minimum(scot.X) && minimum(two.X) ≥ minimum(isl.X)

    # round snaps outwards to the nearest multiple, enclosing the exact box
    rnd = EcoSISTEM.boundingbox("Scotland", coverage = LargestLandmass(),
                                round = 5°)
    @test minimum(rnd.Y) == 50° && maximum(rnd.Y) == 60°
    @test minimum(rnd.X) == -10° && maximum(rnd.X) == 0°
    @test minimum(rnd.Y) ≤ minimum(scot.Y) &&
          maximum(rnd.Y) ≥ maximum(scot.Y)
    @test minimum(rnd.X) ≤ minimum(scot.X) &&
          maximum(rnd.X) ≥ maximum(scot.X)
    @test_throws ErrorException EcoSISTEM.boundingbox("Atlantis")

    # A name meaning different ground at different levels is refused rather than guessed at: as a
    # continent "Africa" is 55 countries, as a UN region 62, and the two boxes disagree.
    @test_throws ErrorException EcoSISTEM.boundingbox("Africa")
    afr = EcoSISTEM.boundingbox("Africa", level = "CONTINENT",
                                coverage = LargestLandmass())
    @test minimum(afr.X) ≈ -17.54° atol=0.01°
    @test maximum(afr.X) ≈ 51.42° atol=0.01°
    # A name that is unambiguous *in effect* - the same box at every level it appears at - needs no
    # level, which is why "Scotland" above did not take one.
    @test EcoSISTEM.boundingbox("Scotland", level = "SUBUNIT") == isl
    @test EcoSISTEM.boundingbox("Scotland", level = "GEOUNIT") == isl

    # An `Extents.Extent` holds an interval, so a region spanning the date line has no box it could
    # be given: refused, rather than handed back inside out.
    @test_throws ErrorException EcoSISTEM.boundingbox("Russia", level = "ADMIN",
                                                      coverage = AllTerritories())
    # ...but its principal landmass is a single component, which cannot cross the line.
    @test maximum(EcoSISTEM.boundingbox("Russia", level = "ADMIN",
                                        coverage = LargestLandmass()).X) ==
          180.0°

    # `round` takes any angular step, not just whole degrees - snapping to a source's own lattice
    # (WorldClim 10 arcmin, EarthEnv/CHELSA 30 arcsec) is what lets a cut layer be aggregated exactly
    # rather than resampled. The result must still come back in `°` whichever unit the step used:
    # `_snapout` returns `n * r`, so before it converted back this handed out `3270.0 ′`, correct but
    # not what the docstring promises and not what anything downstream expects.
    arcmin = EcoSISTEM.boundingbox("Scotland", coverage = LargestLandmass(),
                                   round = 30arcminute)
    @test unit(minimum(arcmin.Y)) === °
    @test minimum(arcmin.Y) == 54.5° && maximum(arcmin.Y) == 59.0°
    @test minimum(arcmin.Y) ≤ minimum(scot.Y)     # still encloses the exact box
    @test maximum(arcmin.Y) ≥ maximum(scot.Y)
    # ...and a degree step expressed as arcminutes agrees exactly with the degree form.
    @test EcoSISTEM.boundingbox("Scotland", coverage = LargestLandmass(),
                                round = 300arcminute).Y == rnd.Y
end

end
