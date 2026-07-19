# SPDX-License-Identifier: LGPL-3.0-or-later

module TestCoordinates

using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
using IntervalSets
import GeoInterface
import Extents
using Test

@testset "LatLong construction & validation" begin
    p = LatLong(50.0°, -3.0°)
    @test p.lat == 50.0° && p.long == -3.0°
    q = LatLong(; lat = 50.0°, long = -3.0°)
    @test q.lat == 50.0° && q.long == -3.0°

    b = LatLong(54.6° .. 58.7°, -6.2° .. -1.8°)
    @test b.lat == 54.6° .. 58.7° && b.long == -6.2° .. -1.8°

    # out-of-bounds points error
    @test_throws ErrorException LatLong(100.0°, 0.0°)     # latitude > 90°
    @test_throws ErrorException LatLong(0.0°, 200.0°)     # longitude > 180°
    # out-of-bounds / reversed boxes error
    @test_throws ErrorException LatLong(0.0° .. 100.0°, 0.0° .. 1.0°)    # lat endpoint > 90°
    @test_throws ErrorException LatLong(58.7° .. 54.6°, -6.2° .. -1.8°)  # reversed latitude
    # a box crossing the antimeridian (west > east) is not supported
    @test_throws ErrorException LatLong(0.0° .. 1.0°, 170.0° .. -170.0°)
end

@testset "LatLong point → GeoInterface PointTrait" begin
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

@testset "LatLong box → extent" begin
    b = LatLong(54.6° .. 58.7°, -6.2° .. -1.8°)
    @test Extents.extent(b) ==
          Extents.Extent(X = (-6.2, -1.8), Y = (54.6, 58.7))
end

end
