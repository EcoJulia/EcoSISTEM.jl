# SPDX-License-Identifier: LGPL-3.0-or-later

module TestExtractClimate

using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
using RasterDataSources
using Test

# A non-square-grid bioclim raster (`:var` third axis): latitude step 20°, longitude step 10° — the
# differing steps guard the separate lat/long half-cell fix.
blat = collect(0.0:20.0:80.0) .* °
blong = collect(0.0:10.0:40.0) .* °
bc = ClimateRaster(WorldClim{BioClim},
                   AxisArray(Float64.(reshape(1:(length(blat) * length(blong) * 3),
                                              length(blat), length(blong), 3)),
                             Axis{:latitude}(blat),
                             Axis{:longitude}(blong),
                             Axis{:var}(1:3)))

# A worldclim monthly raster (`:time` third axis = (1:12)·month).
wlat = collect(0.0:10.0:40.0) .* °
wlong = collect(0.0:10.0:40.0) .* °
wc = ClimateRaster(WorldClim{Climate},
                   AxisArray(Float64.(reshape(1:(length(wlat) * length(wlong) * 12),
                                              length(wlat), length(wlong), 12)),
                             Axis{:latitude}(wlat),
                             Axis{:longitude}(wlong),
                             Axis{:time}((1:12) .* month)))

# An ERA raster with an absolute calendar time axis: Jan 2000 … Mar 2001 (12 + 3 = 15 months).
# Longitudes run out to 120° (> 90°) to exercise the coordinate bounds. 12month == year, so
# month index `24000` is January 2000.
elat = collect(0.0:10.0:40.0) .* °
elong = collect(0.0:30.0:120.0) .* °
etime = (24000:24014) .* month
era = ERA(AxisArray(Float64.(reshape(1:(length(elat) * length(elong) * 15),
                                     length(elat), length(elong), 15)) .* K,
                    Axis{:latitude}(elat),
                    Axis{:longitude}(elong),
                    Axis{:time}(etime)))

ref = create_reference(1.0)

@testset "separate lat/long steps (non-square grid)" begin
    # long step (10°) differs from lat step (20°): using the lat step for longitude would widen the
    # cell to [0°, 20°] and pick the wrong longitude. The correct cell is isolated.
    @test extractvalues(20.0°, 10.0°, bc, 2) ==
          bc.array[atvalue(20.0°), atvalue(10.0°), 2]
end

@testset "coordinate bounds" begin
    @test_throws ErrorException extractvalues(120.0°, 10.0°, bc, 1)      # lat > 90°
    @test_throws ErrorException extractvalues(10.0°, 200.0°, bc, 1)      # long > 180°
    @test_throws ErrorException extractvalues([120.0°], [10.0°], bc, 1)
    # a longitude past ±90° is fine (guards the earlier checkbounds swap)
    @test extractvalues(10.0°, 120.0°, era, 3) ==
          era.array[atvalue(10.0°), atvalue(120.0°), 3]
end

@testset "shape grid: scalar/vector × single/multi/all" begin
    # scalar point, single slice -> scalar
    @test extractvalues(10.0°, 10.0°, wc, 3month) ==
          wc.array[atvalue(10.0°), atvalue(10.0°), 3month]
    # scalar point, several slices -> length-M vector
    v = extractvalues(10.0°, 10.0°, wc, [3month, 6month])
    @test v isa AbstractVector && length(v) == 2
    # scalar point, all slices -> full third axis
    @test length(extractvalues(10.0°, 10.0°, wc)) == 12
    # vector points, single slice -> N×1
    m1 = extractvalues([10.0°, 20.0°], [10.0°, 20.0°], wc, 3month)
    @test size(m1) == (2, 1)
    # vector points, all slices -> N×M
    @test size(extractvalues([10.0°, 20.0°], [10.0°, 20.0°], wc)) == (2, 12)
    # all 3 bioclim variables at N points
    @test size(extractvalues([20.0°, 40.0°], [10.0°, 20.0°], bc)) == (2, 3)
    # mismatched vector lengths error
    @test_throws ErrorException extractvalues([10.0°], [10.0°, 20.0°], wc)
end

@testset "year selection on an absolute time axis" begin
    # a full year present in the data -> its 12 entries
    @test length(extractvalues(10.0°, 10.0°, era; year = 2000)) == 12
    # a partially-present year -> only the entries that exist (Jan–Mar 2001)
    @test length(extractvalues(10.0°, 10.0°, era; year = 2001)) == 3
    # a year outside the data -> empty, no error/NaN
    @test isempty(extractvalues(10.0°, 10.0°, era; year = 1999))
    # a month subset within a year
    @test length(extractvalues(10.0°, 10.0°, era, [6, 7, 8]; year = 2000)) == 3
    # requested months that don't exist are simply absent
    @test length(extractvalues(10.0°, 10.0°, era, [2, 3]; year = 2001)) == 2
    # vector points × a partial year -> N × (count)
    @test size(extractvalues([10.0°, 20.0°], [10.0°, 20.0°], era; year = 2001)) ==
          (2, 3)
    # `year` on a dataset without a time axis is an error
    @test_throws ErrorException extractvalues(10.0°, 10.0°, bc; year = 2000)
end

@testset "Reference (2-D) and its lat/long convention" begin
    a = ref.array
    @test AxisArrays.axisnames(a) == (:latitude, :longitude)
    @test extrema(AxisArrays.axisvalues(a)[1]) == (-90.0°, 90.0°)
    @test extrema(AxisArrays.axisvalues(a)[2]) == (-180.0°, 180.0°)
    # scalar extraction indexes lat -> dim 1, long -> dim 2 (longitude 120° round-trips)
    @test extractvalues(30.0°, 120.0°, ref) ==
          a[atvalue(30.0°), atvalue(120.0°)]
    @test_throws ErrorException extractvalues(120.0°, 30.0°, ref)   # lat > 90°
end

@testset "LatLong location forms (the core methods)" begin
    # a single LatLong location matches the separate-scalar form
    @test extractvalues(LatLong(10.0°, 10.0°), wc, 3month) ==
          extractvalues(10.0°, 10.0°, wc, 3month)
    # a vector of locations matches the separate-vectors form
    pts = [LatLong(10.0°, 10.0°), LatLong(20.0°, 20.0°)]
    @test extractvalues(pts, wc, 3month) ==
          extractvalues([10.0°, 20.0°], [10.0°, 20.0°], wc, 3month)
    @test size(extractvalues(pts, wc)) == (2, 12)
end

end
