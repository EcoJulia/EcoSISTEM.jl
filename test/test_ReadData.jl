# SPDX-License-Identifier: LGPL-3.0-or-later

module TestReadData

using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources
using Rasters
using Statistics
using Test

if !Sys.iswindows()
    # `getraster` returns the full path(s) to the downloaded file(s), so use those directly rather
    # than reconstructing RasterDataSources' folder layout. Pre-fetching here (outside the
    # `@test_nowarn`s) also keeps download messages out of those tests on an empty cache.
    bio1 = getraster(WorldClim{BioClim}, :bio1)               # one tif path
    wind = getraster(WorldClim{Climate}, :wind, month = 1:12) # 12 monthly tif paths
    getraster(EarthEnv{LandCover})
    getraster(CHELSA{BioClim}, 1)
    # A directory holding exactly the 12 downloaded wind tifs, to exercise the directory readers
    # `readCRUTS`/`readCHELSA_monthly` (the variable name only fixes the unit that gets attached).
    winddir = dirname(first(wind))

    @testset "Reading functions" begin
        @test_nowarn read(WorldClim{Climate}, :wind, month = 1:12)
        @test_nowarn readCRUTS(winddir, "tavg")
        @test_nowarn readCHELSA_monthly(winddir, "wind")
        # CHELSA bioclim is a 43200×20880 global grid; reading it at full
        # resolution allocates several ~7 GiB Float64 arrays and OOMs CI.
        # Downsample to WorldClim's 10-arcmin resolution to keep it bounded.
        @test_nowarn read(CHELSA{BioClim}, 1, scale = 20)
        @test_nowarn read(EarthEnv{LandCover})
        @test_nowarn readfile(bio1)
    end

    @testset "Output data" begin
        bioclim = read(WorldClim{BioClim})
        cr = readCRUTS(winddir, "tavg")
        ch_b = read(CHELSA{BioClim}, 1, scale = 20)
        rf = readfile(bio1)

        @test unit(bioclim.array[1]) == unit(rf[1]) == unit(ch_b.array[1]) ==
              NoUnits
    end

    @testset "Output data 2" begin
        landcover = read(EarthEnv{LandCover})
        @test unit(landcover.array[1]) == NoUnits
    end

    @testset "Output data 3" begin
        cr = readCRUTS(winddir, "tavg")
        worldclim = read(WorldClim{Climate}, :wind)
        ch_m = readCHELSA_monthly(winddir, "wind")

        # the directory readers self-attach the actual unit (temperature now in °C, no hidden K shift)
        @test unit(cr.array[1]) == °C
        @test unit(ch_m.array[1]) == m / s
        # `read` returns bare magnitudes — the unit lives in the layer table (`layerunit`)
        @test unit(worldclim.array[1]) == NoUnits
        @test layerunit(WorldClim{Climate}, :wind) == m / s
    end
end

# `boundingbox` only reads the shipped CSV, so it runs on every platform.
@testset "Bounding boxes" begin
    scot = EcoSISTEM.ClimatePref.boundingbox("Scotland")
    @test minimum(scot.lat) == 54.63° && maximum(scot.lat) == 58.68°
    @test minimum(scot.long) == -6.23° && maximum(scot.long) == -1.76°
    # islands = true selects the island-inclusive extent
    isl = EcoSISTEM.ClimatePref.boundingbox("Scotland"; islands = true)
    @test maximum(isl.lat) == 60.86° && minimum(isl.long) == -8.65°
    # round snaps outwards to the nearest multiple, enclosing the exact box
    rnd = EcoSISTEM.ClimatePref.boundingbox("Scotland"; round = 5°)
    @test minimum(rnd.lat) == 50° && maximum(rnd.lat) == 60°
    @test minimum(rnd.long) == -10° && maximum(rnd.long) == 0°
    @test minimum(rnd.lat) ≤ minimum(scot.lat) &&
          maximum(rnd.lat) ≥ maximum(scot.lat)
    @test minimum(rnd.long) ≤ minimum(scot.long) &&
          maximum(rnd.long) ≥ maximum(scot.long)
    @test_throws ErrorException EcoSISTEM.ClimatePref.boundingbox("Atlantis")
end

# `_mask_int_fills` removes the raw integer-band fill sentinels (GDAL `typemax`/`typemin`) that a file's
# declared nodata misses — e.g. CHELSA's `0xffffffff`, which the default scaled read otherwise leaves as a
# spurious ~4.29e8. Synthetic (no download), so it runs on every platform.
@testset "integer-band fill-sentinel masking" begin
    CP = EcoSISTEM.ClimatePref
    dims = (X(1:3), Y(1:3))
    A = Float64[1 2 3; 4 5 6; 7 8 9]
    r = Rasters.Raster(A, dims)

    # unsigned band: only typemax is a fill; ordinary cells are kept
    rawU = Rasters.Raster(UInt16[typemax(UInt16) 2 3; 4 5 6;
                                 7 8 typemax(UInt16)], dims)
    maskedU = CP._mask_int_fills(r, rawU)
    @test ismissing(maskedU[1, 1]) && ismissing(maskedU[3, 3])
    @test maskedU[2, 2] == 5.0
    @test count(ismissing, Array(maskedU)) == 2

    # signed band: both typemin and typemax are fills
    rawS = Rasters.Raster(Int16[typemin(Int16) 2 3; 4 5 6; 7 8 typemax(Int16)],
                          dims)
    maskedS = CP._mask_int_fills(r, rawS)
    @test ismissing(maskedS[1, 1]) && ismissing(maskedS[3, 3])

    # a float band carries no such sentinel → returned unchanged
    @test CP._mask_int_fills(r, Rasters.Raster(Float32.(A), dims)) === r

    # the scale>1 path: aggregate(mean) must propagate the introduced missings, not error
    B = Rasters.Raster(Float64.(reshape(1:16, 4, 4)), (X(1:4), Y(1:4)))
    rawB = Rasters.Raster(reshape(UInt8.(vcat(typemax(UInt8), 2:16)), 4, 4),
                          (X(1:4), Y(1:4)))
    agg = Rasters.aggregate(mean, CP._mask_int_fills(B, rawB), 2)
    @test ismissing(agg[1, 1])                       # the fill-touching block is masked
    @test !ismissing(agg[2, 2])                      # a clean block survives
end

end
