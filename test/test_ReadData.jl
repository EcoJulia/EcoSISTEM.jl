# SPDX-License-Identifier: LGPL-3.0-or-later

module TestReadData

using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources
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
        bc = read(WorldClim{BioClim})
        cr = readCRUTS(winddir, "tavg")
        ch_b = read(CHELSA{BioClim}, 1, scale = 20)
        rf = readfile(bio1)

        @test unit(bc.array[1]) == unit(rf[1]) == unit(ch_b.array[1]) == NoUnits
    end

    @testset "Output data 2" begin
        lc = read(EarthEnv{LandCover})
        @test unit(lc.array[1]) == NoUnits
    end

    @testset "Output data 3" begin
        cr = readCRUTS(winddir, "tavg")
        wc = read(WorldClim{Climate}, :wind)
        ch_m = readCHELSA_monthly(winddir, "wind")

        @test unit(cr.array[1]) == K
        @test unit(wc.array[1]) == unit(ch_m.array[1]) == m / s
    end
end

end
