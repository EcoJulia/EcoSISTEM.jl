# SPDX-License-Identifier: LGPL-3.0-or-later

module TestReadData

using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources
using Test

if !Sys.iswindows()
    ENV["RASTERDATASOURCES_PATH"] = mkpath("assets")
    # Download layers of bioclim data and test on all read functions
    # (essentially all the same file type)
    getraster(WorldClim{BioClim}, :bio1)
    getraster(WorldClim{Climate}, :wind; month = 1:12)
    getraster(EarthEnv{LandCover})

    @testset "Reading functions" begin
        @test_nowarn read(WorldClim{Climate}, :wind)
        @test_nowarn readCRUTS("assets/WorldClim/BioClim/", "tavg")
        @test_nowarn readCHELSA_monthly("assets/WorldClim/Climate/wind/",
                                        "wind")
        @test_nowarn read(CHELSA{BioClim})
        @test_nowarn read(EarthEnv{LandCover})
        @test_nowarn readfile("assets/WorldClim/BioClim/wc2.1_10m_bio_1.tif")
    end

    @testset "Output data" begin
        bc = read(WorldClim{BioClim})
        cr = readCRUTS("assets/WorldClim/BioClim/", "tavg")
        ch_b = read(CHELSA{BioClim})
        rf = readfile("assets/WorldClim/BioClim/wc2.1_10m_bio_1.tif")

        @test unit(bc.array[1]) == unit(rf[1]) == unit(ch_b.array[1]) == NoUnits
    end

    @testset "Output data 2" begin
        lc = read(EarthEnv{LandCover})
        @test unit(lc.array[1]) == NoUnits
    end

    @testset "Output data 3" begin
        cr = readCRUTS("assets/WorldClim/BioClim/", "tavg")
        wc = readworldclim("assets/WorldClim/Climate/wind/")
        ch_m = readCHELSA_monthly("assets/WorldClim/Climate/wind/", "wind")

        @test unit(cr.array[1]) == K
        @test unit(wc.array[1]) == unit(ch_m.array[1]) == m / s
    end
end

end
