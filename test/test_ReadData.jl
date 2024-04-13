module TestReadData

using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources
using Test

if !Sys.iswindows()
    if !isdir("assets")
        mkdir("assets")
        ENV["RASTERDATASOURCES_PATH"] = "assets"
        # Download layers of bioclim data and test on all read functions
        # (essentially all the same file type)
        getraster(WorldClim{BioClim}, :bio1)
        getraster(WorldClim{Climate}, :wind; month = 1:12)
        getraster(EarthEnv{LandCover})
    end

    @testset "Reading functions" begin
        @test_nowarn readbioclim("assets/WorldClim/BioClim/")
        @test_nowarn readworldclim("assets/WorldClim/Climate/wind/")
        @test_nowarn readCRUTS("assets/WorldClim/BioClim/", "tavg")
        @test_nowarn readCHELSA_monthly("assets/WorldClim/Climate/wind/",
                                        "wind")
        @test_nowarn readCHELSA_bioclim("assets/WorldClim/BioClim/")
        @test_nowarn readlc("assets/EarthEnv/LandCover/without_DISCover/")
        @test_nowarn readfile("assets/WorldClim/BioClim/wc2.1_10m_bio_1.tif")
    end

    @testset "Output data" begin
        bc = readbioclim("assets/WorldClim/BioClim/")
        cr = readCRUTS("assets/WorldClim/BioClim/", "tavg")
        ch_b = readCHELSA_bioclim("assets/WorldClim/BioClim/")
        rf = readfile("assets/WorldClim/BioClim/wc2.1_10m_bio_1.tif")

        @test unit(bc.array[1]) == unit(rf[1]) == unit(ch_b.array[1]) == NoUnits
    end

    @testset "Output data 2" begin
        lc = readlc("assets/EarthEnv/LandCover/without_DISCover/")
        @test unit(lc[1]) == NoUnits
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
