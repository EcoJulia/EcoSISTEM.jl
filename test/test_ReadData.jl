using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources
using Test

if !isdir("assets")
    mkdir("assets")
end
ENV["RASTERDATASOURCES_PATH"] = "assets"
# Download layers of bioclim data and test on all read functions
# (essentially all the same file type)
temp = getraster(WorldClim{BioClim}, 1:12)
lc = getraster(EarthEnv{LandCover})

@testset "Reading functions" begin
    @test_nowarn readbioclim("assets/WorldClim/BioClim/")
    @test_nowarn readworldclim("assets/WorldClim/BioClim/")
    @test_nowarn readCRUTS("assets/WorldClim/BioClim/", "tavg")
    @test_nowarn readCHELSA_monthly("assets/WorldClim/BioClim/", "tavg")
    @test_nowarn readCHELSA_bioclim("assets/WorldClim/BioClim/")
    @test_nowarn readlc("assets/EarthEnv/LandCover/without_DISCover/")
    @test_nowarn readfile("assets/WorldClim/BioClim/wc2.1_10m_bio_1.tif")
end

@testset "Output data" begin
    bc = readbioclim("assets/WorldClim/BioClim/")
    rf = readfile("assets/WorldClim/BioClim/wc2.1_10m_bio_1.tif")
    wc = readworldclim("assets/WorldClim/BioClim/")
    cr = readCRUTS("assets/WorldClim/BioClim/", "tavg")
    ch_m = readCHELSA_monthly("assets/WorldClim/BioClim/", "tavg")
    ch_b = readCHELSA_bioclim("assets/WorldClim/BioClim/")
    lc = readlc("assets/EarthEnv/LandCover/without_DISCover/")

    @test unit(bc.array[1]) == unit(wc.array[1]) ==  unit(rf[1]) == unit(ch_b.array[1]) == unit(lc[1])== NoUnits
    @test unit(cr.array[1]) == unit(ch_m.array[1]) == K
end

if isdir("assets/WorldClim")
    rm("assets/WorldClim", recursive = true)
end