using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using SimpleSDMLayers
using Test

@testset "Reading functions" begin
    # Download layers of bioclim data and test on all read functions
    # (essentially all the same file type)
    temp = SimpleSDMPredictor(WorldClim, BioClim, 1:12)
    @test_nowarn readbioclim("assets/WorldClim/BioClim/10/")
    @test_nowarn readworldclim("assets/WorldClim/BioClim/10/")
    @test_nowarn readCRUTS("assets/WorldClim/BioClim/10/", "tavg")
    @test_nowarn readCHELSA_monthly("assets/WorldClim/BioClim/10/", "tavg")
    @test_nowarn readCHELSA_bioclim("assets/WorldClim/BioClim/10/")
    @test_nowarn readfile("assets/WorldClim/BioClim/10/wc2.1_10m_bio_1.tif")
end

@testset "Output data" begin
    bc = readbioclim("assets/WorldClim/BioClim/10/")
    rf = readfile("assets/WorldClim/BioClim/10/wc2.1_10m_bio_1.tif")
    wc = readworldclim("assets/WorldClim/BioClim/10/")
    cr = readCRUTS("assets/WorldClim/BioClim/10/", "tavg")
    ch_m = readCHELSA_monthly("assets/WorldClim/BioClim/10/", "tavg")
    ch_b = readCHELSA_bioclim("assets/WorldClim/BioClim/10/")

    @test unit(bc.array[1]) == unit(wc.array[1]) ==  unit(rf[1]) == unit(ch_b.array[1]) == NoUnits
    @test unit(cr.array[1]) == unit(ch_m.array[1]) == K
end

if isdir("assets")
    rm("assets", recursive = true)
end
