using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using SimpleSDMLayers
using Test

@testset "Reading functions" begin
    temp = SimpleSDMPredictor(WorldClim, BioClim, 1)
    @test_nowarn readbioclim("assets/WorldClim/BioClim/10/")
    @test_nowarn readfile("assets/WorldClim/BioClim/10/wc2.1_10m_bio_1.tif")

    chelsa_temp = SimpleSDMPredictor(CHELSA, BioClim, 1)
    @test_nowarn readCHELSA_bioclim("assets/CHELSA/BioClim/")
end

@testset "Output data" begin
    bc = readbioclim("assets/WorldClim/BioClim/10/")
    rf = readfile("assets/WorldClim/BioClim/10/wc2.1_10m_bio_1.tif")

    @test unit(bc.array[1]) == NoUnits
    @test unit(rf[1]) == NoUnits
end

if isdir("assets")
    rm("assets", recursive = true)
end
