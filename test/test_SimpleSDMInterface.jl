using EcoSISTEM.ClimatePref
using AxisArrays
using SimpleSDMLayers
using Unitful
using Unitful.DefaultSymbols
using Test

africa_temp = SimpleSDMPredictor(WorldClim, BioClim, 1, left = -25, right = 50, bottom = -35, top = 40)
bio_africa = Bioclim(africa_temp, °C)
@test unit(bio_africa.array[1]) == K

africa_temp_chelsa = SimpleSDMPredictor(CHELSA, BioClim, 1, left = -25, right = 50, bottom = -35, top = 40)
chelsa_africa = CHELSA_bioclim(africa_temp, °C)
@test unit(chelsa_africa.array[1]) == K

africa_lc = SimpleSDMPredictor(EarthEnv, LandCover, 1, left = -25, right = 50, bottom = -35, top = 40)
lc_africa = Landcover(africa_temp)
@test unit(lc_africa.array[1]) == NoUnits

if isdir("assets")
    rm("assets", recursive = true)
end