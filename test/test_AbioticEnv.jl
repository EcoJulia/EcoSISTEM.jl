# SPDX-License-Identifier: LGPL-3.0-or-later

module TestAbioticEnv

using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using AxisArrays
using RasterDataSources

if !Sys.iswindows()
    # Download layers of bioclim data and test on all read functions
    # (essentially all the same file type)
    getraster(WorldClim{BioClim}, :bio1)
    getraster(WorldClim{Climate}, :wind; month = 1:12)
    getraster(EarthEnv{LandCover})

    grid = (5, 5)
    area = 25.0km^2
    totalK = 10000.0kJ / km^2
    numNiches = 4
    active = fill(true, grid)

    @testset "simple regime" begin
        # TEST simplehabitatAE
        fillval = 0.0
        habitat = simplehabitatAE(fillval, grid, totalK, area)
        @test_nowarn simplehabitatAE(fillval, grid, totalK, area)
        @test_nowarn simplehabitatAE(fillval, grid, totalK, area, active)
        @test all(habitat.regime.matrix .== fillval)
        @test size(habitat.regime.matrix) == grid
        @test sum(habitat.supply.matrix) == totalK * area
        @test habitat.active == active
        @test all(habitat.active)

        @test EcoSISTEM._getsubcommunitynames(habitat) == habitat.names
        @test EcoSISTEM.getavailablesupply(habitat) ==
              sum(habitat.supply.matrix)
    end

    @testset "temperature gradient" begin
        # TEST tempgradAE
        habitat = tempgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K / month)
        @test_nowarn tempgradAE(-10.0K, 10.0K, grid, totalK, area,
                                0.01K / month)
        @test minimum(habitat.regime.matrix) == -10.0K
        @test maximum(habitat.regime.matrix) == 10.0K
        @test size(habitat.regime.matrix) == grid
        @test sum(habitat.supply.matrix) == totalK * area
        @test habitat.active == active
        @test all(habitat.active)
    end

    @testset "peaked temperature gradient" begin
        # TEST peakedgradAE
        habitat = peakedgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K / month)
        @test_nowarn peakedgradAE(-10.0K, 10.0K, grid, totalK, area,
                                  0.01K / month)
        @test minimum(habitat.regime.matrix) == -10.0K
        @test maximum(habitat.regime.matrix) == 10.0K
        @test size(habitat.regime.matrix) == grid
        @test sum(habitat.supply.matrix) == totalK * area
        @test habitat.active == active
        @test all(habitat.active)
    end

    @testset "rainfall gradient" begin
        # TEST raingradAE
        habitat = raingradAE(0.0mm, 100.0mm, grid, totalK, area, 0.01mm / month)
        @test_nowarn raingradAE(0.0mm, 100.0mm, grid, totalK, area,
                                0.01mm / month)
        @test minimum(habitat.regime.matrix) == 0.0mm
        @test maximum(habitat.regime.matrix) == 100.0mm
        @test size(habitat.regime.matrix) == grid
        @test sum(habitat.supply.matrix) == totalK * area
        @test habitat.active == active
        @test all(habitat.active)
    end

    @testset "niche regime" begin
        # TEST simplenicheAE
        habitat = simplenicheAE(numNiches, grid, totalK, area)
        @test_nowarn simplenicheAE(numNiches, grid, totalK, area)
        @test_nowarn simplenicheAE(numNiches, grid, totalK, area, active)
        @test maximum(habitat.regime.matrix) <= numNiches
        @test size(habitat.regime.matrix) == grid
        @test sum(habitat.supply.matrix) == totalK * area
        @test habitat.active == active
        @test all(habitat.active)
    end

    @testset "ERA data" begin
        # TEST eraAE
        temp = AxisArray(fill(1.0K, 10, 10, 3),
                         Axis{:latitude}((1:10) .* °),
                         Axis{:longitude}((1:10) .* °),
                         Axis{:time}(collect(1:3) .* s))
        eratemp = ERA(temp)
        active = fill(true, 10, 10)
        totalK = 1000.0kJ / m^2
        area = 100.0km^2
        ea = eraAE(eratemp, totalK, area)
        @test eraAE(eratemp, totalK, area) isa GridHabitat
        @test eraAE(eratemp, totalK, area, active) isa GridHabitat
        @test size(ea.regime.matrix) == size(temp)
        @test EcoSISTEM.getavailablesupply(ea) == totalK * area
        @test ea.active == active
        @test all(ea.active)

        solar = SolarTimeSupply(fill(10.0kJ, 10, 10, 3), 1)
        ea = eraAE(eratemp, solar, active)
    end

    @testset "Worldclim data" begin
        # TEST worldclimAE
        temp = AxisArray(fill(1.0K, 10, 10, 12),
                         Axis{:latitude}(collect(1:10) .* km),
                         Axis{:longitude}(collect(1:10) .* km),
                         Axis{:time}(collect(1:12) .* month))
        wctemp = ClimateRaster(WorldClim{Climate}, temp)
        active = fill(true, 10, 10)
        totalK = 1000.0kJ / m^2
        area = 100.0km^2
        wc = worldclimAE(wctemp, totalK, area)
        @test worldclimAE(wctemp, totalK, area) isa GridHabitat
        @test worldclimAE(wctemp, totalK, area, active) isa GridHabitat
        @test size(wc.regime.matrix) == size(temp)
        @test EcoSISTEM.getavailablesupply(wc) == totalK * area
        @test wc.active == active
        @test all(wc.active)

        solar = SolarTimeSupply(fill(10.0kJ, 10, 10, 3), 1)
        wc = worldclimAE(wctemp, solar, active)
    end

    @testset "Bioclim data" begin
        bio_africa = read(WorldClim{BioClim}, 1)
        active = fill(true, size(bio_africa.array))
        totalK = 1000.0kJ / km^2
        area = 100.0km^2
        bc = bioclimAE(bio_africa, totalK, area)
        @test bioclimAE(bio_africa, totalK, area) isa GridHabitat
        @test bioclimAE(bio_africa, totalK, area, active) isa GridHabitat
        @test size(bc.regime.matrix) == size(bio_africa.array)
        @test isapprox(EcoSISTEM.getavailablesupply(bc), totalK * area)
        solar = SolarSupply(fill(10.0kJ, size(bio_africa.array)))
        bc = bioclimAE(bio_africa, solar, active)
    end

    @testset "LandCover data" begin
        world = read(EarthEnv{LandCover})
        world_lc = compressLC(world)
        active = fill(true, size(world_lc.array.data))
        totalK = 1000.0kJ / km^2
        area = 100.0km^2
        lc = lcAE(world_lc, totalK, area)
        @test lcAE(world_lc, totalK, area) isa GridHabitat
        @test lcAE(world_lc, totalK, area, active) isa GridHabitat
        @test size(lc.regime.matrix) == size(world_lc.array)
        @test isapprox(EcoSISTEM.getavailablesupply(lc), totalK * area)
        solar = SolarSupply(fill(10.0kJ, size(world_lc.array)))
        lc = lcAE(world_lc, solar, active)
    end
end

end
