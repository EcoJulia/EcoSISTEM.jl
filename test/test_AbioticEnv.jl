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
        # TEST simplehabitat
        fillval = 0.0
        habitat = simplehabitat(fillval, grid, totalK, area)
        @test_nowarn simplehabitat(fillval, grid, totalK, area)
        @test_nowarn simplehabitat(fillval, grid, totalK, area, active)
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
        # TEST tempgradhabitat
        habitat = tempgradhabitat(-10.0K, 10.0K, grid, totalK, area,
                                  0.01K / month)
        @test_nowarn tempgradhabitat(-10.0K, 10.0K, grid, totalK, area,
                                     0.01K / month)
        @test minimum(habitat.regime.matrix) == -10.0K
        @test maximum(habitat.regime.matrix) == 10.0K
        @test size(habitat.regime.matrix) == grid
        @test sum(habitat.supply.matrix) == totalK * area
        @test habitat.active == active
        @test all(habitat.active)
    end

    @testset "peaked temperature gradient" begin
        # TEST peakedgradhabitat
        habitat = peakedgradhabitat(-10.0K, 10.0K, grid, totalK, area,
                                    0.01K / month)
        @test_nowarn peakedgradhabitat(-10.0K, 10.0K, grid, totalK, area,
                                       0.01K / month)
        @test minimum(habitat.regime.matrix) == -10.0K
        @test maximum(habitat.regime.matrix) == 10.0K
        @test size(habitat.regime.matrix) == grid
        @test sum(habitat.supply.matrix) == totalK * area
        @test habitat.active == active
        @test all(habitat.active)
    end

    @testset "rainfall gradient" begin
        # TEST raingradhabitat
        habitat = raingradhabitat(0.0mm, 100.0mm, grid, totalK, area,
                                  0.01mm / month)
        @test_nowarn raingradhabitat(0.0mm, 100.0mm, grid, totalK, area,
                                     0.01mm / month)
        @test minimum(habitat.regime.matrix) == 0.0mm
        @test maximum(habitat.regime.matrix) == 100.0mm
        @test size(habitat.regime.matrix) == grid
        @test sum(habitat.supply.matrix) == totalK * area
        @test habitat.active == active
        @test all(habitat.active)
    end

    @testset "niche regime" begin
        # TEST simplenichehabitat
        habitat = simplenichehabitat(numNiches, grid, totalK, area)
        @test_nowarn simplenichehabitat(numNiches, grid, totalK, area)
        @test_nowarn simplenichehabitat(numNiches, grid, totalK, area, active)
        @test maximum(habitat.regime.matrix) <= numNiches
        @test size(habitat.regime.matrix) == grid
        @test sum(habitat.supply.matrix) == totalK * area
        @test habitat.active == active
        @test all(habitat.active)
    end

    @testset "ERA data" begin
        # TEST erahabitat
        temp = AxisArray(fill(1.0K, 10, 10, 3),
                         Axis{:latitude}((1:10) .* °),
                         Axis{:longitude}((1:10) .* °),
                         Axis{:time}(collect(1:3) .* s))
        eratemp = ERA(temp)
        active = fill(true, 10, 10)
        totalK = 1000.0kJ / m^2
        area = 100.0km^2
        era = erahabitat(eratemp, totalK, area)
        @test erahabitat(eratemp, totalK, area) isa GridHabitat
        @test erahabitat(eratemp, totalK, area, active) isa GridHabitat
        @test size(era.regime.matrix) == size(temp)
        @test EcoSISTEM.getavailablesupply(era) == totalK * area
        @test era.active == active
        @test all(era.active)
        # the niche `axis` is threaded through onto the built regime (default `Unclassified`)
        @test EcoSISTEM.axisof(era.regime) === Unclassified
        @test EcoSISTEM.axisof(erahabitat(eratemp, totalK, area;
                                          axis = Temperature).regime) ===
              Temperature

        solar = SolarTimeSupply(fill(10.0kJ, 10, 10, 3), 1)
        era = erahabitat(eratemp, solar, active)
    end

    @testset "Worldclim data" begin
        # TEST worldclimhabitat
        temp = AxisArray(fill(1.0K, 10, 10, 12),
                         Axis{:latitude}(collect(1:10) .* km),
                         Axis{:longitude}(collect(1:10) .* km),
                         Axis{:time}(collect(1:12) .* month))
        worldclimtemp = ClimateRaster(WorldClim{Climate}, temp)
        active = fill(true, 10, 10)
        totalK = 1000.0kJ / m^2
        area = 100.0km^2
        worldclim = worldclimhabitat(worldclimtemp, totalK, area)
        @test worldclimhabitat(worldclimtemp, totalK, area) isa GridHabitat
        @test worldclimhabitat(worldclimtemp, totalK, area, active) isa
              GridHabitat
        @test size(worldclim.regime.matrix) == size(temp)
        @test EcoSISTEM.getavailablesupply(worldclim) == totalK * area
        @test worldclim.active == active
        @test all(worldclim.active)
        # the niche `axis` is threaded through onto the built regime (default `Unclassified`)
        @test EcoSISTEM.axisof(worldclim.regime) === Unclassified
        @test EcoSISTEM.axisof(worldclimhabitat(worldclimtemp, totalK, area;
                                                axis = Temperature).regime) ===
              Temperature

        solar = SolarTimeSupply(fill(10.0kJ, 10, 10, 3), 1)
        worldclim = worldclimhabitat(worldclimtemp, solar, active)
    end

    @testset "Bioclim data" begin
        bio_africa = read(WorldClim{BioClim}, 1)
        active = fill(true, size(bio_africa.array))
        totalK = 1000.0kJ / km^2
        area = 100.0km^2
        bioclim = bioclimhabitat(bio_africa, totalK, area)
        @test bioclimhabitat(bio_africa, totalK, area) isa GridHabitat
        @test bioclimhabitat(bio_africa, totalK, area, active) isa GridHabitat
        @test size(bioclim.regime.matrix) == size(bio_africa.array)
        @test isapprox(EcoSISTEM.getavailablesupply(bioclim), totalK * area)
        solar = SolarSupply(fill(10.0kJ, size(bio_africa.array)))
        bioclim = bioclimhabitat(bio_africa, solar, active)
    end

    @testset "LandCover data" begin
        world = read(EarthEnv{LandCover})
        world_landcover = compressLandCover(world)
        active = fill(true, size(world_landcover.array.data))
        totalK = 1000.0kJ / km^2
        area = 100.0km^2
        landcover = landcoverhabitat(world_landcover, totalK, area)
        @test landcoverhabitat(world_landcover, totalK, area) isa GridHabitat
        @test landcoverhabitat(world_landcover, totalK, area, active) isa
              GridHabitat
        @test size(landcover.regime.matrix) == size(world_landcover.array)
        @test isapprox(EcoSISTEM.getavailablesupply(landcover), totalK * area)
        solar = SolarSupply(fill(10.0kJ, size(world_landcover.array)))
        landcover = landcoverhabitat(world_landcover, solar, active)
    end
end

end
