# SPDX-License-Identifier: LGPL-3.0-or-later

module TestAbioticEnv

using EcoSISTEM
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using AxisArrays
using RasterDataSources

if !Sys.iswindows()
    ENV["RASTERDATASOURCES_PATH"] = mkpath("assets")
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

    @testset "simple habitat" begin
        # TEST simplehabitatAE
        fillval = 0.0
        abenv = simplehabitatAE(fillval, grid, totalK, area)
        @test_nowarn simplehabitatAE(fillval, grid, totalK, area)
        @test_nowarn simplehabitatAE(fillval, grid, totalK, area, active)
        @test all(abenv.habitat.matrix .== fillval)
        @test size(abenv.habitat.matrix) == grid
        @test sum(abenv.budget.matrix) == totalK * area
        @test abenv.active == active
        @test all(abenv.active)

        @test EcoSISTEM._getsubcommunitynames(abenv) == abenv.names
        @test EcoSISTEM.getavailableenergy(abenv) == sum(abenv.budget.matrix)
    end

    @testset "temperature gradient" begin
        # TEST tempgradAE
        abenv = tempgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K / month)
        @test_nowarn tempgradAE(-10.0K, 10.0K, grid, totalK, area,
                                0.01K / month)
        @test minimum(abenv.habitat.matrix) == -10.0K
        @test maximum(abenv.habitat.matrix) == 10.0K
        @test size(abenv.habitat.matrix) == grid
        @test sum(abenv.budget.matrix) == totalK * area
        @test abenv.active == active
        @test all(abenv.active)
    end

    @testset "peaked temperature gradient" begin
        # TEST peakedgradAE
        abenv = peakedgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K / month)
        @test_nowarn peakedgradAE(-10.0K, 10.0K, grid, totalK, area,
                                  0.01K / month)
        @test minimum(abenv.habitat.matrix) == -10.0K
        @test maximum(abenv.habitat.matrix) == 10.0K
        @test size(abenv.habitat.matrix) == grid
        @test sum(abenv.budget.matrix) == totalK * area
        @test abenv.active == active
        @test all(abenv.active)
    end

    @testset "rainfall gradient" begin
        # TEST raingradAE
        abenv = raingradAE(0.0mm, 100.0mm, grid, totalK, area, 0.01mm / month)
        @test_nowarn raingradAE(0.0mm, 100.0mm, grid, totalK, area,
                                0.01mm / month)
        @test minimum(abenv.habitat.matrix) == 0.0mm
        @test maximum(abenv.habitat.matrix) == 100.0mm
        @test size(abenv.habitat.matrix) == grid
        @test sum(abenv.budget.matrix) == totalK * area
        @test abenv.active == active
        @test all(abenv.active)
    end

    @testset "niche habitat" begin
        # TEST simplenicheAE
        abenv = simplenicheAE(numNiches, grid, totalK, area)
        @test_nowarn simplenicheAE(numNiches, grid, totalK, area)
        @test_nowarn simplenicheAE(numNiches, grid, totalK, area, active)
        @test maximum(abenv.habitat.matrix) <= numNiches
        @test size(abenv.habitat.matrix) == grid
        @test sum(abenv.budget.matrix) == totalK * area
        @test abenv.active == active
        @test all(abenv.active)
    end

    @testset "ERA data" begin
        # TEST eraAE
        temp = AxisArray(fill(1.0K, 10, 10, 3), Axis{:latitude}(1:10),
                         Axis{:longitude}(1:10), Axis{:time}(collect(1:3) .* s))
        eratemp = ERA(temp)
        active = fill(true, 10, 10)
        totalK = 1000.0kJ / m^2
        area = 100.0km^2
        ea = eraAE(eratemp, totalK, area)
        @test_nowarn eraAE(eratemp, totalK, area)
        @test_nowarn eraAE(eratemp, totalK, area, active)
        @test size(ea.habitat.matrix) == size(temp)
        @test EcoSISTEM.getavailableenergy(ea) == totalK * area
        @test ea.active == active
        @test all(ea.active)

        solar = SolarTimeBudget(fill(10.0kJ, 10, 10, 3), 1)
        ea = eraAE(eratemp, solar, active)
    end

    @testset "Worldclim data" begin
        # TEST worldclimAE
        temp = AxisArray(fill(1.0K, 10, 10, 12),
                         Axis{:latitude}(collect(1:10) .* m),
                         Axis{:longitude}(collect(1:10) .* m),
                         Axis{:time}(collect(1:12) .* month))
        wctemp = Worldclim_monthly(temp)
        active = fill(true, 10, 10)
        totalK = 1000.0kJ / m^2
        area = 100.0km^2
        wc = worldclimAE(wctemp, totalK, area)
        @test_nowarn worldclimAE(wctemp, totalK, area)
        @test_nowarn worldclimAE(wctemp, totalK, area, active)
        @test size(wc.habitat.matrix) == size(temp)
        @test EcoSISTEM.getavailableenergy(wc) == totalK * area
        @test wc.active == active
        @test all(wc.active)

        solar = SolarTimeBudget(fill(10.0kJ, 10, 10, 3), 1)
        wc = worldclimAE(wctemp, solar, active)
    end

    @testset "Bioclim data" begin
        bio_africa = read(WorldClim{BioClim}, 1)
        active = fill(true, size(bio_africa.array))
        totalK = 1000.0kJ / km^2
        area = 100.0km^2
        bc = bioclimAE(bio_africa, totalK, area)
        @test_nowarn bioclimAE(bio_africa, totalK, area)
        @test_nowarn bioclimAE(bio_africa, totalK, area, active)
        @test size(bc.habitat.matrix) == size(bio_africa.array)
        @test isapprox(EcoSISTEM.getavailableenergy(bc), totalK * area)
        solar = SolarBudget(fill(10.0kJ, size(bio_africa.array)))
        bc = bioclimAE(bio_africa, solar, active)
    end

    @testset "LandCover data" begin
        world = read(EarthEnv{LandCover})
        world_lc = compressLC(world)
        active = fill(true, size(world_lc.array.data))
        totalK = 1000.0kJ / km^2
        area = 100.0km^2
        lc = lcAE(world_lc, totalK, area)
        @test_nowarn lcAE(world_lc, totalK, area)
        @test_nowarn lcAE(world_lc, totalK, area, active)
        @test size(lc.habitat.matrix) == size(world_lc.array)
        @test isapprox(EcoSISTEM.getavailableenergy(lc), totalK * area)
        solar = SolarBudget(fill(10.0kJ, size(world_lc.array)))
        lc = lcAE(world_lc, solar, active)
    end
end

end
