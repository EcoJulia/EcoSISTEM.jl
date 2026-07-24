# SPDX-License-Identifier: LGPL-3.0-or-later

module TestHabitatUpdate

using EcoSISTEM
using EcoSISTEM.ClimatePref
using Distributions
using Test
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using RasterDataSources

include("TestCases.jl")

@testset "Condition update" begin
    eco = Test1Ecosystem()
    @test_nowarn EcoSISTEM.LayerUpdate(EcoSISTEM.NoChange, 0.0 / month,
                                       Unitful.Dimensions{()})
    @test_nowarn EcoSISTEM.LayerUpdate(EcoSISTEM.TempChange, 1.0K / month,
                                       typeof(dimension(1.0K)))
    @test_nowarn EcoSISTEM.LayerUpdate(EcoSISTEM.RainfallChange,
                                       1.0mm / month,
                                       typeof(dimension(1.0mm)))
    @test_nowarn EcoSISTEM.LayerUpdate(EcoSISTEM.TempFluct, 1.0K / month,
                                       typeof(dimension(1.0K)))
    @test_nowarn EcoSISTEM.regimeupdate!(eco, 1month)
    @test_nowarn EcoSISTEM.supplyupdate!(eco, 1month)

    eco = TestMultiEcosystem()
    @test_nowarn EcoSISTEM.regimeupdate!(eco, 1month)
    @test_nowarn EcoSISTEM.supplyupdate!(eco, 1month)

    # Test worldclim update
    temp = AxisArray(fill(1.0K, 10, 10, 12),
                     Axis{:latitude}(collect(1:10) .* km),
                     Axis{:longitude}(collect(1:10) .* km),
                     Axis{:time}(collect(1:12) .* month))
    worldclimtemp = ClimateRaster(WorldClim{Climate}, temp)
    active = fill(true, 10, 10)
    solar = SolarTimeSupply(fill(10.0kJ, 10, 10, 3), 1)
    worldclim = worldclimhabitat(worldclimtemp, solar, active)
    eco = TestMultiEcosystem()
    eco = Ecosystem(eco.spplist, worldclim, eco.nichefit)
    @test_nowarn EcoSISTEM.regimeupdate!(eco, 1month)
    @test_nowarn EcoSISTEM.supplyupdate!(eco, 1month)
    @test eco.habitat.regime.time == 2
    @test eco.habitat.supply.time == 2

    # Test era update
    temp = AxisArray(fill(1.0K, 10, 10, 3),
                     Axis{:latitude}((1:10) .* °),
                     Axis{:longitude}((1:10) .* °),
                     Axis{:time}(collect(1:3) .* s))
    eratemp = ERA(temp)
    water = WaterTimeSupply(fill(10.0mm, 10, 10, 3), 1)
    ea = erahabitat(eratemp, water, active)
    eco = TestMultiEcosystem()
    eco = Ecosystem(eco.spplist, ea, eco.nichefit)
    @test_nowarn EcoSISTEM.regimeupdate!(eco, 1month)
    @test_nowarn EcoSISTEM.supplyupdate!(eco, 1month)
    @test eco.habitat.regime.time == 2
    @test eco.habitat.supply.time == 2
end

@testset "Condition loss" begin
    eco = Test1Ecosystem()
    # A regime carrying a HabitatLoss change whose rate destroys every active cell
    # over one timestep (rate * 1month == 1 -> loss probability 1).
    change = EcoSISTEM.LayerUpdate(EcoSISTEM.HabitatLoss, 1.0 / month,
                                   Unitful.Dimensions{()})
    losshab = EcoSISTEM.ContinuousRegime(fill(1.0K, 10, 10), 1.0km, change)
    @test EcoSISTEM.HabitatLoss(eco, losshab, 1month) === eco
    @test all(iszero, eco.habitat.supply.matrix)
    @test all(iszero, eco.abundances.matrix)
end

end
