# SPDX-License-Identifier: LGPL-3.0-or-later

module TestHabitatUpdate

using EcoSISTEM
using Distributions
using Test
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units

include("TestCases.jl")

@testset "Habitat update" begin
    eco = Test1Ecosystem()
    @test_nowarn EcoSISTEM.HabitatUpdate(EcoSISTEM.NoChange, 0.0 / month,
                                         Unitful.Dimensions{()})
    @test_nowarn EcoSISTEM.HabitatUpdate(EcoSISTEM.TempChange, 1.0K / month,
                                         typeof(dimension(1.0K)))
    @test_nowarn EcoSISTEM.HabitatUpdate(EcoSISTEM.RainfallChange,
                                         1.0mm / month,
                                         typeof(dimension(1.0mm)))
    @test_nowarn EcoSISTEM.HabitatUpdate(EcoSISTEM.TempFluct, 1.0K / month,
                                         typeof(dimension(1.0K)))
    @test_nowarn EcoSISTEM.habitatupdate!(eco, 1month)
    @test_nowarn EcoSISTEM.budgetupdate!(eco, 1month)

    eco = TestMultiEcosystem()
    @test_nowarn EcoSISTEM.habitatupdate!(eco, 1month)
    @test_nowarn EcoSISTEM.budgetupdate!(eco, 1month)

    # Test worldclim update
    temp = AxisArray(fill(1.0K, 10, 10, 12),
                     Axis{:latitude}(collect(1:10) .* m),
                     Axis{:longitude}(collect(1:10) .* m),
                     Axis{:time}(collect(1:12) .* month))
    wctemp = Worldclim_monthly(temp)
    active = fill(true, 10, 10)
    solar = SolarTimeBudget(fill(10.0kJ, 10, 10, 3), 1)
    wc = worldclimAE(wctemp, solar, active)
    eco = TestMultiEcosystem()
    eco = Ecosystem(eco.spplist, wc, eco.relationship)
    @test_nowarn EcoSISTEM.habitatupdate!(eco, 1month)
    @test_nowarn EcoSISTEM.budgetupdate!(eco, 1month)
    @test eco.abenv.habitat.time == 2
    @test eco.abenv.budget.time == 2

    # Test era update
    temp = AxisArray(fill(1.0K, 10, 10, 3), Axis{:latitude}(1:10),
                     Axis{:longitude}(1:10), Axis{:time}(collect(1:3) .* s))
    eratemp = ERA(temp)
    water = WaterTimeBudget(fill(10.0mm, 10, 10, 3), 1)
    ea = eraAE(eratemp, water, active)
    eco = TestMultiEcosystem()
    eco = Ecosystem(eco.spplist, ea, eco.relationship)
    @test_nowarn EcoSISTEM.habitatupdate!(eco, 1month)
    @test_nowarn EcoSISTEM.budgetupdate!(eco, 1month)
    @test eco.abenv.habitat.time == 2
    @test eco.abenv.budget.time == 2
    water = VolWaterTimeBudget(fill(10.0m^3, 10, 10, 3), 1)
    ea = eraAE(eratemp, water, active)
    eco = TestMultiEcosystem()
    eco = Ecosystem(eco.spplist, ea, eco.relationship)
    @test_nowarn EcoSISTEM.habitatupdate!(eco, 1month)
    @test_nowarn EcoSISTEM.budgetupdate!(eco, 1month)
    @test eco.abenv.habitat.time == 2
    @test eco.abenv.budget.time == 2
end

end
