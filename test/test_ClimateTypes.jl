# SPDX-License-Identifier: LGPL-3.0-or-later

module TestClimateTypes

using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using AxisArrays
using Test

@testset "ECMWF types" begin
    temp = AxisArray(fill(1.0K, 10, 10, 3), Axis{:latitude}(1:10),
                     Axis{:longitude}(1:10), Axis{:time}(collect(1:3) .* s))
    @test_nowarn ERA(temp)
    @test_nowarn CERA(temp)
    @test_nowarn CRUTS(temp)
end

@testset "Worldclim types" begin
    temp = AxisArray(fill(1.0K, 10, 10, 12), Axis{:latitude}(1:10),
                     Axis{:longitude}(1:10), Axis{:time}(collect(1:12) .* s))
    @test_nowarn Worldclim_monthly(temp)
    @test_nowarn CHELSA_monthly(temp)
end

@testset "Bioclim types" begin
    temp = AxisArray(fill(1.0K, 10, 10, 19), Axis{:latitude}(1:10),
                     Axis{:longitude}(1:10), Axis{:vars}(collect(1:19)))
    @test_nowarn ClimateRaster(WorldClim{BioClim}, temp)
end

@testset "Reference types" begin
    ref = AxisArray(fill(1, 10, 10), Axis{:latitude}(1:10),
                    Axis{:longitude}(1:10))
    @test_nowarn Reference(ref)
end

end
