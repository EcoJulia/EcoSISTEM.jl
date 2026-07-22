# SPDX-License-Identifier: LGPL-3.0-or-later

module TestEnergy

using EcoSISTEM
using EcoSISTEM.ClimatePref
using AxisArrays
using Test
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using RasterDataSources

@testset "Requirement and supply types" begin
    numspecies = 10
    abun = fill(10, numspecies)
    # Test SimpleRequirement
    energy_vec = SimpleRequirement(fill(2.0, numspecies))
    @test_nowarn energy_vec = SimpleRequirement(fill(2.0, numspecies))
    @test_nowarn energy_vec = SimpleRequirement(collect(1.0:10.0))
    @test EcoSISTEM._getenergyusage(abun, energy_vec) ==
          sum(abun .* energy_vec.energy)
    @test eltype(energy_vec) == typeof(energy_vec.energy[1])
    @test length(energy_vec) == length(energy_vec.energy)
    @test EcoSISTEM.numrequirements(typeof(energy_vec)) == 1

    # Test SizeRequirement
    energy_vec = SizeRequirement(fill(0.2, 10), -0.1, 1000.0km^2)
    @test_nowarn energy_vec = SizeRequirement(fill(0.2, 10), -0.1, 1000.0km^2)
    @test EcoSISTEM._getenergyusage(abun, energy_vec) ==
          sum(abun .* energy_vec.energy)
    @test EcoSISTEM.numrequirements(typeof(energy_vec)) == 1
    @test eltype(energy_vec) == typeof(energy_vec.energy[1])
    @test length(energy_vec) == length(energy_vec.energy)

    numSpecies = 10
    # Test SolarRequirement
    energy_vec = SolarRequirement(fill(0.2 * kJ, numSpecies))
    @test_nowarn energy_vec = SolarRequirement(fill(0.2 * kJ, numSpecies))
    @test EcoSISTEM._getenergyusage(abun, energy_vec) ==
          sum(abun .* energy_vec.energy)
    @test EcoSISTEM.numrequirements(typeof(energy_vec)) == 1
    @test eltype(energy_vec) == typeof(energy_vec.energy[1])
    @test length(energy_vec) == length(energy_vec.energy)

    # Test WaterRequirement
    energy_vec = WaterRequirement(fill(0.2 * mm, numSpecies))
    @test_nowarn energy_vec = WaterRequirement(fill(0.2 * mm, numSpecies))
    @test EcoSISTEM._getenergyusage(abun, energy_vec) ==
          sum(abun .* energy_vec.energy)
    @test EcoSISTEM.numrequirements(typeof(energy_vec)) == 1
    @test eltype(energy_vec) == typeof(energy_vec.energy[1])
    @test length(energy_vec) == length(energy_vec.energy)

    # Test VolWaterRequirement
    energy_vec = VolWaterRequirement(fill(20.0m^3, numSpecies))
    @test_nowarn energy_vec = VolWaterRequirement(fill(20.0m^3, numSpecies))
    @test EcoSISTEM._getenergyusage(abun, energy_vec) ==
          sum(abun .* energy_vec.energy)
    @test EcoSISTEM.numrequirements(typeof(energy_vec)) == 1
    @test eltype(energy_vec) == typeof(energy_vec.energy[1])
    @test length(energy_vec) == length(energy_vec.energy)

    energy_vec1 = SolarRequirement(fill(0.2 * kJ, numSpecies))
    energy_vec2 = WaterRequirement(fill(0.2 * mm, numSpecies))
    req = ReqCollection2(energy_vec1, energy_vec2)
    @test EcoSISTEM.numrequirements(typeof(req)) == 2
    @test eltype(req) ==
          [typeof(energy_vec1.energy[1]), typeof(energy_vec2.energy[1])]
    @test length(req) == length(energy_vec1.energy)

    # Test SimpleSupply
    sup = Matrix{Float64}(undef, 2, 2)
    fill!(sup, 100.0)
    @test_nowarn SimpleSupply(sup)
    sup = SimpleSupply(sup)
    @test EcoSISTEM._countsubcommunities(sup) == 4
    @test EcoSISTEM._getsupply(sup) == sup.matrix
    @test eltype(sup) == typeof(sup.matrix[1])
    @test EcoSISTEM._getavailableenergy(sup) == sum(sup.matrix)

    # Test SolarSupply
    sol = fill(200.0 * kJ, 100, 100)
    @test_nowarn SolarSupply(sol)
    sup1 = SolarSupply(sol)
    @test EcoSISTEM._countsubcommunities(sup1) == 100 * 100
    @test EcoSISTEM._getsupply(sup1) == sup1.matrix[:, :, 1]
    @test eltype(sup1) == typeof(sup1.matrix[1])
    @test EcoSISTEM._getavailableenergy(sup1) == sum(sup1.matrix)

    # Test WaterSupply
    water = fill(2000.0 * mm, 100, 100)
    @test_nowarn WaterSupply(water)
    sup2 = WaterSupply(water)
    @test EcoSISTEM._countsubcommunities(sup2) == 100 * 100
    @test EcoSISTEM._getsupply(sup2) == sup2.matrix[:, :, 1]
    @test eltype(sup2) == typeof(sup2.matrix[1])
    @test EcoSISTEM._getavailableenergy(sup2) == sum(sup2.matrix)

    # Test VolWaterSupply
    water = fill(200.0 * m^3, 100, 100)
    @test_nowarn VolWaterSupply(water)
    sup2 = VolWaterSupply(water)
    @test EcoSISTEM._countsubcommunities(sup2) == 100 * 100
    @test EcoSISTEM._getsupply(sup2) == sup2.matrix[:, :, 1]
    @test eltype(sup2) == typeof(sup2.matrix[1])
    @test EcoSISTEM._getavailableenergy(sup2) == sum(sup2.matrix)

    # Test SolarTimeSupply
    sol = fill(200.0 * kJ, 100, 100, 10)
    @test_nowarn SolarTimeSupply(sol, 1)
    sup1 = SolarTimeSupply(sol, 1)
    @test EcoSISTEM._countsubcommunities(sup1) == 100 * 100
    @test EcoSISTEM._getsupply(sup1) == sup1.matrix[:, :, 1]
    @test eltype(sup1) == typeof(sup1.matrix[1])
    @test EcoSISTEM._getavailableenergy(sup1) == sum(sup1.matrix)

    # Test WaterTimeSupply
    water = fill(2000.0 * mm, 100, 100, 10)
    @test_nowarn WaterTimeSupply(water, 1)
    sup2 = WaterTimeSupply(water, 1)
    @test EcoSISTEM._countsubcommunities(sup2) == 100 * 100
    @test EcoSISTEM._getsupply(sup2) == sup2.matrix[:, :, 1]
    @test eltype(sup2) == typeof(sup2.matrix[1])
    @test EcoSISTEM._getavailableenergy(sup2) == sum(sup2.matrix)

    # Test VolWaterTimeSupply
    water = fill(2000.0 * m^3, 100, 100, 10)
    @test_nowarn VolWaterTimeSupply(water, 1)
    sup3 = VolWaterTimeSupply(water, 1)
    @test EcoSISTEM._countsubcommunities(sup3) == 100 * 100
    @test EcoSISTEM._getsupply(sup3) == sup3.matrix[:, :, 1]
    @test eltype(sup3) == typeof(sup3.matrix[1])
    @test EcoSISTEM._getavailableenergy(sup3) == sum(sup3.matrix)

    # Test SupplyCollection
    sup = SupplyCollection2(sup1, sup2)
    @test_nowarn SupplyCollection2(sup1, sup2)
    @test EcoSISTEM._countsubcommunities(sup) == 100 * 100
    @test EcoSISTEM._getsupply(sup, :one) == sup1.matrix[:, :, 1]
    @test eltype(sup) == [typeof(sup1.matrix[1]), typeof(sup2.matrix[1])]
    @test EcoSISTEM._getavailableenergy(sup) ==
          [sum(sup1.matrix), sum(sup2.matrix)]
end

@testset "Worldclim/Bioclim supplies" begin
    water = AxisArray(fill(1.0mm, 10, 10, 12),
                      Axis{:latitude}(collect(1:10) .* m),
                      Axis{:longitude}(collect(1:10) .* m),
                      Axis{:time}(collect(1:12) .* month))
    wc = ClimateRaster(WorldClim{Climate}, water)
    sup = WaterTimeSupply(wc, 1)
    @test_nowarn WaterTimeSupply(wc, 1)
    @test EcoSISTEM._countsubcommunities(sup) == 100
    @test EcoSISTEM._getsupply(sup) == sup.matrix[:, :, 1]
    @test eltype(sup) == typeof(sup.matrix[1])
    @test EcoSISTEM._getavailableenergy(sup) == sum(sup.matrix)

    solar = AxisArray(fill(1.0kJ, 10, 10, 12),
                      Axis{:latitude}(collect(1:10) .* m),
                      Axis{:longitude}(collect(1:10) .* m),
                      Axis{:time}(collect(1:12) .* month))
    wc = ClimateRaster(WorldClim{Climate}, solar)
    sup = SolarTimeSupply(wc, 1)
    @test_nowarn SolarTimeSupply(wc, 1)
    @test EcoSISTEM._countsubcommunities(sup) == 100
    @test EcoSISTEM._getsupply(sup) == sup.matrix[:, :, 1]
    @test eltype(sup) == eltype(sup.matrix)
    @test EcoSISTEM._getavailableenergy(sup) == sum(sup.matrix)

    water = AxisArray(fill(1.0mm, 10, 10), Axis{:latitude}(collect(1:10) .* m),
                      Axis{:longitude}(collect(1:10) .* m))
    wc = ClimateRaster(WorldClim{BioClim}, water)
    sup = WaterSupply(wc)
    @test_nowarn WaterSupply(wc)
    @test EcoSISTEM._countsubcommunities(sup) == 100
    @test EcoSISTEM._getsupply(sup) == sup.matrix
    @test eltype(sup) == eltype(sup.matrix)
    @test EcoSISTEM._getavailableenergy(sup) == sum(sup.matrix)
end

end
