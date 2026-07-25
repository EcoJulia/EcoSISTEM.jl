# SPDX-License-Identifier: LGPL-3.0-or-later

module TestEnergy

using EcoSISTEM
using EcoSISTEM.ClimatePref
using AxisArrays
using InteractiveUtils
using Test
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using RasterDataSources

@testset "Demand and supply types" begin
    numspecies = 10
    abun = fill(10, numspecies)
    # Test SimpleDemand
    resource_vec = SimpleDemand(fill(2.0, numspecies))
    @test_nowarn resource_vec = SimpleDemand(fill(2.0, numspecies))
    @test_nowarn resource_vec = SimpleDemand(collect(1.0:10.0))
    @test EcoSISTEM._getdemand(abun, resource_vec) ==
          sum(abun .* resource_vec.resource)
    @test eltype(resource_vec) == typeof(resource_vec.resource[1])
    @test length(resource_vec) == length(resource_vec.resource)
    @test EcoSISTEM.numdemands(typeof(resource_vec)) == 1

    # Test SizeDemand
    resource_vec = SizeDemand(fill(0.2, 10), -0.1, 1000.0km^2)
    @test_nowarn resource_vec = SizeDemand(fill(0.2, 10), -0.1, 1000.0km^2)
    @test EcoSISTEM._getdemand(abun, resource_vec) ==
          sum(abun .* resource_vec.resource)
    @test EcoSISTEM.numdemands(typeof(resource_vec)) == 1
    @test eltype(resource_vec) == typeof(resource_vec.resource[1])
    @test length(resource_vec) == length(resource_vec.resource)

    numSpecies = 10
    # Test SolarDemand
    resource_vec = SolarDemand(fill(0.2 * kJ / day, numSpecies))
    @test_nowarn resource_vec = SolarDemand(fill(0.2 * kJ / day, numSpecies))
    @test EcoSISTEM._getdemand(abun, resource_vec) ==
          sum(abun .* resource_vec.resource)
    @test EcoSISTEM.numdemands(typeof(resource_vec)) == 1
    @test eltype(resource_vec) == typeof(resource_vec.resource[1])
    @test length(resource_vec) == length(resource_vec.resource)

    # Test WaterDemand
    resource_vec = WaterDemand(fill(0.2 * Unitful.L / day, numSpecies))
    @test_nowarn resource_vec = WaterDemand(fill(0.2 * Unitful.L / day,
                                                 numSpecies))
    @test EcoSISTEM._getdemand(abun, resource_vec) ==
          sum(abun .* resource_vec.resource)
    @test EcoSISTEM.numdemands(typeof(resource_vec)) == 1
    @test eltype(resource_vec) == typeof(resource_vec.resource[1])
    @test length(resource_vec) == length(resource_vec.resource)

    resource_vec1 = SolarDemand(fill(0.2 * kJ / day, numSpecies))
    resource_vec2 = WaterDemand(fill(0.2 * Unitful.L / day, numSpecies))
    demand = DemandCollection2(resource_vec1, resource_vec2)
    @test EcoSISTEM.numdemands(typeof(demand)) == 2
    @test eltype(demand) ==
          [typeof(resource_vec1.resource[1]), typeof(resource_vec2.resource[1])]
    @test length(demand) == length(resource_vec1.resource)

    # Test SimpleSupply
    supply = Matrix{Float64}(undef, 2, 2)
    fill!(supply, 100.0)
    @test_nowarn SimpleSupply(supply)
    supply = SimpleSupply(supply)
    @test EcoSISTEM._countsubcommunities(supply) == 4
    @test EcoSISTEM._getsupply(supply) == supply.matrix
    @test eltype(supply) == typeof(supply.matrix[1])
    @test EcoSISTEM._getavailablesupply(supply) == sum(supply.matrix)

    # Test SolarSupply
    sol = fill(200.0 * kJ / day, 100, 100)
    @test_nowarn SolarSupply(sol)
    supply1 = SolarSupply(sol)
    @test EcoSISTEM._countsubcommunities(supply1) == 100 * 100
    @test EcoSISTEM._getsupply(supply1) == supply1.matrix[:, :, 1]
    @test eltype(supply1) == typeof(supply1.matrix[1])
    @test EcoSISTEM._getavailablesupply(supply1) == sum(supply1.matrix)

    # Test WaterSupply
    water = fill(2000.0 * Unitful.L / day, 100, 100)
    @test_nowarn WaterSupply(water)
    supply2 = WaterSupply(water)
    @test EcoSISTEM._countsubcommunities(supply2) == 100 * 100
    @test EcoSISTEM._getsupply(supply2) == supply2.matrix[:, :, 1]
    @test eltype(supply2) == typeof(supply2.matrix[1])
    @test EcoSISTEM._getavailablesupply(supply2) == sum(supply2.matrix)

    # Test SolarTimeSupply
    sol = fill(200.0 * kJ / day, 100, 100, 10)
    @test_nowarn SolarTimeSupply(sol, 1)
    supply1 = SolarTimeSupply(sol, 1)
    @test EcoSISTEM._countsubcommunities(supply1) == 100 * 100
    @test EcoSISTEM._getsupply(supply1) == supply1.matrix[:, :, 1]
    @test eltype(supply1) == typeof(supply1.matrix[1])
    @test EcoSISTEM._getavailablesupply(supply1) == sum(supply1.matrix)

    # Test WaterTimeSupply
    water = fill(2000.0 * Unitful.L / day, 100, 100, 10)
    @test_nowarn WaterTimeSupply(water, 1)
    supply2 = WaterTimeSupply(water, 1)
    @test EcoSISTEM._countsubcommunities(supply2) == 100 * 100
    @test EcoSISTEM._getsupply(supply2) == supply2.matrix[:, :, 1]
    @test eltype(supply2) == typeof(supply2.matrix[1])
    @test EcoSISTEM._getavailablesupply(supply2) == sum(supply2.matrix)

    # Test SupplyCollection
    supply = SupplyCollection2(supply1, supply2)
    @test_nowarn SupplyCollection2(supply1, supply2)
    @test EcoSISTEM._countsubcommunities(supply) == 100 * 100
    @test EcoSISTEM._getsupply(supply, :one) == supply1.matrix[:, :, 1]
    @test eltype(supply) ==
          [typeof(supply1.matrix[1]), typeof(supply2.matrix[1])]
    @test EcoSISTEM._getavailablesupply(supply) ==
          [sum(supply1.matrix), sum(supply2.matrix)]
end

@testset "Worldclim/Bioclim supplies" begin
    # The hand constructors take a raw ClimateRaster whose element type is the areal rate
    # (e.g. `mm/day`) at its own native resolution — `cancel` (dispatched on dimension)
    # converts it to the absolute `L/day`/`kJ/day` supply value against the grid's cell area,
    # so the resulting `supply.matrix` values are scaled by that area, not equal to the raw
    # input; the assertions below compare the built supply against itself, not a hardcoded
    # number, so they hold regardless of the exact scaling.
    water = AxisArray(fill(1.0mm / day, 10, 10, 12),
                      Axis{:latitude}(collect(1:10) .* m),
                      Axis{:longitude}(collect(1:10) .* m),
                      Axis{:time}(collect(1:12) .* month))
    worldclim = ClimateRaster(WorldClim{Climate}, water)
    supply = WaterTimeSupply(worldclim, 1)
    @test_nowarn WaterTimeSupply(worldclim, 1)
    @test EcoSISTEM._countsubcommunities(supply) == 100
    @test EcoSISTEM._getsupply(supply) == supply.matrix[:, :, 1]
    @test eltype(supply) == typeof(supply.matrix[1])
    @test EcoSISTEM._getavailablesupply(supply) == sum(supply.matrix)

    solar = AxisArray(fill(1.0kJ / m^2 / day, 10, 10, 12),
                      Axis{:latitude}(collect(1:10) .* m),
                      Axis{:longitude}(collect(1:10) .* m),
                      Axis{:time}(collect(1:12) .* month))
    worldclim = ClimateRaster(WorldClim{Climate}, solar)
    supply = SolarTimeSupply(worldclim, 1)
    @test_nowarn SolarTimeSupply(worldclim, 1)
    @test EcoSISTEM._countsubcommunities(supply) == 100
    @test EcoSISTEM._getsupply(supply) == supply.matrix[:, :, 1]
    @test eltype(supply) == eltype(supply.matrix)
    @test EcoSISTEM._getavailablesupply(supply) == sum(supply.matrix)

    water = AxisArray(fill(1.0mm / day, 10, 10),
                      Axis{:latitude}(collect(1:10) .* m),
                      Axis{:longitude}(collect(1:10) .* m))
    worldclim = ClimateRaster(WorldClim{BioClim}, water)
    supply = WaterSupply(worldclim)
    @test_nowarn WaterSupply(worldclim)
    @test EcoSISTEM._countsubcommunities(supply) == 100
    @test EcoSISTEM._getsupply(supply) == supply.matrix
    @test eltype(supply) == eltype(supply.matrix)
    @test EcoSISTEM._getavailablesupply(supply) == sum(supply.matrix)
end

# All concrete leaves of an abstract type, recursing through intermediate abstract types
# (mirrors `ClimatePref.LayerUnits._leafaxes`'s pattern for `NicheAxis`).
function _leaftypes(T)
    subs = subtypes(T)
    return isempty(subs) ? [T] : reduce(vcat, _leaftypes.(subs); init = Type[])
end

# guard: every concrete single-resource Demand type is genuinely rate-dimensioned —
# `dimension(element) / _basedimension(element) == 𝐓⁻¹` exactly (not merely "any negative
# time exponent" — energy already embeds 𝐓⁻² in its own definition, see `_basedimension`'s
# docstring in `NicheInfo.jl`). Catches a future Demand type that doesn't follow the pattern.
@testset "guard: every Demand type is substance × 𝐓⁻¹" begin
    sample(::Type{SimpleDemand}) = SimpleDemand(fill(1.0, 3))
    sample(::Type{SizeDemand}) = SizeDemand(fill(1.0, 3), 1.0, 1.0km^2)
    sample(::Type{SolarDemand}) = SolarDemand(fill(1.0kJ / day, 3))
    sample(::Type{WaterDemand}) = WaterDemand(fill(1.0 * Unitful.L / day, 3))
    concrete = filter(!isabstracttype, _leaftypes(EcoSISTEM.Abstract1Demand))
    @test !isempty(concrete)
    for T in concrete
        d = sample(T)
        elt = eltype(d)
        @test dimension(elt) / EcoSISTEM._basedimension(elt) == Unitful.𝐓^-1
    end
end

end
