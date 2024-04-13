module TestEnergy

using EcoSISTEM
using EcoSISTEM.ClimatePref
using AxisArrays
using Test
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units

@testset "Requirement and budget types" begin
    numspecies = 10
    abun = fill(10, numspecies)
    # Test SimpleRequirement
    energy_vec = SimpleRequirement(fill(2.0, numspecies))
    @test_nowarn energy_vec = SimpleRequirement(fill(2.0, numspecies))
    @test_nowarn energy_vec = SimpleRequirement(collect(1.0:10.0))
    @test EcoSISTEM._getenergyusage(abun, energy_vec) == sum(abun .* energy_vec.energy)
    @test eltype(energy_vec) == typeof(energy_vec.energy[1])
    @test length(energy_vec) == length(energy_vec.energy)
    @test EcoSISTEM.numrequirements(typeof(energy_vec)) == 1

    # Test SizeRequirement
    energy_vec = SizeRequirement(fill(0.2, 10), -0.1, 1000.0km^2)
    @test_nowarn energy_vec = SizeRequirement(fill(0.2, 10), -0.1, 1000.0km^2)
    @test EcoSISTEM._getenergyusage(abun, energy_vec) == sum(abun .* energy_vec.energy)
    @test EcoSISTEM.numrequirements(typeof(energy_vec)) == 1
    @test eltype(energy_vec) == typeof(energy_vec.energy[1])
    @test length(energy_vec) == length(energy_vec.energy)

    numSpecies = 10
    # Test SolarRequirement
    energy_vec = SolarRequirement(fill(0.2*kJ, numSpecies))
    @test_nowarn energy_vec = SolarRequirement(fill(0.2*kJ, numSpecies))
    @test EcoSISTEM._getenergyusage(abun, energy_vec) == sum(abun .* energy_vec.energy)
    @test EcoSISTEM.numrequirements(typeof(energy_vec)) == 1
    @test eltype(energy_vec) == typeof(energy_vec.energy[1])
    @test length(energy_vec) == length(energy_vec.energy)

    # Test WaterRequirement
    energy_vec = WaterRequirement(fill(0.2*mm, numSpecies))
    @test_nowarn energy_vec =  WaterRequirement(fill(0.2*mm, numSpecies))
    @test EcoSISTEM._getenergyusage(abun, energy_vec) == sum(abun .* energy_vec.energy)
    @test EcoSISTEM.numrequirements(typeof(energy_vec)) == 1
    @test eltype(energy_vec) == typeof(energy_vec.energy[1])
    @test length(energy_vec) == length(energy_vec.energy)

    # Test VolWaterRequirement
    energy_vec = VolWaterRequirement(fill(20.0m^3, numSpecies))
    @test_nowarn energy_vec =  VolWaterRequirement(fill(20.0m^3, numSpecies))
    @test EcoSISTEM._getenergyusage(abun, energy_vec) == sum(abun .* energy_vec.energy)
    @test EcoSISTEM.numrequirements(typeof(energy_vec)) == 1
    @test eltype(energy_vec) == typeof(energy_vec.energy[1])
    @test length(energy_vec) == length(energy_vec.energy)

    energy_vec1 = SolarRequirement(fill(0.2*kJ, numSpecies))
    energy_vec2 = WaterRequirement(fill(0.2*mm, numSpecies))
    req = ReqCollection2(energy_vec1, energy_vec2)
    @test EcoSISTEM.numrequirements(typeof(req)) == 2
    @test eltype(req) == [typeof(energy_vec1.energy[1]), typeof(energy_vec2.energy[1])]
    @test length(req) == length(energy_vec1.energy)

    energy_vec1 = SolarRequirement(fill(0.2*kJ, numSpecies))
    energy_vec2 = WaterRequirement(fill(0.2*mm, numSpecies))
    energy_vec3 = VolWaterRequirement(fill(0.2*m^3, numSpecies))
    req = ReqCollection3(energy_vec1, energy_vec2, energy_vec3)
    @test EcoSISTEM.numrequirements(typeof(req)) == 3
    @test eltype(req) == [typeof(energy_vec1.energy[1]), typeof(energy_vec2.energy[1]), typeof(energy_vec3.energy[1])]
    @test length(req) == length(energy_vec1.energy)

    # Test SimpleBudget
    bud = Matrix{Float64}(undef, 2, 2)
    fill!(bud, 100.0)
    @test_nowarn SimpleBudget(bud)
    bud = SimpleBudget(bud)
    @test EcoSISTEM._countsubcommunities(bud) == 4
    @test EcoSISTEM._getbudget(bud) ==  bud.matrix
    @test eltype(bud) == typeof(bud.matrix[1])
    @test EcoSISTEM._getavailableenergy(bud) == sum(bud.matrix)

    # Test SolarBudget
    sol = fill(200.0*kJ, 100, 100)
    @test_nowarn SolarBudget(sol)
    bud1 = SolarBudget(sol)
    @test EcoSISTEM._countsubcommunities(bud1) == 100 * 100
    @test EcoSISTEM._getbudget(bud1) ==  bud1.matrix[:, :, 1]
    @test eltype(bud1) == typeof(bud1.matrix[1])
    @test EcoSISTEM._getavailableenergy(bud1) == sum(bud1.matrix)

    # Test WaterBudget
    water = fill(2000.0*mm, 100, 100)
    @test_nowarn WaterBudget(water)
    bud2 = WaterBudget(water)
    @test EcoSISTEM._countsubcommunities(bud2) == 100 * 100
    @test EcoSISTEM._getbudget(bud2) ==  bud2.matrix[:, :, 1]
    @test eltype(bud2) == typeof(bud2.matrix[1])
    @test EcoSISTEM._getavailableenergy(bud2) == sum(bud2.matrix)

    # Test VolWaterBudget
    water = fill(200.0*m^3, 100, 100)
    @test_nowarn VolWaterBudget(water)
    bud2 = VolWaterBudget(water)
    @test EcoSISTEM._countsubcommunities(bud2) == 100 * 100
    @test EcoSISTEM._getbudget(bud2) ==  bud2.matrix[:, :, 1]
    @test eltype(bud2) == typeof(bud2.matrix[1])
    @test EcoSISTEM._getavailableenergy(bud2) == sum(bud2.matrix)

    # Test SolarTimeBudget
    sol = fill(200.0*kJ, 100, 100, 10)
    @test_nowarn SolarTimeBudget(sol, 1)
    bud1 = SolarTimeBudget(sol, 1)
    @test EcoSISTEM._countsubcommunities(bud1) == 100 * 100
    @test EcoSISTEM._getbudget(bud1) ==  bud1.matrix[:, :, 1]
    @test eltype(bud1) == typeof(bud1.matrix[1])
    @test EcoSISTEM._getavailableenergy(bud1) == sum(bud1.matrix)

    # Test WaterTimeBudget
    water = fill(2000.0*mm, 100, 100, 10)
    @test_nowarn WaterTimeBudget(water, 1)
    bud2 = WaterTimeBudget(water, 1)
    @test EcoSISTEM._countsubcommunities(bud2) == 100 * 100
    @test EcoSISTEM._getbudget(bud2) ==  bud2.matrix[:, :, 1]
    @test eltype(bud2) == typeof(bud2.matrix[1])
    @test EcoSISTEM._getavailableenergy(bud2) == sum(bud2.matrix)

    # Test VolWaterTimeBudget
    water = fill(2000.0*m^3, 100, 100, 10)
    @test_nowarn VolWaterTimeBudget(water, 1)
    bud3 = VolWaterTimeBudget(water, 1)
    @test EcoSISTEM._countsubcommunities(bud3) == 100 * 100
    @test EcoSISTEM._getbudget(bud3) ==  bud3.matrix[:, :, 1]
    @test eltype(bud3) == typeof(bud3.matrix[1])
    @test EcoSISTEM._getavailableenergy(bud3) == sum(bud3.matrix)

    # Test BudgetCollection
    bud = BudgetCollection2(bud1, bud2)
    @test_nowarn BudgetCollection2(bud1, bud2)
    @test EcoSISTEM._countsubcommunities(bud) == 100 * 100
    @test EcoSISTEM._getbudget(bud, :b1) ==  bud1.matrix[:, :, 1]
    @test eltype(bud) == [typeof(bud1.matrix[1]), typeof(bud2.matrix[1])]
    @test EcoSISTEM._getavailableenergy(bud) == [sum(bud1.matrix), sum(bud2.matrix)]

    bud = BudgetCollection3(bud1, bud2, bud3)
    @test_nowarn BudgetCollection3(bud1, bud2, bud3)
    @test EcoSISTEM._countsubcommunities(bud) == 100 * 100
    @test EcoSISTEM._getbudget(bud, :b1) ==  bud1.matrix[:, :, 1]
    @test eltype(bud) == [typeof(bud1.matrix[1]), typeof(bud2.matrix[1]), typeof(bud3.matrix[1])]
    @test EcoSISTEM._getavailableenergy(bud) == [sum(bud1.matrix), sum(bud2.matrix), sum(bud3.matrix)]

end

@testset "Worldclim/Bioclim budgets" begin
    water = AxisArray(fill(1.0mm, 10, 10, 12), Axis{:latitude}(collect(1:10) .* m), Axis{:longitude}(collect(1:10) .* m), Axis{:time}(collect(1:12) .* month))
    wc = Worldclim_monthly(water)
    bud = WaterTimeBudget(wc, 1)
    @test_nowarn WaterTimeBudget(wc, 1)
    @test EcoSISTEM._countsubcommunities(bud) == 100
    @test EcoSISTEM._getbudget(bud) ==  bud.matrix[:, :, 1]
    @test eltype(bud) == typeof(bud.matrix[1])
    @test EcoSISTEM._getavailableenergy(bud) == sum(bud.matrix)

    solar = AxisArray(fill(1.0kJ, 10, 10, 12), Axis{:latitude}(collect(1:10) .* m), Axis{:longitude}(collect(1:10) .* m), Axis{:time}(collect(1:12) .* month))
    wc = Worldclim_monthly(solar)
    bud = SolarTimeBudget(wc, 1)
    @test_nowarn SolarTimeBudget(wc, 1)
    @test EcoSISTEM._countsubcommunities(bud) == 100
    @test EcoSISTEM._getbudget(bud) ==  bud.matrix[:, :, 1]
    @test eltype(bud) == eltype(bud.matrix)
    @test EcoSISTEM._getavailableenergy(bud) == sum(bud.matrix)
    
    water = AxisArray(fill(1.0mm, 10, 10), Axis{:latitude}(collect(1:10) .* m), Axis{:longitude}(collect(1:10) .* m))
    wc = Worldclim_bioclim(water)
    bud = WaterBudget(wc)
    @test_nowarn WaterBudget(wc)
    @test EcoSISTEM._countsubcommunities(bud) == 100
    @test EcoSISTEM._getbudget(bud) ==  bud.matrix
    @test eltype(bud) == eltype(bud.matrix)
    @test EcoSISTEM._getavailableenergy(bud) == sum(bud.matrix)

end

end
