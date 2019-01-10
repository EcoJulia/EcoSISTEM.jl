using Simulation
using Compat.Test
using Unitful
using Unitful.DefaultSymbols
using MyUnitful

abun = fill(10, 10)
# Test SimpleRequirement
energy_vec = SimpleRequirement(fill(2, 10))
@test_nowarn energy_vec = SimpleRequirement(fill(2, 10))
@test_nowarn energy_vec = SimpleRequirement(2:10)
@test_nowarn energy_vec = SimpleRequirement(2.0:10.0)
@test Simulation._getenergyusage(abun, energy_vec) == sum(abun .* energy_vec.energy)
bud = Array{Float64, 2}(undef, 2, 2)
fill!(bud, 100.0)
@test_nowarn Simulation.SimpleBudget(bud)

# Test SizeRequirement
energy_vec = SizeRequirement(fill(0.2, 10), -0.1, 1000.0km^2)
@test_nowarn energy_vec = SizeRequirement(fill(0.2, 10), -0.1, 1000.0km^2)
@test Simulation._getenergyusage(abun, energy_vec) == sum(abun .* energy_vec.energy)

numSpecies = 10
# Test SolarRequirement
energy_vec = SolarRequirement(fill(0.2*day^-1*kJ*m^-2, numSpecies))
@test_nowarn energy_vec = SolarRequirement(fill(0.2*day^-1*kJ*m^-2, numSpecies))
@test Simulation._getenergyusage(abun, energy_vec) == sum(abun .* energy_vec.energy)

# Test WaterRequirement
energy_vec = WaterRequirement(fill(0.2*mm, numSpecies))
@test_nowarn energy_vec =  WaterRequirement(fill(0.2*mm, numSpecies))
@test Simulation._getenergyusage(abun, energy_vec) == sum(abun .* energy_vec.energy)

#Test SolarTimeBudget
sol = fill(200.0*day^-1*kJ*m^-2, 100, 100, 10)
@test_nowarn SolarTimeBudget(sol, 1)
bud1 = SolarTimeBudget(sol, 1)
@test Simulation._countsubcommunities(bud1) == 100 * 100
@test Simulation._getbudget(bud1) ==  bud1.matrix[:, :, 1]

#Test WaterTimeBudget
water = fill(2000.0*mm, 100, 100, 10)
@test_nowarn WaterTimeBudget(water, 1)
bud2 = WaterTimeBudget(water, 1)
@test Simulation._countsubcommunities(bud2) == 100 * 100
@test Simulation._getbudget(bud2) ==  bud2.matrix[:, :, 1]

# Test BudgetCollection
bud = BudgetCollection2(bud1, bud2)
@test_nowarn BudgetCollection2(bud1, bud2)
@test Simulation._countsubcommunities(bud) == 100 * 100
@test Simulation._getbudget(bud, :b1) ==  bud1.matrix[:, :, 1]
