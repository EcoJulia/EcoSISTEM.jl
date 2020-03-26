using Simulation
using Unitful.DefaultSymbols
using Compat.Test
using Simulation.Units
using Simulation.ClimatePref
using AxisArrays

grid = (5, 5)
area = 25.0km^2
totalK = 10000.0
numNiches = 4
active = fill(true, grid)

## TEST simplehabitatAE
@test_nowarn abenv = simplehabitatAE(0.0, grid, totalK, area)
@test_nowarn abenv = simplehabitatAE(0.0, grid, totalK, area, active)

## TEST tempgradAE - temp offline until Unitful is sorted
#@test_nowarn abenv = tempgradAE(-10.0K, 10.0K, grid, totalK, area,
#    0.01K/month)
## TEST simplenicheAE
@test_nowarn abenv = simplenicheAE(numNiches, grid, totalK, area)
@test_nowarn abenv = simplenicheAE(numNiches, grid, totalK, area, active)

## TEST eraAE and worldclimAE
era = ClimatePref.TestERA()
wc = ClimatePref.TestWorldclim()
wc = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3}, wc.array[-10° .. 60°, 35° .. 80°,:])
bud = SolarTimeBudget(wc, 1)
@test_nowarn abenv = eraAE(era, totalK)
@test_nowarn abenv = eraAE(era, totalK, active)
@test_nowarn abenv = eraAE(era, bud, active)
