using Simulation.ClimatePref
using Unitful.DefaultSymbols
using Simulation.Units

dir = "test/examples/dailyrain/"
file = "ukv"
dim = collect(1day:1day:3days)
param = "stratiform_rainfall_flux"
ukrain = readMet(dir, file, param, dim)

using Plots
using Unitful
display(heatmap(ustrip.(ukrain[:, :, 1]')))

dir = "test/examples/temp/"
file = "ukv"
dim = collect(1day:1day:3days)
param = "air_temperature"
uktemp = readMet(dir, file, param, dim)

display(heatmap(ustrip.(uktemp[:, :, 1]')))
