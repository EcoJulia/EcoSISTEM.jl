using Simulation
using Simulation.Units
using Simulation.ClimatePref
using Dates
using Unitful
using HTTP
using AxisArrays
using Unitful.DefaultSymbols


# Show available data
getMetdata()
getMetparams()

# Test with different variables
startDate = DateTime("2020-01-01")
endDate = DateTime("2020-01-02")
outputfolder = "test/examples/rain"
ukrain = MetOfficeDownload(:uk_daily, :rain_mean, outputfolder, startDate, endDate)

outputfolder = "test/examples/humidity"
ukhumid = MetOfficeDownload(:uk_daily, :humidity_mean, outputfolder, startDate, endDate)

outputfolder = "test/examples/sunshine"
uksun = MetOfficeDownload(:uk_daily, :sunshine_mean, outputfolder, startDate, endDate)

outputfolder = "test/examples/temp"
uktemp = MetOfficeDownload(:uk_daily, :temp_mean, outputfolder, startDate, endDate)

# Set up simple gridded environment
area = 501_768.0km^2
grid = (621, 808)
active = fill(true, grid)
active[1:10, :] .= false
epienv = ukclimateAE(ukrain, grid, area, active, NoControl())
