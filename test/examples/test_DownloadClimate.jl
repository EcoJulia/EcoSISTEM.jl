using RCall
using EcoSISTEM
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using Unitful
using Unitful.DefaultSymbols
import Dates.DateTime
using Test

# Download climate data and write to HDF5
startDate = DateTime("2020-01-01")
endDate = DateTime("2020-01-02")
uktemp = MetOfficeDownload(:uk_daily, :temp_mean, "TEMP", startDate, endDate)
@test typeof(uktemp[1]) <: Unitful.Temperature
@test length(size(uktemp)) == 3
@test size(uktemp, 3) == 2

outputfolder = "TEMP/climate.h5"
writeMet(uktemp, "temp", outputfolder)
uktemp2 = readMet(outputfolder, "temp")
@test all(uktemp[.!isnan.(uktemp)] .== uktemp2[.!isnan.(uktemp2)])
rm("TEMP", recursive = true)
