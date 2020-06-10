using Unitful
using Dates
using Simulation.Units
using Simulation.ClimatePref

import Unitful.hr

pathdict = (rain_max = "rain_max/ukv_daily_rain_max_", rain_mean = "rain_mean/ukv_daily_rain_mean_", rain_min = "rain_min/ukv_daily_rain_min_")
datadict = (uk_hourly = "metoffice_ukv_hourly/", uk_daily = "metoffice_ukv_daily/")
paramdict = (rain_max = "stratiform_rainfall_flux", rain_mean = "stratiform_rainfall_flux", rain_min = "stratiform_rainfall_flux")

"""
    MetOfficeDownload(dataset::Symbol, param::Symbol, outputfolder::String, startDate::DateTime, endDate::DateTime, read::Bool=true)

Function to download climate data from Met Office covid response, for a particular dataset (hourly or daily uk) and parameter. An output folder should be specified as well as date range. Optionally can choose to read in the data with `read`.

Parameter names can be accessed with the `getMetparams()` function and dataset names similarly with `getMetdata()`.
"""
function MetOfficeDownload(dataset::Symbol, param::Symbol, outputfolder::String, startDate::DateTime, endDate::DateTime, read::Bool=true)
    # Format url
    coreurl = "https://metdatasa.blob.core.windows.net/covid19-response/"
    dataurl = datadict[dataset]
    endurl = pathdict[param]

    # Format dates
    df = DateFormat("yyyymmdd")
    dates = collect(startDate:Dates.Day(1):endDate)
    formatted_dates = Dates.format.(dates, df)
    for dt in formatted_dates
        download(joinpath(coreurl, dataurl, endurl * dt *".nc"), joinpath(outputfolder, splitpath(endurl * dt *".nc")[2]))
    end

    # Read data from file
    !read && return "Files stored in $outputfolder"
    if dataset == :uk_daily
        times = collect(1days:1days:length(dates)*days)
    elseif dataset == :uk_hourly
        times = collect(1hr:hr:length(dates)*hr)
    else
        error("No other datasets currently implemented")
    end
    return readMet(outputfolder, splitpath(endurl)[2], paramdict[param], times)
end

getMetparams() = keys(pathdict)
getMetdata() = keys(datadict)
