using Unitful
using Simulation.Units
using Simulation.ClimatePref

import Unitful: hr
import Dates: DateTime, DateFormat, format, Day

pathdict = (rain_max = "rain_max/ukv_daily_rain_max_", rain_mean = "rain_mean/ukv_daily_rain_mean_", rain_min = "rain_min/ukv_daily_rain_min_", temp_max = "t1o5m_max/ukv_daily_t1o5m_max_", temp_mean = "t1o5m_mean/ukv_daily_t1o5m_mean_", temp_min = "t1o5m_min/ukv_daily_t1o5m_min_", humidity_max = "sh_max/ukv_daily_sh_max_", humidity_mean = "sh_mean/ukv_daily_sh_mean_", humidity_min = "sh_min/ukv_daily_sh_min_", sunshine_max = "sw_max/ukv_daily_sw_max_", sunshine_mean = "sw_mean/ukv_daily_sw_mean_", sunshine_min = "sw_min/ukv_daily_sw_min_")
datadict = (uk_hourly = "metoffice_ukv_hourly/", uk_daily = "metoffice_ukv_daily/")
paramdict = (rain_max = "stratiform_rainfall_flux", rain_mean = "stratiform_rainfall_flux", rain_min = "stratiform_rainfall_flux", temp_max = "air_temperature", temp_mean = "air_temperature", temp_min = "air_temperature", humidity_max = "specific_humidity", humidity_mean = "specific_humidity", humidity_min = "specific_humidity", sunshine_max = "m01s01i202", sunshine_mean = "m01s01i202", sunshine_min = "m01s01i202")

"""
    MetOfficeDownload(dataset::Symbol, param::Symbol, outputfolder::String, startDate::DateTime, endDate::DateTime, read::Bool=true)

Function to download climate data from Met Office covid response, for a particular dataset (hourly or daily uk) and parameter. An output folder should be specified as well as date range. Optionally can choose to read in the data with `read`. Files will only be downloaded if they don't already exist in the `outputfolder` path.

Parameter names can be accessed with the `getMetparams()` function and dataset names similarly with `getMetdata()`.
"""
function MetOfficeDownload(dataset::Symbol, param::Symbol, outputfolder::String, startDate::DateTime, endDate::DateTime; read::Bool=true, process::Bool=true)
    # Check outputfolder exists and make folder if it doesn't
    isdir(outputfolder) || mkdir(outputfolder)

    # Format url
    coreurl = "https://metdatasa.blob.core.windows.net/covid19-response/"
    dataurl = datadict[dataset]
    endurl = pathdict[param]

    # Format dates
    df = DateFormat("yyyymmdd")
    dates = collect(startDate:Day(1):endDate)
    formatted_dates = format.(dates, df)
    for dt in formatted_dates
        file = joinpath(outputfolder, splitpath(endurl * dt *".nc")[2])
        if !isfile(file)
            download(joinpath(coreurl, dataurl, endurl * dt *".nc"), file)
        end
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
    if process
        return processMet(outputfolder, splitpath(endurl)[2], paramdict[param], times)
    else
        return readMet_raw(outputfolder, splitpath(endurl)[2], paramdict[param], times)
    end
end

getMetparams() = keys(pathdict)
getMetdata() = keys(datadict)
