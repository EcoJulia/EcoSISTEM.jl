using SimpleSDMLayers
using AxisArrays
using Unitful
using Unitful.DefaultSymbols
import AxisArrays.Axis

function Bioclim(sdm::SimpleSDMPredictor, units::Unitful.Units)
    sdm = convert(Float64, sdm)
    sdm.grid[isnothing.(sdm.grid)] .= NaN
    unit_grid = sdm.grid .* units  
    if units == °C
        unit_grid = uconvert.(K, unit_grid)
    end
    lat = collect(latitudes(sdm)) .* °
    lon = collect(longitudes(sdm)) .* °
    ax_array = AxisArray(unit_grid, Axis{:latitude}(lat), Axis{:longitude}(lon))
    return Bioclim(ax_array)
end

function CHELSA_bioclim(sdm::SimpleSDMPredictor, units::Unitful.Units)
    sdm = convert(Float64, sdm)
    sdm.grid[isnothing.(sdm.grid)] .= NaN
    unit_grid = sdm.grid .* units  
    if units == °C
        unit_grid = uconvert.(K, unit_grid)
    end
    lat = collect(latitudes(sdm)) .* °
    lon = collect(longitudes(sdm)) .* °
    ax_array = AxisArray(unit_grid, Axis{:latitude}(lat), Axis{:longitude}(lon))
    return CHELSA_bioclim(ax_array)
end

function Landcover(sdm::SimpleSDMPredictor)
    sdm = convert(Float64, sdm)
    sdm.grid[isnothing.(sdm.grid)] .= NaN
    lat = collect(latitudes(sdm)) .* °
    lon = collect(longitudes(sdm)) .* °
    ax_array = AxisArray(sdm.grid, Axis{:latitude}(lat), Axis{:longitude}(lon))
    return Landcover(ax_array)
end