# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using NetCDF
using IndexedTables

function checkbounds(lat::Vector{typeof(1.0°)}, long::Vector{typeof(1.0°)})
    return all(-180.0° .< lat .<= 180.0°) && all(-90.0° .<= long .<= 90.0°)
end

function calculatestep(dat::AbstractClimate, dim::Int64)
    return AxisArrays.axes(dat.array, dim).val[2] -
           AxisArrays.axes(dat.array, dim).val[1]
end

function _extract(lat, long, dat, dim::Unitful.Time, thisstep)
    res = map((i, j) -> dat[(i - thisstep / 2) .. (i + thisstep / 2),
                            (j - thisstep / 2) .. (j + thisstep / 2),
                            dim][1, 1, 1],
              lat, long)
    return transpose(hcat(res...))
end

"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        wc::Worldclim_monthly, dim::Unitful.Time)

Function to extract values from a worldclim object, at specified x, y locations and
time, `dim`.
"""
function extractvalues(lat::Vector{typeof(1.0°)}, long::Vector{typeof(1.0°)},
                       wc::Worldclim_monthly, dim::Unitful.Time)
    checkbounds(lat, long) || error("Coordinates out of bounds")
    thisstep = calculatestep(wc, 1)
    res = map((i, j) -> wc.array[(i - thisstep / 2) .. (i + thisstep / 2),
                                 (j - thisstep / 2) .. (j + thisstep / 2),
                                 dim][1, 1, 1], lat, long)
    return transpose(hcat(res...))
end

"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        wc::Worldclim_monthly, dim::Unitful.Time)

Function to extract values from a worldclim object, at specified x, y locations and
over a range of times, `dim`.
"""
function extractvalues(lat::Vector{typeof(1.0°)}, long::Vector{typeof(1.0°)},
                       wc::Worldclim_monthly, dim::StepRange{typeof(1month)})
    checkbounds(lat, long) || error("Coordinates out of bounds")
    res = map((i, j) -> wc.array[(i - thisstep / 2) .. (i + thisstep / 2),
                                 (j - thisstep / 2) .. (j + thisstep / 2),
                                 dim], lat, long)
    return transpose(hcat(res...))
end

"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        bc::ClimateRaster, dim::Unitful.Time)

Function to extract values from a bioclim object, at specified lat, long locations and time, `dim`.
"""
function extractvalues(lat::Vector{typeof(1.0°)}, long::Vector{typeof(1.0°)},
                       bc::ClimateRaster, val::Int64)
    checkbounds(lat, long) || error("Coordinates out of bounds")
    thisstep = axes(bc.array, 1).val[2] - axes(bc.array, 1).val[1]
    res = map((i, j) -> bc.array[(i - thisstep / 2) .. (i + thisstep / 2),
                                 (j - thisstep / 2) .. (j + thisstep / 2),
                                 val], lat, long)
    return transpose(hcat(res...))
end

"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        bc::ClimateRaster, dim::Unitful.Time)

Function to extract values from a bioclim object, at specified x, y locations and over a time period, `dim`.
"""
function extractvalues(lat::Vector{typeof(1.0°)}, long::Vector{typeof(1.0°)},
                       bc::ClimateRaster, vals::StepRange{Int64})
    checkbounds(lat, long) || error("Coordinates out of bounds")
    thisstep = axes(bc.array, 1).val[2] - axes(bc.array, 1).val[1]
    res = map((i, j) -> bc.array[(i - thisstep / 2) .. (i + thisstep / 2),
                                 (j - thisstep / 2) .. (j + thisstep / 2),
                                 start(vals) .. last(vals)][1, 1, :], x, y)
    return transpose(hcat(res...))
end

function extractvalues(lat::Union{Missing, typeof(1.0°)},
                       long::Union{Missing, typeof(1.0°)},
                       cal_yr::Union{Missing, Int64}, era::ERA)
    startyr = ustrip(uconvert(year, axes(era.array)[3].val[1]))
    endyr = ustrip(uconvert(year, axes(era.array)[3].val[end]))
    if any(ismissing.([lat, long, cal_yr])) || cal_yr < startyr ||
       cal_yr > endyr
        return fill(NaN, 12) .* unit(era.array[1, 1, 1])
    else
        -180.0° ≤ lat ≤ 180.0° || error("Latitude coordinate is out of bounds")
        -90.0° < long < 90.0° || error("Longitude coordinate is out of bounds")
        thisstep1 = AxisArrays.axes(era.array, 1).val[2] -
                    AxisArrays.axes(era.array, 1).val[1]
        thisstep2 = AxisArrays.axes(era.array, 2).val[2] -
                    AxisArrays.axes(era.array, 2).val[1]
        time = cal_yr * year
        long += 0.375°
        return era.array[(long - thisstep1 / 2) .. (long + thisstep1 / 2),
                         (lat - thisstep2 / 2) .. (lat + thisstep2 / 2),
                         time .. (time + 11month)][1, 1, :]
    end
end

function extractvalues(lat::Vector{typeof(1.0°)}, long::Vector{typeof(1.0°)},
                       years::Vector{Int64}, era::ERA)
    checkbounds(lat, long) || error("Coordinates out of bounds")
    startyr = ustrip(uconvert(year, axes(era.array)[3].val[1]))
    endyr = ustrip(uconvert(year, axes(era.array)[3].val[end]))
    thisstep1 = axes(era.array, 1).val[2] - axes(era.array, 1).val[1]
    thisstep2 = axes(era.array, 2).val[2] - axes(era.array, 2).val[1]
    map(x, y, years) do lat, long, cal_yr
        if cal_yr < startyr || cal_yr > endyr
            return fill(NaN, 12) .* unit(era.array[1, 1, 1])
        else
            time = cal_yr * 1year
            return Array(era.array[(lat - thisstep1 / 2) .. (lat + thisstep1 / 2),
                                   (long - thisstep2 / 2) .. (long + thisstep2 / 2),
                                   time .. (time + 11month)][1, 1, :])
        end
    end
end

function extractvalues(lat::Union{Missing, typeof(1.0°)},
                       long::Union{Missing, typeof(1.0°)}, ref::Reference)
    if any(ismissing.([lat, long]))
        return NaN .* unit(ref.array[1, 1])
    else
        -180.0° ≤ lat ≤ 180.0° || error("Latitude coordinate is out of bounds")
        -90.0° < long < 90.0° || error("Longitude coordinate is out of bounds")
        thisstep1 = AxisArrays.axes(ref.array, 1).val[2] -
                    AxisArrays.axes(ref.array, 1).val[1]
        thisstep2 = AxisArrays.axes(ref.array, 2).val[2] -
                    AxisArrays.axes(ref.array, 2).val[1]
        return ref.array[(long - thisstep1 / 2) .. (long + thisstep1 / 2),
                         (lat - thisstep2 / 2) .. (lat + thisstep2 / 2)][1, 1]
    end
end

function extractvalues(lat::Vector{typeof(1.0°)},
                       long::Vector{typeof(1.0°)},
                       ref::Reference)
    checkbounds(lat, long) || error("Coordinates out of bounds")
    thisstep1 = step(axes(ref.array, 1).val)
    thisstep2 = step(axes(ref.array, 2).val)
    map(x, y) do long, lat
        return ref.array[(long - thisstep1 / 2) .. (long + thisstep1 / 2),
                         (lat - thisstep2 / 2) .. (lat + thisstep2 / 2)][1, 1]
    end
end
