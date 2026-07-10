# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using NetCDF
using IndexedTables

function checkbounds(x::Vector{typeof(1.0°)}, y::Vector{typeof(1.0°)})
    all(-180.0 .< x .<= 180.0) || error("X coordinate is out of bounds")
    return all(-90.0 .<= x .<= 90.0) || error("Y coordinate is out of bounds")
end

function calculatestep(dat::D, dim::Int64) where {D <: AbstractClimate}
    return AxisArrays.axes(dat.array, dim).val[2] -
           AxisArrays.axes(dat.array, dim).val[1]
end

function _extract(x, y, dat, dim::Unitful.Time)
    res = map((i, j) -> wc.array[(i - thisstep / 2)..(i + thisstep / 2),
                                 (j - thisstep / 2)..(j + thisstep / 2),
                                 dim][1, 1, 1], x, y)
    return transpose(hcat(res...))
end

function _extract(x, y, dat, dim::Unitful.Time, step)
    res = map((i, j) -> wc.array[(i - thisstep / 2)..(i + thisstep / 2),
                                 (j - thisstep / 2)..(j + thisstep / 2),
                                 dim][1, 1, 1], x, y)
    return transpose(hcat(res...))
end

"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        wc::Worldclim_monthly, dim::Unitful.Time)

Extract values from a worldclim object, at specified x, y locations and time,
`dim`.
"""
function extractvalues(x::Vector{typeof(1.0°)},
                       y::Vector{typeof(1.0°)},
                       wc::Worldclim_monthly,
                       dim::Unitful.Time)
    checkbounds(x, y)
    thisstep = calculatestep(wc, 1)
    res = map((i, j) -> wc.array[(i - thisstep / 2)..(i + thisstep / 2),
                                 (j - thisstep / 2)..(j + thisstep / 2),
                                 dim][1, 1, 1], x, y)
    return transpose(hcat(res...))
end

"""
    extractvalues(x::Vector{typeof(1.0°)}, y::Vector{typeof(1.0°)},
        wc::Worldclim_monthly, dim::StepRange{typeof(1month)})

Extract values from a Worldclim object, at specified x, y locations and over a
range of times, `dim`.
"""
function extractvalues(x::Vector{typeof(1.0°)},
                       y::Vector{typeof(1.0°)},
                       wc::Worldclim_monthly,
                       dim::StepRange{typeof(1month)})
    all(x .<= 180.0) && all(x .>= -180.0) ||
        error("X coordinate is out of bounds")
    all(y .< 90.0) && all(y .> -90.0) || error("Y coordinate is out of bounds")
    res = map((i, j) -> wc.array[(i - thisstep / 2)..(i + thisstep / 2),
                                 (j - thisstep / 2)..(j + thisstep / 2),
                                 dim], x, y)
    return transpose(hcat(res...))
end

"""
    extractvalues(x::Vector{typeof(1.0°)}, y::Vector{typeof(1.0°)},
        bc::Worldclim_bioclim, val::Int64)

Extract values from a bioclim object, at specified x, y locations and bioclim
variable index, `val`.
"""
function extractvalues(x::Vector{typeof(1.0°)},
                       y::Vector{typeof(1.0°)},
                       bc::Worldclim_bioclim,
                       val::Int64)
    all(x .<= 180.0) && all(x .>= -180.0) ||
        error("X coordinate is out of bounds")
    all(y .< 90.0) && all(y .> -90.0) || error("Y coordinate is out of bounds")
    thisstep = axes(bc.array, 1).val[2] - axes(bc.array, 1).val[1]
    res = map((i, j) -> bc.array[(i - thisstep / 2)..(i + thisstep / 2),
                                 (j - thisstep / 2)..(j + thisstep / 2),
                                 val], x, y)
    return transpose(hcat(res...))
end

"""
    extractvalues(x::Vector{typeof(1.0°)}, y::Vector{typeof(1.0°)},
        bc::Worldclim_bioclim, vals::StepRange{Int64})

Extract values from a bioclim object, at specified x, y locations and over a
range of bioclim variable indices, `vals`.
"""
function extractvalues(x::Vector{typeof(1.0°)},
                       y::Vector{typeof(1.0°)},
                       bc::Worldclim_bioclim,
                       vals::StepRange{Int64})
    all(x .<= 180.0) && all(x .>= -180.0) ||
        error("X coordinate is out of bounds")
    all(y .< 90.0) && all(y .> -90.0) || error("Y coordinate is out of bounds")
    thisstep = axes(bc.array, 1).val[2] - axes(bc.array, 1).val[1]
    res = map((i, j) -> bc.array[(i - thisstep / 2)..(i + thisstep / 2),
                                 (j - thisstep / 2)..(j + thisstep / 2),
                                 start(vals)..last(vals)][1, 1, :], x, y)
    return transpose(hcat(res...))
end

"""
    extractvalues(lat::Union{Missing, typeof(1.0°)}, lon::Union{Missing, typeof(1.0°)},
        yr::Union{Missing, Int64}, era::ERA)

Extract a 12-month ERA time series for a single location (`lat`, `lon`) and year
`yr`. Returns a vector of 12 values with appropriate units, or a vector of `NaN`
if any of the inputs are `missing` or out of the data range.
"""
function extractvalues(lat::Union{Missing, typeof(1.0°)},
                       lon::Union{Missing, typeof(1.0°)},
                       yr::Union{Missing, Int64},
                       era::ERA)
    startyr = ustrip(uconvert(year, axes(era.array)[3].val[1]))
    endyr = ustrip(uconvert(year, axes(era.array)[3].val[end]))
    if any(ismissing.([lat, lon, yr])) || yr < startyr || yr > endyr
        return fill(NaN, 12) .* unit(era.array[1, 1, 1])
    else
        lon <= 180.0° && lon >= -180.0° ||
            error("Longitude coordinate is out of bounds")
        lat <= 90.0° && lat >= -90.0° ||
            error("Latitude coordinate is out of bounds")
        thisstep1 = AxisArrays.axes(era.array, 1).val[2] -
                    AxisArrays.axes(era.array, 1).val[1]
        thisstep2 = AxisArrays.axes(era.array, 2).val[2] -
                    AxisArrays.axes(era.array, 2).val[1]
        time = yr * 1year
        lon += 0.375°
        return era.array[(lon - thisstep1 / 2)..(lon + thisstep1 / 2),
                         (lat - thisstep2 / 2)..(lat + thisstep2 / 2),
                         time..(time + 11month)][1, 1, :]
    end
end

"""
    extractvalues(x::Vector{typeof(1.0°)}, y::Vector{typeof(1.0°)},
        years::Vector{Int64}, era::ERA)

Extract 12-month ERA time series for vectors of locations `x`, `y` and
corresponding years. Returns `NaN` for any location-year combination outside the
data range.
"""
function extractvalues(x::Vector{typeof(1.0°)},
                       y::Vector{typeof(1.0°)},
                       years::Vector{Int64},
                       era::ERA)
    all(x .<= 180.0°) && all(x .>= -180.0°) ||
        error("X coordinate is out of bounds")
    all(y .< 90.0°) && all(y .> -90.0°) ||
        error("Y coordinate is out of bounds")
    startyr = ustrip(uconvert(year, axes(era.array)[3].val[1]))
    endyr = ustrip(uconvert(year, axes(era.array)[3].val[end]))
    thisstep1 = axes(era.array, 1).val[2] - axes(era.array, 1).val[1]
    thisstep2 = axes(era.array, 2).val[2] - axes(era.array, 2).val[1]
    map(x, y, years) do lat, lon, yr
        if yr < startyr || yr > endyr
            return fill(NaN, 12) .* unit(era.array[1, 1, 1])
        else
            time = yr * 1year
            return Array(era.array[(lat - thisstep1 / 2)..(lat + thisstep1 / 2),
                                   (lon - thisstep2 / 2)..(lon + thisstep2 / 2),
                                   time..(time + 11month)][1, 1, :])
        end
    end
end

"""
    extractvalues(lat::Union{Missing, typeof(1.0°)}, lon::Union{Missing, typeof(1.0°)},
        ref::Reference)

Extract a single value from a [`Reference`](@ref) grid at location (`lat`,
`lon`). Returns `NaN` if either coordinate is `missing`.
"""
function extractvalues(lat::Union{Missing, typeof(1.0°)},
                       lon::Union{Missing, typeof(1.0°)},
                       ref::Reference)
    if any(ismissing.([lat, lon]))
        return NaN .* unit(ref.array[1, 1])
    else
        lon <= 180.0° && lon >= -180.0° ||
            error("Longitude coordinate is out of bounds")
        lat <= 90.0° && lat >= -90.0° ||
            error("Latitude coordinate is out of bounds")
        thisstep1 = AxisArrays.axes(ref.array, 1).val[2] -
                    AxisArrays.axes(ref.array, 1).val[1]
        thisstep2 = AxisArrays.axes(ref.array, 2).val[2] -
                    AxisArrays.axes(ref.array, 2).val[1]
        return ref.array[(lon - thisstep1 / 2)..(lon + thisstep1 / 2),
                         (lat - thisstep2 / 2)..(lat + thisstep2 / 2)][1, 1]
    end
end

"""
    extractvalues(x::Vector{typeof(1.0°)}, y::Vector{typeof(1.0°)}, ref::Reference)

Extract values from a [`Reference`](@ref) grid at vectors of locations `x`, `y`.
"""
function extractvalues(x::Vector{typeof(1.0°)}, y::Vector{typeof(1.0°)},
                       ref::Reference)
    all(x .<= 180.0°) && all(x .>= -180.0°) ||
        error("X coordinate is out of bounds")
    all(y .<= 90.0°) && all(y .>= -90.0°) ||
        error("Y coordinate is out of bounds")
    thisstep1 = step(axes(ref.array, 1).val)
    thisstep2 = step(axes(ref.array, 2).val)
    map(x, y) do lon, lat
        return ref.array[(lon - thisstep1 / 2)..(lon + thisstep1 / 2),
                         (lat - thisstep2 / 2)..(lat + thisstep2 / 2)][1, 1]
    end
end
