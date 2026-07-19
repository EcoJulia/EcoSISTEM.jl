# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM: LatLong
using AxisArrays
using NetCDF
using IndexedTables

# Grid step of a climate dataset's `dim`-th coordinate axis.
function calculatestep(dat::AbstractClimate, dim::Int64)
    return AxisArrays.axes(dat.array, dim).val[2] -
           AxisArrays.axes(dat.array, dim).val[1]
end

"""
    extractvalues(locs::AbstractVector{<:LatLong}, wc::ClimateRaster{WorldClim{Climate}},
        dim::Unitful.Time)

Extract values from a worldclim object at the [`LatLong`](@ref) locations `locs` and time `dim`.
"""
function extractvalues(locs::AbstractVector{<:LatLong{typeof(1.0°)}},
                       wc::ClimateRaster{WorldClim{Climate}}, dim::Unitful.Time)
    thisstep = calculatestep(wc, 1)
    res = map(locs) do loc
        return wc.array[(loc.lat - thisstep / 2) .. (loc.lat + thisstep / 2),
                        (loc.long - thisstep / 2) .. (loc.long + thisstep / 2),
                        dim][1, 1, 1]
    end
    return transpose(hcat(res...))
end

"""
    extractvalues(locs::AbstractVector{<:LatLong}, bc::ClimateRaster, val::Int64)

Extract values from a bioclim object at the [`LatLong`](@ref) locations `locs` and variable `val`.
"""
function extractvalues(locs::AbstractVector{<:LatLong{typeof(1.0°)}},
                       bc::ClimateRaster, val::Int64)
    thisstep = AxisArrays.axes(bc.array, 1).val[2] -
               AxisArrays.axes(bc.array, 1).val[1]
    res = map(locs) do loc
        return bc.array[(loc.lat - thisstep / 2) .. (loc.lat + thisstep / 2),
                        (loc.long - thisstep / 2) .. (loc.long + thisstep / 2),
                        val]
    end
    return transpose(hcat(res...))
end

"""
    extractvalues(loc::LatLong, cal_yr::Int64, era::ERA)

Extract the year `cal_yr` of monthly values from an ERA object at the [`LatLong`](@ref) location
`loc`, or a NaN-filled vector if the year is outside the data.
"""
function extractvalues(loc::LatLong{typeof(1.0°)}, cal_yr::Int64, era::ERA)
    startyr = ustrip(uconvert(year, axes(era.array)[3].val[1]))
    endyr = ustrip(uconvert(year, axes(era.array)[3].val[end]))
    if cal_yr < startyr || cal_yr > endyr
        return fill(NaN, 12) .* unit(era.array[1, 1, 1])
    end
    thisstep1 = AxisArrays.axes(era.array, 1).val[2] -
                AxisArrays.axes(era.array, 1).val[1]
    thisstep2 = AxisArrays.axes(era.array, 2).val[2] -
                AxisArrays.axes(era.array, 2).val[1]
    time = cal_yr * year
    long = loc.long + 0.375°
    return era.array[(long - thisstep1 / 2) .. (long + thisstep1 / 2),
                     (loc.lat - thisstep2 / 2) .. (loc.lat + thisstep2 / 2),
                     time .. (time + 11month)][1, 1, :]
end

"""
    extractvalues(loc::LatLong, ref::Reference)

Extract the value from a reference object at the [`LatLong`](@ref) location `loc`.
"""
function extractvalues(loc::LatLong{typeof(1.0°)}, ref::Reference)
    thisstep1 = AxisArrays.axes(ref.array, 1).val[2] -
                AxisArrays.axes(ref.array, 1).val[1]
    thisstep2 = AxisArrays.axes(ref.array, 2).val[2] -
                AxisArrays.axes(ref.array, 2).val[1]
    return ref.array[(loc.long - thisstep1 / 2) .. (loc.long + thisstep1 / 2),
                     (loc.lat - thisstep2 / 2) .. (loc.lat + thisstep2 / 2)][1,
                                                                             1]
end
