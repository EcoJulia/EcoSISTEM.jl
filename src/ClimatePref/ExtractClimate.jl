using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using AxisArrays
using NetCDF

import JuliaDB.DIndexedTable

function checkbounds(x::Vector{typeof(1.0°)}, y::Vector{typeof(1.0°)})
    -180.0 .<= x .<= 180.0 || error("X coordinate is out of bounds")
    -90.0 .<= x .<= 90.0 || error("Y coordinate is out of bounds")
end

function calculatestep(dat::D, dim::Int64) where {D <: AbstractClimate}
    return axes(dat.array, dim).val[2] - axes(dat.array, dim).val[1]
end

function _extract(x, y, dat, dim::Unitful.Time)
    res =  map((i, j) -> wc.array[(i-thisstep/2)..(i+thisstep/2),
        (j-thisstep/2)..(j+thisstep/2), dim][1,1,1], x, y)
    return transpose(hcat(res...))
end
function _extract(x, y, dat, dim::Unitful.Time, step)
    res =  map((i, j) -> wc.array[(i-thisstep/2)..(i+thisstep/2),
        (j-thisstep/2)..(j+thisstep/2), dim][1,1,1], x, y)
    return transpose(hcat(res...))
end

"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        wc::Worldclim, dim::Unitful.Time)

Function to extract values from a worldclim object, at specified x, y locations and
time, `dim`.
"""
function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
   wc::Worldclim, dim::Unitful.Time)
   checkbounds(x, y)
   thisstep = calculatestep(wc, 1)
   res = map((i, j) -> wc.array[(i-thisstep/2)..(i+thisstep/2),
                              (j-thisstep/2)..(j+thisstep/2),
                              dim][1,1,1], x, y)
   return transpose(hcat(res...))
end

"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        wc::Worldclim, dim::Unitful.Time)

Function to extract values from a worldclim object, at specified x, y locations and
over a range of times, `dim`.
"""
function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
   wc::Worldclim, dim::StepRange{typeof(1month)})
   all(x .<= 180.0) && all(x .>= -180.0) ||
   error("X coordinate is out of bounds")
   all(y .< 90.0) && all(y .> -90.0) ||
   error("Y coordinate is out of bounds")
   res = map((i, j) -> wc.array[(i-thisstep/2)..(i+thisstep/2),
                              (j-thisstep/2)..(j+thisstep/2),
                              dim], x, y)
   return transpose(hcat(res...))
end
"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        bc::Bioclim, dim::Unitful.Time)

Function to extract values from a bioclim object, at specified x, y locations and time, `dim`.
"""
function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
   bc::Bioclim, val::Int64)
   all(x .<= 180.0) && all(x .>= -180.0) ||
   error("X coordinate is out of bounds")
   all(y .< 90.0) && all(y .> -90.0) ||
   error("Y coordinate is out of bounds")
   thisstep = axes(bc.array, 1).val[2] - axes(bc.array, 1).val[1]
   res = map((i, j) -> bc.array[(i - thisstep/2)..(i + thisstep/2),
                              (j - thisstep/2)..(j + thisstep/2),
                              val], x, y)
   return transpose(hcat(res...))
end
"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        bc::Bioclim, dim::Unitful.Time)

Function to extract values from a bioclim object, at specified x, y locations and over a time period, `dim`.
"""
function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
   bc::Bioclim, vals::StepRange{Int64})
   all(x .<= 180.0) && all(x .>= -180.0) ||
   error("X coordinate is out of bounds")
   all(y .< 90.0) && all(y .> -90.0) ||
   error("Y coordinate is out of bounds")
   thisstep = axes(bc.array, 1).val[2] - axes(bc.array, 1).val[1]
   res = map((i, j) -> bc.array[(i - thisstep/2)..(i + thisstep/2),
                              (j - thisstep/2)..(j + thisstep/2),
                              start(vals)..last(vals)][1,1,:], x, y)
   return transpose(hcat(res...))
end
"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        years::Vector{Int64}, era::ERA, start::Int64)

Function to extract values from an ERA object, at specified x, y locations and
years. Must be given a starting year of the dataset.
"""
function extractvalues(tab::Union{IndexedTable, DIndexedTable}, era::ERA, varname::Symbol)
    vals = map(t -> extractvalues(t.decimallatitude * °, t.decimallongitude * °, t.year, era), tab)
    tab = pushcol(tab, varname, vals)
    return tab
end


function extractvalues(lat::Union{Missing, typeof(1.0°)}, lon::Union{Missing, typeof(1.0°)}, yr::Union{Missing, Int64}, era::ERA)
    startyr = ustrip(uconvert(year, axes(era.array)[3].val[1]))
    endyr = ustrip(uconvert(year, axes(era.array)[3].val[end]))
    if any(ismissing.([lat, lon, yr])) || yr < startyr || yr > endyr
        return fill(NaN, 12) .* unit(era.array[1,1,1])
    else
        lon <= 180.0° && lon >= -180.0° || error("Longitude coordinate is out of bounds")
        lat <= 90.0° && lat >= -90.0° || error("Latitude coordinate is out of bounds")
        thisstep1 = AxisArrays.axes(era.array, 1).val[2] - AxisArrays.axes(era.array, 1).val[1]
        thisstep2 = AxisArrays.axes(era.array, 2).val[2] - AxisArrays.axes(era.array, 2).val[1]
        time = yr * 1year
        lon += 0.375°
        return era.array[(lon - thisstep1/2)..(lon + thisstep1/2),
              (lat - thisstep2/2)..(lat + thisstep2/2),
              time .. (time + 11month)][1,1,:]
    end
end




function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
    years::Vector{Int64}, era::ERA)
    all(x .<= 180.0°) && all(x .>= -180.0°) ||
    error("X coordinate is out of bounds")
    all(y .< 90.0°) && all(y .> -90.0°) ||
    error("Y coordinate is out of bounds")
    startyr = ustrip(uconvert(year, axes(era.array)[3].val[1]))
    endyr = ustrip(uconvert(year, axes(era.array)[3].val[end]))
    thisstep1 = axes(era.array, 1).val[2] - axes(era.array, 1).val[1]
    thisstep2 = axes(era.array, 2).val[2] - axes(era.array, 2).val[1]
    map(x , y, years) do lat, lon, yr
        if yr < startyr || yr > endyr
            return fill(NaN, 12) .* unit(era.array[1,1,1])
        else
            time = yr * 1year
            return Array(era.array[(lat - thisstep1/2)..(lat + thisstep1/2),
                  (lon - thisstep2/2)..(lon + thisstep2/2),
                  time .. (time + 11month)][1,1,:])
        end
    end
end
function extractvalues(tab::Union{IndexedTable, DIndexedTable}, ref::Reference, varname::Symbol)
    vals = map(t -> extractvalues(t.decimallatitude * °, t.decimallongitude * °, ref), tab)
    tab = pushcol(tab, varname, vals)
    return tab
end
function extractvalues(lat::Union{Missing, typeof(1.0°)}, lon::Union{Missing, typeof(1.0°)}, ref::Reference)
    if any(ismissing.([lat, lon]))
        return NaN .* unit(ref.array[1,1])
    else
        lon <= 180.0° && lon >= -180.0° || error("Longitude coordinate is out of bounds")
        lat <= 90.0° && lat >= -90.0° || error("Latitude coordinate is out of bounds")
        thisstep1 = AxisArrays.axes(ref.array, 1).val[2] - AxisArrays.axes(ref.array, 1).val[1]
        thisstep2 = AxisArrays.axes(ref.array, 2).val[2] - AxisArrays.axes(ref.array, 2).val[1]
        return ref.array[(lon - thisstep1/2)..(lon + thisstep1/2),
              (lat - thisstep2/2)..(lat + thisstep2/2)][1,1]
    end
end
function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
    ref::Reference)
    all(x .<= 180.0°) && all(x .>= -180.0°) ||
    error("X coordinate is out of bounds")
    all(y .<= 90.0°) && all(y .>= -90.0°) ||
    error("Y coordinate is out of bounds")
    thisstep1 = step(axes(ref.array, 1).val)
    thisstep2 = step(axes(ref.array, 2).val)
    map(x , y) do lon, lat
        return ref.array[(lon - thisstep1/2)..(lon + thisstep1/2),
              (lat - thisstep2/2)..(lat + thisstep2/2)][1,1]
    end
end
