using IndexedTables
using AxisArrays
using Unitful.DefaultSymbols
using Statistics

"""
    upresolution(data::Union{ERA, Worldclim, Bioclim}, rescale::Int64)

Function to increase the resolution of a climate dataset, by a factor, `rescale`.
"""
function upresolution(era::ERA, rescale::Int64)
    array = upresolution(era.array, rescale)
    return ERA(array)
end
function upresolution(wc::Worldclim, rescale::Int64)
    array = upresolution(wc.array, rescale)
    return Worldclim(array)
end
function upresolution(bc::Bioclim, rescale::Int64)
    array = upresolution(bc.array, rescale)
    return Bioclim(array)
end
function upresolution(aa::AxisArray{T, 3} where T, rescale::Int64)
    grid = size(aa)
    grid = (grid[1] .* rescale, grid[2] .* rescale, grid[3])
    array = Array{typeof(aa[1]), 3}(undef, grid)
    map(1:grid[3]) do time
        for x in 1:size(aa, 1)
            for y in 1:size(aa, 2)
        array[(rescale*x-(rescale-1)):(rescale*x),
            (rescale*y-(rescale - 1)):(rescale*y), time] .= aa[x, y, time]
            end
        end
    end
    lon = aa.axes[1].val
    smallstep = (lon[2] - lon[1]) / rescale
    if lon[1] == -180°
        newlon = collect(lon[1]:smallstep:(lon[end]+smallstep))
    else
        newlon = collect((lon[1] -smallstep):smallstep:lon[end])
    end
    lat = aa.axes[2].val
    smallstep = (lat[2] - lat[1]) / rescale
    if lat[1] == -90°
        newlat = collect(lat[1]:smallstep:(lat[end]+smallstep))
    else
        newlat = collect((lat[1]-smallstep):smallstep:lat[end])
    end
    return AxisArray(array,
        Axis{:longitude}(newlon),
        Axis{:latitude}(newlat),
        Axis{:time}(aa.axes[3].val))
end
function upresolution(aa::AxisArray{T, 2} where T, rescale::Int64)
    grid = size(aa) .* rescale
    array = Array{typeof(aa[1]), 2}(undef, grid)
    for x in 1:size(aa, 1)
        for y in 1:size(aa, 2)
            array[(rescale*x-(rescale-1)):(rescale*x),
            (rescale*y-(rescale - 1)):(rescale*y)] .= aa[x, y]
        end
    end
    lon = aa.axes[1].val
    smallstep = (lon[2] - lon[1]) / rescale
    if lon[1] == -180°
        newlon = collect(lon[1]:smallstep:(lon[end]+smallstep))
    else
        newlon = collect((lon[1] -smallstep):smallstep:lon[end])
    end
    lat = aa.axes[2].val
    smallstep = (lat[2] - lat[1]) / rescale
    if lat[1] == -90°
        newlat = collect(lat[1]:smallstep:(lat[end]+smallstep))
    else
        newlat = collect((lat[1]-smallstep):smallstep:lat[end])
    end
    return AxisArray(array,
        Axis{:longitude}(newlon),
        Axis{:latitude}(newlat))
end

"""
    downresolution(data::Union{ERA, Worldclim, Bioclim}, rescale::Int64)

Function to decrease the resolution of a climate dataset, by a factor, `rescale`, and aggregation function, `fn`. The aggregation function has a default setting of taking the mean value.
"""

function downresolution(era::ERA, rescale::Int64; fn::F=mean) where {F<:Function}
    array = downresolution(era.array, rescale, fn)
    return ERA(array)
end
function downresolution(wc::Worldclim, rescale::Int64; fn::F=mean) where {F<:Function}
    array = downresolution(wc.array, rescale, fn)
    return Worldclim(array)
end
function downresolution(bc::Bioclim, rescale::Int64; fn::F=mean) where {F<:Function}
    array = downresolution(bc.array, rescale, fn)
    return Bioclim(array)
end
function downresolution(aa::AxisArray{T, 3} where T, rescale::Int64, fn::F) where {F<:Function}
    grid = size(aa)
    grid = ceil.(Int64, (grid[1] ./ rescale, grid[2] ./ rescale, grid[3]))
    array = Array{typeof(aa[1]), 3}(undef, grid)
    map(1:grid[3]) do tm
        for x in 1:size(array, 1)
            for y in 1:size(array, 2)
                xcoords = filter(x -> x .<= size(aa, 1), (rescale*x-(rescale-1)):(rescale*x))
                ycoords = filter(y -> y .<= size(aa, 2), (rescale*y-(rescale - 1)):(rescale*y))
                array[x, y, tm] = fn(filter(!isnan, aa[xcoords, ycoords, tm]))
            end
        end
    end
    lon = aa.axes[1].val
    bigstep = (lon[2] - lon[1]) * rescale
    newlon = collect(lon[1]:bigstep:lon[end])
    lat = aa.axes[2].val
    bigstep = (lat[2] - lat[1]) * rescale
    newlat = collect(lat[1]:bigstep:lat[end])
    return AxisArray(array,
        Axis{:longitude}(newlon),
        Axis{:latitude}(newlat),
        Axis{:time}(aa.axes[3].val))
end
function downresolution(aa::AxisArray{T, 2} where T, rescale::Int64, fn)
    grid = size(aa)
    grid = ceil.(Int64, (grid[1] ./ rescale, grid[2] ./ rescale))
    array = Array{typeof(aa[1]), 2}(undef, grid)
    for x in 1:size(array, 1)
        for y in 1:size(array, 2)
            xcoords = filter(x -> x .<= size(aa, 1), (rescale*x-(rescale-1)):(rescale*x))
            ycoords = filter(y -> y .<= size(aa, 2), (rescale*y-(rescale - 1)):(rescale*y))
            array[x, y] = fn(filter(!isnan, aa[xcoords, ycoords]))
        end
    end
    lon = aa.axes[1].val
    bigstep = (lon[2] - lon[1]) * rescale
    newlon = collect(lon[1]:bigstep:lon[end])
    lat = aa.axes[2].val
    bigstep = (lat[2] - lat[1]) * rescale
    newlat = collect(lat[1]:bigstep:lat[end])
    return AxisArray(array,
        Axis{:longitude}(newlon),
        Axis{:latitude}(newlat))
end

function downresolution!(resized_array::Array{T, 2}, array::Array{T, 2}, rescale::Int64, fn) where T
    grid = size(array)
    grid = ceil.(Int64, (grid[1] ./ rescale, grid[2] ./ rescale))
    Threads.@threads for i in 1:length(resized_array)
        x, y = convert_coords(i, size(resized_array, 2))
        xcoords = filter(x -> x .<= size(array, 1), (rescale*x-(rescale-1)):(rescale*x))
        ycoords = filter(y -> y .<= size(array, 2), (rescale*y-(rescale - 1)):(rescale*y))
        resized_array[x, y] = fn(filter(!isnan, array[xcoords, ycoords]))
    end
end
