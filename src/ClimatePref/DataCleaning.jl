# SPDX-License-Identifier: LGPL-3.0-or-later

using IndexedTables
using AxisArrays
using Unitful.DefaultSymbols
using Statistics

"""
    convert_coords(i::Int64, width::Int64)
    convert_coords(x::Int64, y::Int64, width::Int64)

Function to convert coordinates from two-dimensional (`x`,`y`) format to one dimension (`i`), or vice versa, using the `width` of the grid. This function can also be applied to arrays of coordinates.
"""
function convert_coords end

function convert_coords(i::Int64, width::Int64)
    x = ((i - 1) % width) + 1
    y = div((i - 1), width) + 1
    return (x, y)
end

function convert_coords(i::Vector{Int64}, width::Int64)
    x = ((i .- 1) .% width) .+ 1
    y = div.((i .- 1), width) .+ 1
    return (x, y)
end

function convert_coords(x::Int64, y::Int64, width::Int64)
    i = x + width * (y - 1)
    return i
end

function convert_coords(x::Vector{Int64}, y::Vector{Int64}, width::Int64)
    i = x .+ (width .* (y .- 1))
    return i
end

"""
    create_reference(gridsize::Float64)

Function to create a reference grid array of type `Reference`.
"""
function create_reference(gridsize::Float64)
    x = 360 * (1 / gridsize) + 1
    y = 180 * (1 / gridsize) + 1
    gridsize = gridsize * °
    refarray = AxisArray(Matrix{Int64}(undef, Int(floor(x)), Int(floor(y))),
                         Axis{:longitude}((-180.0°):gridsize:(180.0°)),
                         Axis{:latitude}((-90.0°):gridsize:(90.0°)))
    vec(refarray) .= eachindex(refarray)
    return Reference(refarray)
end

"""
    upresolution(data::Union{ERA, Worldclim_monthly, ClimateRaster}, rescale::Int64; fn)

Function to increase the resolution of a climate dataset, by a factor, `rescale`.
"""
function upresolution end

function upresolution(era::ERA, rescale::Int64)
    array = upresolution(era.array, rescale)
    return ERA(array)
end

function upresolution(wc::Worldclim_monthly, rescale::Int64)
    array = upresolution(wc.array, rescale)
    return Worldclim_monthly(array)
end

function upresolution(bc::ClimateRaster{T, A}, rescale::Int64) where {T, A}
    array = upresolution(bc.array, rescale)
    return ClimateRaster(T, array)
end

function upresolution(aa::AxisArray{T, 3} where {T}, rescale::Int64)
    grid = size(aa)
    grid = (grid[1] .* rescale .- (rescale - 1),
            grid[2] .* rescale .- (rescale - 1), grid[3])
    array = Array{typeof(aa[1]), 3}(undef, grid)
    aa_ax = Base.axes(aa)
    map(aa_ax[3]) do t
        for x in 1:(size(aa, 1) - 1)
            for y in 1:(size(aa, 2) - 1)
                for dx in 0:rescale
                    fx = dx / rescale
                    ix = rescale * x - (rescale - 1)
                    for dy in 0:rescale
                        fy = dy / rescale
                        iy = rescale * y - (rescale - 1)
                        array[ix + dx, iy + dy, t] = aa[x, y, t] * (1 - fx) *
                                                     (1 - fy) +
                                                     aa[x + 1, y, t] * fx *
                                                     (1 - fy) +
                                                     aa[x, y + 1, t] *
                                                     (1 - fx) * fy +
                                                     aa[x + 1, y + 1, t] * fx *
                                                     fy
                    end
                end
            end
        end
    end

    long = aa.axes[1].val
    newlong = range(long[1], long[end], grid[1])
    lat = aa.axes[2].val
    newlat = range(lat[1], lat[end], grid[2])
    return AxisArray(array,
                     Axis{:longitude}(newlong),
                     Axis{:latitude}(newlat),
                     Axis{:time}(aa.axes[3].val))
end

function upresolution(aa::AxisArray{T, 2} where {T}, rescale::Int64)
    grid = size(aa) .* rescale .- (rescale - 1)
    array = Matrix{typeof(aa[1])}(undef, grid)
    for x in 1:(size(aa, 1) - 1)
        for y in 1:(size(aa, 2) - 1)
            for dx in 0:rescale
                fx = dx / rescale
                ix = rescale * x - (rescale - 1)
                for dy in 0:rescale
                    fy = dy / rescale
                    iy = rescale * y - (rescale - 1)
                    array[ix + dx, iy + dy] = aa[x, y] * (1 - fx) * (1 - fy) +
                                              aa[x + 1, y] * fx * (1 - fy) +
                                              aa[x, y + 1] * (1 - fx) * fy +
                                              aa[x + 1, y + 1] * fx * fy
                end
            end
        end
    end

    long = aa.axes[1].val
    newlong = range(long[1], long[end], grid[1])
    lat = aa.axes[2].val
    newlat = range(lat[1], lat[end], grid[2])

    return AxisArray(array,
                     Axis{:longitude}(newlong),
                     Axis{:latitude}(newlat))
end

"""
    downresolution(data::Union{ERA, Worldclim_monthly, ClimateRaster{WorldClim{BioClim}, <: AxisArray}}, rescale::Int64; fn)

Function to decrease the resolution of a climate dataset, by a factor, `rescale`, and aggregation function, `fn`. The aggregation function has a default setting of taking the mean value.
"""
function downresolution end

function downresolution(era::ERA, rescale::Int64; fn::Function = mean)
    array = downresolution(era.array, rescale, fn = fn)
    return ERA(array)
end

function downresolution(wc::Worldclim_monthly, rescale::Int64;
                        fn::Function = mean)
    array = downresolution(wc.array, rescale, fn = fn)
    return Worldclim_monthly(array)
end

function downresolution(bc::ClimateRaster{T}, rescale::Int64;
                        fn::Function = mean) where {T}
    array = downresolution(bc.array, rescale, fn = fn)
    return ClimateRaster(T, array)
end

function downresolution(aa::AxisArray{T, 3} where {T}, rescale::Int64,
                        fn::Function = mean)
    grid = size(aa)
    grid = ceil.(Int64, (grid[1] / rescale, grid[2] / rescale, grid[3]))
    array = Array{typeof(aa[1]), 3}(undef, grid)
    map(1:grid[3]) do t
        for x in Base.axes(array, 1)
            for y in Base.axes(array, 2)
                xcoords = filter(k -> 1 ≤ rescale * (x - 1) + 1 + k ≤
                                      size(aa, 1),
                                 round(Int, -rescale / 2):round(Int,
                                                                rescale / 2))
                ycoords = filter(k -> 1 ≤ rescale * (y - 1) + 1 + k ≤
                                      size(aa, 2),
                                 round(Int, -rescale / 2):round(Int,
                                                                rescale / 2))
                xrange = min(-minimum(xcoords), maximum(xcoords))
                yrange = min(-minimum(ycoords), maximum(ycoords))
                xcoords = (rescale * (x - 1) + 1) .+ ((-xrange):xrange)
                ycoords = (rescale * (y - 1) + 1) .+ ((-yrange):yrange)
                array[x, y, t] = fn(filter(!isnan, aa[xcoords, ycoords, t]))
            end
        end
    end
    long = aa.axes[1].val
    newlong = range(long[1], long[end], grid[1])
    lat = aa.axes[2].val
    newlat = range(lat[1], lat[end], grid[2])
    return AxisArray(array,
                     Axis{:longitude}(newlong),
                     Axis{:latitude}(newlat),
                     Axis{:time}(aa.axes[3].val))
end

function downresolution(aa::AxisArray{T, 2} where {T}, rescale::Int64;
                        fn::Function = mean)
    grid = size(aa)
    grid = ceil.(Int64, (grid[1] / rescale, grid[2] / rescale))
    array = Matrix{typeof(aa[1])}(undef, grid)
    for x in Base.axes(array, 1)
        for y in Base.axes(array, 2)
            xcoords = filter(k -> 1 ≤ rescale * (x - 1) + 1 + k ≤ size(aa, 1),
                             round(Int, -rescale / 2):round(Int, rescale / 2))
            ycoords = filter(k -> 1 ≤ rescale * (y - 1) + 1 + k ≤ size(aa, 2),
                             round(Int, -rescale / 2):round(Int, rescale / 2))
            xrange = min(-minimum(xcoords), maximum(xcoords))
            yrange = min(-minimum(ycoords), maximum(ycoords))
            xcoords = (rescale * (x - 1) + 1) .+ ((-xrange):xrange)
            ycoords = (rescale * (y - 1) + 1) .+ ((-yrange):yrange)
            array[x, y] = fn(filter(!isnan, aa[xcoords, ycoords]))
        end
    end

    long = aa.axes[1].val
    newlong = range(long[1], long[end], grid[1])
    lat = aa.axes[2].val
    newlat = range(lat[1], lat[end], grid[2])
    return AxisArray(array,
                     Axis{:longitude}(newlong),
                     Axis{:latitude}(newlat))
end

"""
    downresolution!(resized_array::Matrix{T}, array::Matrix{T}, rescale::Int64, fn)
    downresolution!(resized_array::Array{T, 3}, array::Matrix{T}, dim::Int64, rescale::Int64, fn)

Function to decrease the resolution of a climate dataset in place, by a factor, `rescale`, and aggregation function, `fn`. The aggregation function has a default setting of taking the mean value.
"""
function downresolution! end

function downresolution!(resized_array::Matrix{T}, array::Matrix{T},
    dim::Int64, rescale::Int64; fn::Function = mean) where {T}
    dim == 1 || error("Accessing invalid 3rd dimension of 2d array")
    return downresolution!(resized_array, array, rescale, fn = fn)
end

function downresolution!(resized_array::Matrix{T}, array::Matrix{T},
                         rescale::Int64; fn::Function = mean) where {T}
    Threads.@threads for i in eachindex(resized_array)
        x, y = convert_coords(i, size(resized_array, 1))
        xcoords = filter(k -> 1 ≤ rescale * (x - 1) + 1 + k ≤ size(array, 1),
                         round(Int, -rescale / 2):round(Int, rescale / 2))
        ycoords = filter(k -> 1 ≤ rescale * (y - 1) + 1 + k ≤ size(array, 2),
                         round(Int, -rescale / 2):round(Int, rescale / 2))
        xrange = min(-minimum(xcoords), maximum(xcoords))
        yrange = min(-minimum(ycoords), maximum(ycoords))
        xcoords = (rescale * (x - 1) + 1) .+ ((-xrange):xrange)
        ycoords = (rescale * (y - 1) + 1) .+ ((-yrange):yrange)
        resized_array[x, y] = fn(filter(!isnan, array[xcoords, ycoords]))
    end
end

function downresolution!(resized_array::Array{T, 3}, array::Matrix{T},
                         dim::Int64, rescale::Int64;
                         fn::Function = mean) where {T}
    new_dims = size(resized_array, 1) * size(resized_array, 2)
    Threads.@threads for i in 1:new_dims
        x, y = convert_coords(i, size(resized_array, 1))
        xcoords = filter(x -> x .<= size(array, 1),
                         (rescale * x - (rescale - 1)):(rescale * x))
        ycoords = filter(y -> y .<= size(array, 2),
                         (rescale * y - (rescale - 1)):(rescale * y))
        resized_array[x, y, dim] = fn(filter(!isnan, array[xcoords, ycoords]))
    end
end

function compressLC(lc::ClimateRaster{T}) where
         {T <: EarthEnv{<:LandCover}}
    newaa = AxisArray(zeros(Int64, size(lc.array, 1), size(lc.array, 2)),
                      AxisArrays.axes(lc.array, 1),
                      AxisArrays.axes(lc.array, 2))
    Threads.@threads for i in Base.axes(lc.array, 1)
        for j in Base.axes(lc.array, 2)
            newaa[i, j] = findmax(lc.array[i, j, :])[2]
        end
    end

    return ClimateRaster(T, newaa)
end
