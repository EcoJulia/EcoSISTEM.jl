using IndexedTables
using AxisArrays
using Unitful.DefaultSymbols
using Plots
using Statistics

"""
    convert_coords(i::Int64, width::Int64)
    convert_coords(x::Int64, y::Int64, width::Int64)
Function to convert coordinates from two-dimensional (`x`,`y`) format to one dimension (`i`), or vice versa, using the `width` of the grid. This function can also be applied to arrays of coordinates.
"""
function convert_coords(i::Int64, width::Int64)
  x = ((i - 1) % width) + 1
  y = div((i - 1), width)  + 1
  return (x, y)
end
function convert_coords(i::Array{Int64, 1}, width::Int64)
  x = ((i .- 1) .% width) .+ 1
  y = div.((i .- 1), width)  .+ 1
  return (x, y)
end
function convert_coords(x::Int64, y::Int64, width::Int64)
  i = x + width * (y - 1)
  return i
end

function convert_coords(x::Array{Int64, 1}, y::Array{Int64, 1}, width::Int64)
  i = x .+ (width .* (y .- 1))
  return i
end
"""
    create_reference(gridsize::Float64)

Function to create a reference grid array of type `Reference`.
"""
function create_reference(gridsize::Float64)
    x = 360 * (1/gridsize) + 1
    y = 180 * (1/gridsize) + 1
    gridsize = gridsize * °
    refarray = AxisArray(Array{Int64, 2}(undef, Int(floor(x)), Int(floor(y))),Axis{:longitude}(-180.0°:gridsize:180.0°),Axis{:latitude}(-90.0°:gridsize:90.0°))
    refarray[1:length(refarray)]= collect(1:length(refarray))
    ref = Reference(refarray)
end
"""
    gardenmask(occ::IndexedTables.IndexedTable, gard::IndexedTables.IndexedTable,
     masksize::Float64)

Function to mask botanic garden locations found in `gard` from occurrence
records found in `occ`, with an optional size for the mask, `masksize`.
"""
function gardenmask(occ::IndexedTable, gard::IndexedTable,
     masksize::Float64)
    coords = hcat(select(gard, :Latitude), select(gard, :Longitude))
    ref = create_reference(masksize)

    # Add in refval column to garden info of which reference grid square it falls in
    gard = pushcol(gard, :refval, extractvalues(coords[:, 2] * °, coords[:, 1] * °, ref))

    coords_occ = hcat(select(occ, :decimallatitude), select(occ, :decimallongitude))
    occ = pushcol(occ, :refval, extractvalues(coords_occ[:, 2] * °, coords_occ[:, 1] * °, ref))
    # Use anti-join to filter out those that have the same reference value
    occ = join(occ, gard, how=:anti, lkey=:refval, rkey =:refval)
    occ = popcol(occ, :refval)
    return occ
end
"""
    genus_clean(genus::IndexedTables.IndexedTable)

Function to clean occurrence data of botanic garden information.
"""
function genus_clean(genus::IndexedTable)
    gardens = load("data/garden_table")
    genus = gardenmask(genus, gardens, 0.02)
    return genus
end

"""
    genus_worldclim_average(genus::IndexedTables.IndexedTable)

Function to clean occurrence data of botanic garden information, and
join with worldclim data.
"""
function genus_worldclim_average(genus::IndexedTable,
    worldclim::IndexedTable)
    worldclim_names = [:prec, :srad, :tavg, :tmax, :tmin, :vapr, :wind]
    genus = genus_clean(genus)
    ref = create_reference(1/12)
    coords = hcat(select(genus, :decimallatitude), select(genus, :decimallongitude))
    genus = pushcol(genus, :refval, extractvalues(coords[:, 2] * °,
     coords[:, 1] * °, ref))
    genus = join(genus, worldclim_dat,  how=:left, lkey=:refval, rkey =:refval,
        rselect=Tuple(Symbol.(worldclim_names)))
    genus = popcol(genus, :refval)
    return genus
end
"""
    genus_worldclim_monthly(genus::IndexedTables.IndexedTable)

Function to clean occurrence data of botanic garden information, and
join with monthly worldclim data.
"""
function genus_worldclim_monthly(genus::IndexedTable,
    worldclim::IndexedTable)
    worldclim_names = [:prec, :srad, :tavg, :tmax, :tmin, :vapr, :wind]
    genus = genus_clean(genus)
    ref = create_reference(1/12)
    coords = hcat(select(genus, :decimallatitude), select(genus, :decimallongitude))
    genus = pushcol(genus, :refval, extractvalues(coords[:, 2] * °,
     coords[:, 1] * °, ref))
    genus = join(genus, worldclim_dat,  how=:left, lkey=:refval, rkey =:refval,
        rselect=Tuple(Symbol.(worldclim_names)))
    genus = popcol(genus, :refval)
    genus = filter(p-> !isnull(p.tavg), genus)
    if length(genus)== 0
     len = 0
    else
     len = Int(length(genus)/12)
    end
    months = repmat(collect(1:12), len)
    genus = pushcol(genus, :worldclim_month, months)
    return genus
end

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

function downresolution(era::ERA, rescale::Int64; fn::Function = mean)
    array = downresolution(era.array, rescale, fn)
    return ERA(array)
end
function downresolution(wc::Worldclim, rescale::Int64; fn::Function = mean)
    array = downresolution(wc.array, rescale, fn)
    return Worldclim(array)
end
function downresolution(bc::Bioclim, rescale::Int64; fn::Function = mean)
    array = downresolution(bc.array, rescale, fn)
    return Bioclim(array)
end
function downresolution(aa::AxisArray{T, 3} where T, rescale::Int64, fn::Function)
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
