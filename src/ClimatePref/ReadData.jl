# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using NetCDF

import Unitful.°, Unitful.°C, Unitful.mm
import ArchGDAL
import Base.read
const AG = ArchGDAL

const VARDICT = Dict("bio" => NaN, "prec" => mm,
                     "srad" => u"kJ" * u"m"^-2 * day^-1, "tavg" => K,
                     "tmax" => K, "tmin" => K, "vapr" => u"kPa",
                     "wind" => u"m" * u"s"^-1)
const UNITDICT = Dict("K" => K, "m" => m, "J m**-2" => J / m^2,
                      "m**3 m**-3" => m^3)
const BIODICT = Dict(zip(1:19, [fill(K, 11); fill(kg / m^2, 8)]))
"""
    read(f, filename)

Function to read raster file into julia.
"""
function read(f, filename)
    return AG.environment() do
        AG.read(filename) do dataset
            return f(dataset)
        end
    end
end

"""
    searchdir(path,key)

Function to search a directory `path` using a given `key` string.
"""
searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))

"""
    readfile(file::String)

Function to import a selected file from a path string.
"""
function readfile(file::String, xmin::Unitful.Quantity{Float64} = -180.0°,
                  xmax::Unitful.Quantity{Float64} = 180.0°,
                  ymin::Unitful.Quantity{Float64} = -90.0°,
                  ymax::Unitful.Quantity{Float64} = 90.0°)
    txy = [Float64, Int64(1), Int64(1), Float64(1)]
    #
    read(file) do dataset
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        return print(dataset)
    end

    a = Matrix{txy[1]}(undef, txy[2], txy[3])
    read(file) do dataset
        bd = AG.getband(dataset, 1)
        return AG.read!(bd, a)
    end
    lat, long = size(a, 1), size(a, 2)
    step_lat = (xmax - xmin) / lat
    step_long = (ymax - ymin) / long

    world = AxisArray(a[:, long:-1:1],
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_long:(ymax - step_long / 2.0)))

    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN
    end
    return world
end

"""
    readworldclim(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function readworldclim(dir::String, xmin::Unitful.Quantity{Float64} = -180.0°,
                       xmax::Unitful.Quantity{Float64} = 180.0°,
                       ymin::Unitful.Quantity{Float64} = -90.0°,
                       ymax::Unitful.Quantity{Float64} = 90.0°)
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Int32(1)]

    read(files[1]) do dataset
        txy[1] = AG.pixeltype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        return print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, Int64(txy[2]), Int64(txy[3]), numfiles)
    map(eachindex(files)) do count
        a = Matrix{txy[1]}(undef, txy[2], txy[3])
        read(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return b[:, :, count] = a
    end
    lat, long = size(b, 1), size(b, 2)
    variables = split(first(files), "_")
    if any(map(v -> v ∈ keys(VARDICT), variables))
        findvar = findfirst(map(v -> v ∈ keys(VARDICT), variables))
        unit = VARDICT[variables[findvar]]
    else
        unit = 1.0
    end
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long

    world = AxisArray(b[:, long:-1:1, :] * unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon / 2.0)),
                      Axis{:time}(1month * (1:numfiles)))
    if unit == K
        # bugfix
        world .+= uconvert(K, 0.0°C)
    end
    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN
    end

    return Worldclim_monthly(world)
end

"""
    readbioclim(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function readbioclim(dir::String, xmin::Unitful.Quantity{Float64} = -180.0°,
                     xmax::Unitful.Quantity{Float64} = 180.0°,
                     ymin::Unitful.Quantity{Float64} = -90.0°,
                     ymax::Unitful.Quantity{Float64} = 90.0°)
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float32, Int32(1), Int32(1), Float64(1), ""]

    read(files[1]) do dataset
        txy[1] = AG.pixeltype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        txy[5] = AG.getunittype(AG.getband(dataset, 1))
        return print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, Int64(txy[2]), Int64(txy[3]), numfiles)
    map(eachindex(files)) do count
        a = Matrix{txy[1]}(undef, txy[2], txy[3])
        read(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return b[:, :, count] = a
    end
    lat, long = size(b, 1), size(b, 2)
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long
    unit = 1.0
    world = AxisArray(b[:, long:-1:1, :] * unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon / 2.0)),
                      Axis{:var}(1:numfiles))
    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN
    end
    return Worldclim_bioclim(world)
end

"""
    readERA(dir::String, param::String, dim::StepRange(typeof(1month)))

Function to extract a certain parameter, `param`, from an ERA netcdf file,
for a certain timerange, `dim`, and convert into an axis array.
"""
function readERA(dir::String, param::String,
                 dim::Vector{T}) where {T <: Unitful.Time}
    lat = reverse(ncread(dir, "latitude"))
    lon = ncread(dir, "longitude")
    units = ncgetatt(dir, param, "units")
    units = UNITDICT[units]
    array = ncread(dir, param)
    array = array * 1.0
    array[array .≈ ncgetatt(dir, param, "_FillValue")] .= NaN
    scale_factor = ncgetatt(dir, param, "scale_factor") * units
    add_offset = ncgetatt(dir, param, "add_offset") * units
    array = array .* scale_factor .+ add_offset

    # If temperature param, need to convert from Kelvin
    #if typeof(units) <: Unitful.TemperatureUnits
    #    array = uconvert.(°C, array)
    #end
    if any(lon .== 180)
        splitval = findall(lon .== 180)[1]
        firstseg = collect((splitval + 1):size(array, 1))
        secondseg = collect(1:splitval)
        array = array[vcat(firstseg, secondseg), :, :]
        lon[firstseg] .= 180 .- lon[firstseg]
        lon = lon[vcat(reverse(firstseg), secondseg)]
    end
    world = AxisArray(array[:, end:-1:1, :], Axis{:longitude}(lon * °),
                      Axis{:latitude}(lat * °),
                      Axis{:time}(collect(dim)))
    return ERA(world)
end

"""
    readERA(dir::String, file::String, param::String, dim::Vector{Vector{T}})
        where T<: Unitful.Time

Function to extract a certain parameter, `param`, from a directory, `dir`, containing ERA netcdf files,
for a certain timerange, `dim`, and convert into an axis array.
"""
function readERA(dir::String, file::String, param::String,
                 dim::Vector{Vector{T}}) where {T <: Unitful.Time}
    filenames = searchdir(dir, file)
    newera = Vector{AxisArray}(undef, length(filenames))
    for i in eachindex(filenames)
        newera[i] = readERA(joinpath(dir, filenames[i]), param, dim[i]).array
    end
    catera = cat(dims = 3, newera...)
    return ERA(catera)
end

"""
    readCERA(dir::String, file::String, params::String)

Function to extract a certain parameter, `param`, from an CERA-20C netcdf file,
and convert into an axis array.
"""
function readCERA(dir::String, file::String, params::String)
    filenames = searchdir(dir, file)
    times = collect((1901year + 1month):(1month):(1910year))
    cera = readERA(joinpath(dir, filenames[1]), params,
                   times)
    for i in 2:12
        times = 1900year .+ ifelse(i == 12,
                       collect(((i - 1) * 120month + 1month):(1month):((i - 1) * 120month + 1year)),
                       collect(((i - 1) * 120month + 1month):(1month):(i * 10year)))
        newcera = readERA(joinpath(dir, filenames[i]), params,
                          times)
        cera.array = cat(dims = 3, cera.array, newcera.array)
    end
    return CERA(cera.array)
end

"""
    readCRUTS(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function readCRUTS(dir::String, var_name::String)
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Float64(1)]

    read(files[1]) do dataset
        txy[1] = Float64
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        return print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, Int64(txy[2]), Int64(txy[3]), numfiles)
    map(eachindex(files)) do count
        a = Matrix{txy[1]}(undef, txy[2], txy[3])
        read(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return b[:, :, count] = a
    end
    lat, long = size(b, 1), size(b, 2)
    unit = VARDICT[var_name]
    step_lat = 360.0° / lat
    step_lon = 180.0° / long

    numfiles = length(files)

    world = AxisArray(b .* unit,
                      Axis{:latitude}((-180.0°):step_lat:(180.0° - step_lat / 2.0)),
                      Axis{:longitude}((-90.0°):step_lon:(90.0° - step_lat / 2.0)),
                      Axis{:time}(1month * (1:numfiles)))
    if unit == K
        # bugfix
        world .+= uconvert(K, 0.0°C)
    end
    if txy[1] <: AbstractFloat && !isnothing(txy[4])
        world[isapprox.(world, txy[4])] *= NaN
    end

    return CRUTS(world)
end

"""
    readCHELSA_monthly(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function readCHELSA_monthly(dir::String, var_name::String,
                            xmin::Unitful.Quantity{Float64} = -180.0°,
                            xmax::Unitful.Quantity{Float64} = 180.0°,
                            ymin::Unitful.Quantity{Float64} = -90.0°,
                            ymax::Unitful.Quantity{Float64} = 90.0°;
                            res = 1, fn = mean)
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Float64(1)]
    read(files[1]) do dataset
        txy[1] = Float64
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        return print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, ceil(Int64, txy[2] / res),
                         ceil(Int64, txy[3] / res), numfiles)
    a = Matrix{txy[1]}(undef, txy[2], txy[3])
    map(eachindex(files)) do count
        read(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return downresolution!(b, a, count, res, fn = fn)
    end
    lat, long = size(b, 1), size(b, 2)
    unit = VARDICT[var_name]
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long

    world = AxisArray(b[:, long:-1:1, :] * unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon)),
                      Axis{:time}((1month):(1month):(numfiles * 1month)))
    if unit == K
        # bugfix
        world .+= uconvert(K, 0.0°C)
    end
    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN
    end

    return CHELSA_monthly(world)
end

function readCHELSA_bioclim(dir::String,
                            xmin::Unitful.Quantity{Float64} = -180.0°,
                            xmax::Unitful.Quantity{Float64} = 180.0°,
                            ymin::Unitful.Quantity{Float64} = -90.0°,
                            ymax::Unitful.Quantity{Float64} = 90.0°;
                            res = 1, fn = mean)
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Float64(1)]
    read(files[1]) do dataset
        txy[1] = Float64
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        return print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, ceil(Int64, txy[2] / res),
                         ceil(Int64, txy[3] / res), numfiles)
    a = Matrix{txy[1]}(undef, txy[2], txy[3])
    map(eachindex(files)) do count
        read(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return downresolution!(b, a, count, res, fn = fn)
    end
    lat, long = size(b, 1), size(b, 2)
    unit = 1.0
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long

    world = AxisArray(b[:, long:-1:1, :] * unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon / 2.0)),
                      Axis{:var}(1:numfiles))
    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN
    end

    return CHELSA_bioclim(world)
end

"""
    readlc(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function readlc(dir::String, xmin::Unitful.Quantity{Float64} = -180.0°,
                xmax::Unitful.Quantity{Float64} = 180.0°,
                ymin::Unitful.Quantity{Float64} = -90.0°,
                ymax::Unitful.Quantity{Float64} = 90.0°;
                res = 10, fn = x -> round(mean(x)))
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float32, Int32(1), Int32(1), Float64(1), ""]

    read(files[1]) do dataset
        txy[1] = AG.pixeltype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        txy[5] = AG.getunittype(AG.getband(dataset, 1))
        return print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, ceil(Int64, txy[2] / res),
                         ceil(Int64, txy[3] / res), numfiles)
    a = Matrix{txy[1]}(undef, txy[2], txy[3])
    map(eachindex(files)) do count
        read(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return downresolution!(b, a, count, res, fn = fn)
    end
    lat, long = size(b, 1), size(b, 2)
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long
    unit = 1.0
    world = AxisArray(b[:, long:-1:1, :] * unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon / 2.0)),
                      Axis{:var}(1:numfiles))
    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN
    end

    return world
end
