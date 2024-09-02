# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using NetCDF
using RasterDataSources
const RDS = RasterDataSources

import Unitful.°, Unitful.°C, Unitful.mm
import ArchGDAL
const AG = ArchGDAL

import Base.read

function Base.read(T::Type{<:RDS.RasterDataSource}; kw...)
    return read(T, RDS.layers(T); kw...)
end

function Base.read(T::Type{WorldClim{BioClim}}, layers;
                   cut = nothing, kw...)
    getraster(T, layers; kw...)
    return readbioclim(T, RDS.rasterpath(T), cut = cut)
end

function Base.read(T::Type{<:EarthEnv{<:LandCover}}, layers; res = 10,
                   fn = x -> round(mean(x)), cut = nothing, kw...)
    getraster(T, layers; kw...)
    return readlc(T, RDS.rasterpath(T), res = res, fn = fn,
                  cut = cut)
end

const VARDICT = Dict("bio" => NaN, "prec" => mm,
                     "srad" => u"kJ" * u"m"^-2 * day^-1, "tavg" => K,
                     "tmax" => K, "tmin" => K, "vapr" => u"kPa",
                     "wind" => u"m" * u"s"^-1)
const UNITDICT = Dict("K" => K, "m" => m, "J m**-2" => J / m^2,
                      "m**3 m**-3" => m^3)
const BIODICT = Dict(zip(1:19, [fill(K, 11); fill(kg / m^2, 8)]))

"""
    readag(f, filename)

Function to read raster file into julia.
"""
function readag(f, filename)
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
function readfile(file::String; xmin::Unitful.Quantity{Float64} = -180.0°,
                  xmax::Unitful.Quantity{Float64} = 180.0°,
                  ymin::Unitful.Quantity{Float64} = -90.0°,
                  ymax::Unitful.Quantity{Float64} = 90.0°, cut = nothing)
    txy = [Float64, Int64(1), Int64(1), Float64(1)]

    readag(file) do dataset
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        return print(dataset)
    end

    a = Matrix{txy[1]}(undef, txy[2], txy[3])
    readag(file) do dataset
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
        @view(world.data[world.data .=== txy[4]]) .*= NaN
        @view(world.data[isapprox.(world.data, txy[4])]) .*= NaN
    end

    if !isnothing(cut)
        world = world[cut.lat, cut.long, :]
    end

    return world
end

"""
    readworldclim(dir::String; cut = nothing)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function readworldclim(dir::String; cut = nothing)
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Int32(1)]

    readag(files[1]) do dataset
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
        readag(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return b[:, :, count] = a
    end

    variables = split(first(files), "_")
    if any(map(v -> v ∈ keys(VARDICT), variables))
        findvar = findfirst(map(v -> v ∈ keys(VARDICT), variables))
        unit = VARDICT[variables[findvar]]
    else
        unit = 1.0
    end

    lat, long = size(b, 1), size(b, 2)
    xmin = -180.0°
    xmax = 180.0°
    ymin = -90.0°
    ymax = 90.0°
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long
    world = AxisArray(b[:, long:-1:1, :] * unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon / 2.0)),
                      Axis{:time}((1:numfiles) * month))
    if unit == K
        # bugfix
        world .+= uconvert(K, 0.0°C)
    end
    if txy[1] <: AbstractFloat
        @view(world.data[world.data .=== txy[4]]) .*= NaN
        @view(world.data[isapprox.(world.data, txy[4])]) .*= NaN
    end

    if !isnothing(cut)
        world = world[cut.lat, cut.long, :]
    end

    return Worldclim_monthly(world)
end

"""
    readbioclim(dir::String; cut = nothing)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function readbioclim(T::Type{WorldClim{BioClim}}, dir::String; cut = nothing)
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float32, Int32(1), Int32(1), Float64(1), ""]

    readag(files[1]) do dataset
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
        readag(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return b[:, :, count] = a
    end
    unit = 1.0

    lat, long = size(b, 1), size(b, 2)
    xmin = -180.0°
    xmax = 180.0°
    ymin = -90.0°
    ymax = 90.0°
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long
    world = AxisArray(b[:, long:-1:1, :] * unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon / 2.0)),
                      Axis{:var}(1:numfiles))
    if txy[1] <: AbstractFloat
        @view(world.data[world.data .=== txy[4]]) .*= NaN
        @view(world.data[isapprox.(world.data, txy[4])]) .*= NaN
    end

    if !isnothing(cut)
        world = world[cut.lat, cut.long, :]
    end

    return ClimateRaster(T, world)
end

"""
    readERA(dir::String, param::String, dim::StepRange(typeof(1month)); cut = nothing)

Function to extract a certain parameter, `param`, from an ERA netcdf file,
for a certain timerange, `dim`, and convert into an axis array.
"""
function readERA(dir::String, param::String,
                 dim::Vector{<:Unitful.Time}; cut = nothing)
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

    if !isnothing(cut)
        world = world[cut.lat, cut.long, :]
    end

    return ERA(world)
end

"""
    readERA(dir::String, file::String, param::String, dim::Vector{Vector{<: Unitful.Time}}; cut = nothing)

Function to extract a certain parameter, `param`, from a directory, `dir`, containing ERA netcdf files,
for a certain timerange, `dim`, and convert into an axis array.
"""
function readERA(dir::String, file::String, param::String,
                 dim::Vector{Vector{<:Unitful.Time}}; cut = nothing)
    filenames = searchdir(dir, file)
    newera = Vector{AxisArray}(undef, length(filenames))
    for i in eachindex(filenames)
        newera[i] = readERA(joinpath(dir, filenames[i]), param, dim[i],
                            cut = cut).array
    end
    catera = cat(dims = 3, newera...)

    return ERA(catera)
end

"""
    readCERA(dir::String, file::String, params::String)

Function to extract a certain parameter, `param`, from an CERA-20C netcdf file,
and convert into an axis array.
"""
function readCERA(dir::String, file::String, params::String; cut = nothing)
    filenames = searchdir(dir, file)
    times = collect((1901year + 1month):(1month):(1910year))
    cera = readERA(joinpath(dir, filenames[1]), params,
                   times)
    for i in 2:12
        times = 1900year .+ ifelse(i == 12,
                       collect(((i - 1) * 120month + 1month):(1month):((i - 1) * 120month + 1year)),
                       collect(((i - 1) * 120month + 1month):(1month):(i * 10year)))
        newcera = readERA(joinpath(dir, filenames[i]), params,
                          times, cut = cut)
        cera.array = cat(dims = 3, cera.array, newcera.array)
    end

    return CERA(cera.array)
end

"""
    readCRUTS(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function readCRUTS(dir::String, var_name::String; cut = nothing)
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Float64(1)]

    readag(files[1]) do dataset
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
        readag(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return b[:, :, count] = a
    end

    numfiles = length(files)
    unit = VARDICT[var_name]

    lat, long = size(b, 1), size(b, 2)
    xmin = -180.0°
    xmax = 180.0°
    ymin = -90.0°
    ymax = 90.0°
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long
    world = AxisArray(b .* unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon / 2.0)),
                      Axis{:time}((1:numfiles) * month))
    if unit == K
        # bugfix
        world .+= uconvert(K, 0.0°C)
    end
    if txy[1] <: AbstractFloat && !isnothing(txy[4])
        @view(world.data[world.data .=== txy[4]]) .*= NaN
        @view(world.data[isapprox.(world.data, txy[4])]) .*= NaN
    end

    if !isnothing(cut)
        world = world[cut.lat, cut.long, :]
    end

    return CRUTS(world)
end

"""
    readCHELSA_monthly(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function readCHELSA_monthly(dir::String, var_name::String; res = 1, fn = mean,
                            cut = nothing)
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Float64(1)]
    readag(files[1]) do dataset
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
        readag(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return downresolution!(b, a, count, res, fn = fn)
    end
    unit = VARDICT[var_name]

    lat, long = size(b, 1), size(b, 2)
    xmin = -180.0°
    xmax = 180.0°
    ymin = -90.0°
    ymax = 90.0°
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long
    world = AxisArray(b[:, long:-1:1, :] * unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon / 2.0)),
                      Axis{:time}((1:numfiles) * month))
    if unit == K
        # bugfix
        world .+= uconvert(K, 0.0°C)
    end
    if txy[1] <: AbstractFloat
        @view(world.data[world.data .=== txy[4]]) .*= NaN
        @view(world.data[isapprox.(world.data, txy[4])]) .*= NaN
    end

    if !isnothing(cut)
        world = world[cut.lat, cut.long, :]
    end

    return CHELSA_monthly(world)
end

function readCHELSA_bioclim(T::Type{CHELSA{BioClim}}, dir::String; res = 1,
                            fn = mean, cut = nothing)
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Float64(1)]
    readag(files[1]) do dataset
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
        readag(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return downresolution!(b, a, count, res, fn = fn)
    end
    unit = 1.0

    lat, long = size(b, 1), size(b, 2)
    xmin = -180.0°
    xmax = 180.0°
    ymin = -90.0°
    ymax = 90.0°
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long
    world = AxisArray(b[:, long:-1:1, :] * unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon / 2.0)),
                      Axis{:var}(1:numfiles))

    if txy[1] <: AbstractFloat
        @view(world.data[world.data .=== txy[4]]) .*= NaN
        @view(world.data[isapprox.(world.data, txy[4])]) .*= NaN
    end

    if !isnothing(cut)
        world = world[cut.lat, cut.long, :]
    end

    return ClimateRaster(T, world)
end

"""
    readlc(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""

function readlc(::Type{T}, dir::String; res = 10, fn = x -> round(mean(x)),
                cut = nothing) where {T <: EarthEnv{<:LandCover}}
    files = map(searchdir(dir, ".tif")) do files
        return joinpath(dir, files)
    end
    txy = [Float32, Int32(1), Int32(1), Float64(1), ""]

    readag(files[1]) do dataset
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
        readag(files[count]) do dataset
            bd = AG.getband(dataset, 1)
            return AG.read!(bd, a)
        end
        return downresolution!(b, a, count, res, fn = fn)
    end
    unit = 1.0

    lat, long = size(b, 1), size(b, 2)
    xmin = -180.0°
    xmax = 180.0°
    ymin = -56.0°
    ymax = 90.0°
    step_lat = (xmax - xmin) / lat
    step_lon = (ymax - ymin) / long
    world = AxisArray(b[:, long:-1:1, :] * unit,
                      Axis{:latitude}(xmin:step_lat:(xmax - step_lat / 2.0)),
                      Axis{:longitude}(ymin:step_lon:(ymax - step_lon / 2.0)),
                      Axis{:var}(1:numfiles))
    if txy[1] <: AbstractFloat
        @view(world.data[world.data .=== txy[4]]) .*= NaN
        @view(world.data[isapprox.(world.data, txy[4])]) .*= NaN
    end

    if !isnothing(cut)
        world = world[cut.lat, cut.long, :]
    end

    return ClimateRaster(T, world)
end
