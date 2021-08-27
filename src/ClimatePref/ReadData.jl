using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using NetCDF

import Unitful.°, Unitful.°C, Unitful.mm
import ArchGDAL
import Base.read
const AG = ArchGDAL

vardict = Dict("bio" => NaN, "prec" => mm, "srad" => u"kJ"* u"m"^-2 * day^-1, "tavg" => K, "tmax" => K, "tmin" => K, "vapr" => u"kPa", "wind" => u"m" * u"s"^-1)
unitdict = Dict("K" => K, "m" => m, "J m**-2" => J/m^2, "m**3 m**-3" => m^3)
biodict = Dict(zip(1:19, [fill(K, 11); fill(kg/m^2, 8)]))
"""
    read(f, filename)

Function to read raster file into julia.
"""
function read(f, filename)
    return AG.environment() do
        AG.read(filename) do dataset
            f(dataset)
        end
    end
end

"""
    searchdir(path,key)

Function to search a directory `path` using a given `key` string.
"""
searchdir(path,key) = filter(x->occursin(key, x), readdir(path))

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
        #txy[1] = AG.pixeltype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        print(dataset)
    end

    a = Array{txy[1], 2}(undef, txy[2], txy[3])
    read(file) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    lat, long = size(a, 1), size(a, 2);
    step1 = (xmax - xmin) / lat;
    step2 = (ymax - ymin) / long;

    world = AxisArray(a[:, long:-1:1],
                           Axis{:latitude}((xmin + step1):step1:(xmax)),
                           Axis{:longitude}((ymin + step2):step2:ymax));

    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN;
    end;
    world
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
        joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1)];

    read(files[1]) do dataset
        txy[1] = AG.pixeltype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, Int64(txy[2]), Int64(txy[3]), numfiles);
    map(eachindex(files)) do count
    a = Array{txy[1], 2}(undef, txy[2], txy[3]);
    read(files[count]) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    b[:, :, count] = a
    end
    lat, long = size(b, 1), size(b, 2);
    variable = split(dir, "wc2.0_5m_")[2]
    unit = vardict[variable]
    step1 = (xmax - xmin) / lat;
    step2 = (ymax - ymin) / long;

    world = AxisArray(b[:, long:-1:1, :] * unit,
                                            Axis{:latitude}((xmin + step1):step1:(xmax)),
                                            Axis{:longitude}((ymin + step2):step2:ymax),
                                            Axis{:time}(1month:1month:12month));
    if unit == K
        # bugfix
        world .+= uconvert(K, 0.0°C)
    end
    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN;
    end;

    Worldclim(world)
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
        joinpath(dir, files)
    end
    txy = [Float32, Int32(1), Int32(1), Float64(1), ""];

    read(files[1]) do dataset
        txy[1] = AG.pixeltype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        txy[5] = AG.getunittype(AG.getband(dataset, 1))
        print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, Int64(txy[2]), Int64(txy[3]), numfiles);
    map(eachindex(files)) do count
    a = Array{txy[1], 2}(undef, txy[2], txy[3]);
    read(files[count]) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    b[:, :, count] = a
    end
    lat, long = size(b, 1), size(b, 2);
    step1 = (xmax - xmin) / lat;
    step2 = (ymax - ymin) / long;
    unit = 1.0
    world = AxisArray(b[:, long:-1:1, :] * unit,
                            Axis{:latitude}((xmin + step1):step1:(xmax)),
                            Axis{:longitude}((ymin + step2):step2:ymax),
                            Axis{:var}(1:1:numfiles));
    if txy[1] <: AbstractFloat  
        world[isapprox.(world, txy[4])] *= NaN;
    end;
    Bioclim(world)
end

"""
    readERA(dir::String, param::String, dim::StepRange(typeof(1month)))

Function to extract a certain parameter, `param`, from an ERA netcdf file,
for a certain timerange, `dim`, and convert into an axis array.
"""
function readERA(dir::String, param::String, dim::Vector{T}) where T<: Unitful.Time
    lat = reverse(ncread(dir, "latitude"))
    lon = ncread(dir, "longitude")
    units = ncgetatt(dir, param, "units")
    units = unitdict[units]
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
        firstseg = collect((splitval+1):size(array,1))
        secondseg = collect(1:splitval)
        array = array[vcat(firstseg ,secondseg), :, :]
        lon[firstseg] .= 180 .- lon[firstseg]
        lon = lon[vcat(reverse(firstseg), secondseg)]
    end
    world = AxisArray(array[:, end:-1:1, :], Axis{:longitude}(lon * °), Axis{:latitude}(lat * °),
                        Axis{:time}(collect(dim)))
    ERA(world)
end

"""
    readERA(dir::String, file::String, param::String, dim::Vector{Vector{T}})
        where T<: Unitful.Time

Function to extract a certain parameter, `param`, from a directory, `dir`, containing ERA netcdf files,
for a certain timerange, `dim`, and convert into an axis array.
"""
function readERA(dir::String, file::String, param::String, dim::Vector{Vector{T}}) where T<: Unitful.Time
    filenames = searchdir(dir, file)
    newera = Array{AxisArray, 1}(undef, length(filenames))
    for i in eachindex(filenames)
        newera[i] = readERA(joinpath(dir, filenames[i]), param, dim[i]).array
    end
    catera = cat(dims=3, newera ...)
    return ERA(catera)
end

"""
    readCERA(dir::String, file::String, params::String)

Function to extract a certain parameter, `param`, from an CERA-20C netcdf file,
and convert into an axis array.
"""
function readCERA(dir::String, file::String, params::String)
    filenames = searchdir(dir, file)
    times = collect((1901year+1month):1month:1910year)
    cera = readERA(joinpath(dir, filenames[1]),params,
        times)
    for i in 2:12
        times = 1900year .+ ifelse(i == 12, collect(((i-1)*120month +1month):1month:((i-1)*120month +1year)),
                    collect(((i-1)*120month +1month):1month:(i*10year)))
        newcera = readERA(joinpath(dir, filenames[i]),params,
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
function readCRUTS(dir::String)
    files = map(searchdir(dir, ".tif")) do files
        joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Float64(1)];
    files = files[1:10]
    read(files[1]) do dataset
        txy[1] = Float64
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, Int64(txy[2]), Int64(txy[3]), numfiles);
    map(eachindex(files)) do count
        a = Array{txy[1], 2}(undef, txy[2], txy[3]);
        read(files[count]) do dataset
            bd = AG.getband(dataset, 1);
            AG.read!(bd, a);
        end;
        b[:, :, count] = a
    end
    lat, long = size(b, 1), size(b, 2);
    variable = split(dir, "wc2.0_5m_")[2]
    unit = vardict[variable]
    step = 180.0° / long;

    world = AxisArray(b[:, long:-1:1, :] * unit,
                           Axis{:latitude}((-180.0°+ step):step:180.0°),
                           Axis{:longitude}((-90.0°+step):step:90.0°),
                           Axis{:time}(1month:1month:12month));
    if unit == K
        # bugfix
        world .+= uconvert(K, 0.0°C)
    end
    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN;
    end;

    Worldclim(world)
end

"""
    readCHELSA(dir::String)

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
        joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Float64(1)];
    read(files[1]) do dataset
        txy[1] = Float64
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, ceil(Int64, txy[2]/res),  ceil(Int64, txy[3]/res), numfiles);
    a = Array{txy[1], 2}(undef, txy[2], txy[3]);
    map(eachindex(files)) do count
        read(files[count]) do dataset
            bd = AG.getband(dataset, 1);
            AG.read!(bd, a);
        end;
        downresolution!(b[:, :, count], a, res, fn)
    end
    lat, long = size(b, 1), size(b, 2);
    unit = vardict[var_name]
    step1 = (xmax - xmin) / lat;
    step2 = (ymax - ymin) / long;

    world = AxisArray(b[:, long:-1:1, :] * unit,
                            Axis{:latitude}((xmin + step1):step1:(xmax)),
                            Axis{:longitude}((ymin + step2):step2:ymax),
                            Axis{:time}(1month:1month:35years));
    if unit == K
        # bugfix
        world .+= uconvert(K, 0.0°C)
    end
    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN;
    end;

    CHELSA_monthly(world)
end

function readCHELSA_bioclim(dir::String,
    xmin::Unitful.Quantity{Float64} = -180.0°, 
    xmax::Unitful.Quantity{Float64} = 180.0°,
    ymin::Unitful.Quantity{Float64} = -90.0°, 
    ymax::Unitful.Quantity{Float64} = 90.0°;
    res = 1, fn = mean)
    files = map(searchdir(dir, ".tif")) do files
        joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1), Float64(1)];
    read(files[1]) do dataset
        txy[1] = Float64
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(undef, ceil(Int64, txy[2]/res),  ceil(Int64, txy[3]/res), numfiles);
    a = Array{txy[1], 2}(undef, txy[2], txy[3]);
    map(eachindex(files)) do count
        read(files[count]) do dataset
            bd = AG.getband(dataset, 1);
            AG.read!(bd, a);
        end;
        downresolution!(b[:, :, count], a, res, fn)
    end
    lat, long = size(b, 1), size(b, 2);
    unit = 1.0
    step1 = (xmax - xmin) / lat;
    step2 = (ymax - ymin) / long;

    world = AxisArray(b[:, long:-1:1, :] * unit,
                            Axis{:latitude}((xmin + step1):step1:(xmax)),
                            Axis{:longitude}((ymin + step2):step2:ymax),
                            Axis{:var}(1:1:numfiles));
    if unit == K
        # bugfix
        world .+= uconvert(K, 0.0°C)
    end
    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN;
    end;

    CHELSA_bioclim(world)
end
