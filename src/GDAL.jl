#using BinDeps

#@BinDeps.setup

#function version_check(name, handle)
#    fptr = Libdl.dlsym_e(handle, :GDALVersionInfo)
#    if fptr == C_NULL   #lookup failure
#        return false
#    end
#    versionptr = ccall(fptr,Cstring,(Cstring,),"RELEASE_NAME")
#    versionstring = unsafe_string(versionptr)
#    gdalversion = convert(VersionNumber, versionstring)
#    gdalversion >= v"2.1.0"
#end
#Libdl.dlopen("/home/claireh/.julia/v0.6/Conda/deps/usr/lib/libgdal.so", Libdl.RTLD_GLOBAL)
#libgdal = library_dependency("libgdal",
#                             aliases=["gdal","gdal201", "gdal202",
#                             "gdal_w32","gdal_w64","libgdal.20",
#                             "libgdal.a", "libgdal.la",
#                             "libgdal.so", "libgdal.so.20",
#                             "libgdal.so.20.3.2"],
#                             validate=version_check)
#using Conda
#Conda.add_channel("mychannel")
#provides(Conda.Manager, "gdal", libgdal, installed_libpath="/home/claireh/.julia/v0.6/Conda/deps/usr/lib/")
#@BinDeps.install Dict(:libgdal => :libgdal)
 #Only needs to be done once

#Pkg.clone("https://github.com/visr/GDAL.jl.git")
#Pkg.build("GDAL")
#Pkg.clone("https://github.com/yeesian/ArchGDAL.jl.git")
#Pkg.build("ArchGDAL")

using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using AxisArrays
using NetCDF
using Interpolations

import Unitful: °, °C, mm
import ArchGDAL
import Base.read
const AG = ArchGDAL

vardict = Dict("bio" => NaN, "prec" => mm, "srad" => u"kJ"* u"m"^-2 * day^-1,
"tavg" => °C, "tmax" => °C, "tmin" => °C, "vapr" => u"kPa", "wind" => u"m" * u"s"^-1)
unitdict = Dict("K" => K, "m" => m)
"""
    read(f, filename)

Function to read raster file into julia.
"""
function read(f, filename)
    return AG.registerdrivers() do
        AG.read(filename) do dataset
            f(dataset)
        end
    end
end
"""
    searchdir(path,key)

Function to search a directory `path` using a given `key` string.
"""
searchdir(path,key) = filter(x->contains(x,key), readdir(path))
"""
    extractfile(dir::String)

Function to import a selected file from a path string.
"""
function extractfile(file::String)
    txy = [Float64, Int64(1), Int64(1)]

    read(file) do dataset
        #txy[1] = AG.getdatatype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        print(dataset)
    end

    a = Array{txy[1], 2}(txy[2], txy[3])
    read(file) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    lat, long = size(a, 1), size(a, 2);
    step = 360.0° / lat;
    #step = 180.0° / long;

    world = AxisArray(a[:, long:-1:1],
                           Axis{:latitude}((-180.0°+ step):step:(180.0°)),
                           Axis{:longitude}((-90.0°):step:90.0°));

    if txy[1] <: AbstractFloat
        world[world .== world[Axis{:latitude}(0°),
                                             Axis{:longitude}(step)]] *= NaN;
    end;
    world
end

"""
    extractworldclim(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function extractworldclim(dir::String)
    files = map(searchdir(dir, ".tif")) do files
        joinpath(dir, files)
    end
    txy = [Float64, Int32(1), Int32(1)];

    read(files[1]) do dataset
        txy[1] = Float64
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(Int64(txy[2]), Int64(txy[3]), numfiles);
    map(eachindex(files)) do count
    a = Array{txy[1], 2}(txy[2], txy[3]);
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

    if txy[1] <: AbstractFloat
        world[world .== world[Axis{:latitude}(0°),
                                             Axis{:longitude}(0°),
                                             Axis{:time}(1month)]] *= NaN;
    end;

    Worldclim(world)
end

"""
    extractbioclim(dir::String)

Function to extract all raster files from a specified folder directory,
and convert into an axis array.
"""
function extractbioclim(dir::String)
    files = map(searchdir(dir, ".tif")) do files
        joinpath(dir, files)
    end
    txy = [Float32, Int32(1), Int32(1)];

    read(files[1]) do dataset
        txy[1] = AG.getdatatype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        print(dataset)
    end

    numfiles = length(files)
    b = Array{txy[1], 3}(Int64(txy[2]), Int64(txy[3]), numfiles);
    map(eachindex(files)) do count
    a = Array{txy[1], 2}(txy[2], txy[3]);
    read(files[count]) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    b[:, :, count] = a
    end
    lat, long = size(b, 1), size(b, 2);
    step = 180.0° / long;
    unit = 1.0

    world = AxisArray(b[:, long:-1:1, :] * unit,
                           Axis{:latitude}(-180.0°:step:(180.0° - step / 2)),
                           Axis{:longitude}(-90.0°:step:(90.0°-step/2)),
                           Axis{:var}(1:1:numfiles));

    if txy[1] <: AbstractFloat
        world[world .== world[Axis{:latitude}(0°),
                                             Axis{:longitude}(0°),
                                             Axis{:var}(1)]] *= NaN;
    end;
    Bioclim(world)
end
"""
    extractERA(dir::String, param::String, dim::StepRange(typeof(1month)))

Function to extract a certain parameter, `param`, from an ERA netcdf file,
for a certain timerange, `dim`, and convert into an axis array.
"""
function extractERA(dir::String, param::String, dim::Vector{T}) where T<: Unitful.Time
    lat = reverse(ncread(dir, "latitude"))
    lon = ncread(dir, "longitude")
    units = ncgetatt(dir, param, "units")
    units = unitdict[units]
    array = ncread(dir, param)
    array = array * 1.0
    array[array .≈ ncgetatt(dir, param, "_FillValue")] = NaN
    scale_factor = ncgetatt(dir, param, "scale_factor") * units
    add_offset = ncgetatt(dir, param, "add_offset") * units
    array = array .* scale_factor .+ add_offset

    # If temperature param, need to convert from Kelvin
    if typeof(units) <: Unitful.TemperatureUnits
        array = uconvert.(°C, array)
    end
    if any(lon .== 180)
        splitval = find(lon .== 180)[1]
        firstseg = collect((splitval+1):size(array,1))
        secondseg = collect(1:splitval)
        array = array[vcat(firstseg ,secondseg), :, :]
        lon[firstseg] = 180 - lon[firstseg]
        lon = lon[vcat(reverse(firstseg), secondseg)]
    end
    world = AxisArray(array[:, end:-1:1, :], Axis{:longitude}(lon * °), Axis{:latitude}(lat * °),
                        Axis{:time}(collect(dim)))
    ERA(world)
end
function catERA()
end
function extractCERA(dir::String, file::String, params::String)
    filenames = searchdir(dir, file)
    times = collect((1year+1month):1month:10year)
    cera = extractERA(joinpath(dir, filenames[1]),params,
        times)
    for i in 2:12
        times = ifelse(i == 12, collect(((i-1)*120month +1month):1month:((i-1)*120month +1year)),
                    collect(((i-1)*120month +1month):1month:(i*10year)))
        newcera = extractERA(joinpath(dir, filenames[i]),params,
            times)
        cera.array = cat(3, cera.array, newcera.array)
    end
    return cera
end

function findXY(wc::Worldclim, x::typeof(1.0°), y::typeof(1.0°), m::Int64)
    a = wc.array
    itp = interpolate(a, BSpline(Constant()), OnCell())
    sc = scale(itp, a.axes[1].val, a.axes[2].val, 1:12)
    return sc(x, y, m)
end
function findXY(wc::Worldclim, x::typeof(1.0°), y::typeof(1.0°))
    a = wc.array
    itp = interpolate(a, BSpline(Constant()), OnCell())
    sc = scale(itp, a.axes[1].val, a.axes[2].val, 1:12)
    return sc.(x, y, 1:12)
end
"""
    extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
        wc::Worldclim, dim::Unitful.Time)

Function to extract values from a worldclim object, at specified x, y locations and
time period, `dim`.
"""
function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
   wc::Worldclim, dim::Unitful.Time)
   all(x .<= 180.0) && all(x .>= -180.0) ||
   error("X coordinate is out of bounds")
   all(y .< 90.0) && all(y .> -90.0) ||
   error("Y coordinate is out of bounds")
   thisstep = axes(wc.array, 1).val[2] - axes(wc.array, 1).val[1]
   res = map((i, j) -> wc.array[(i-thisstep/2)..(i+thisstep/2),
                              (j-thisstep/2)..(j+thisstep/2),
                              dim][1,1,1], x, y)
   return transpose(hcat(res...))
end
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

Function to extract values from a bioclim object, at specified x, y locations and
time period, `dim`.
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
function extractvalues(x::Vector{typeof(1.0°)},y::Vector{typeof(1.0°)},
    years::Vector{Int64}, era::ERA, startyr::Int64, endyr::Int64)
    all(x .<= 180.0°) && all(x .>= -180.0°) ||
    error("X coordinate is out of bounds")
    all(y .< 90.0°) && all(y .> -90.0°) ||
    error("Y coordinate is out of bounds")
    thisstep1 = axes(era.array, 1).val[2] - axes(era.array, 1).val[1]
    thisstep2 = axes(era.array, 2).val[2] - axes(era.array, 2).val[1]
    map(x , y, years) do lat, lon, yr
        if yr < startyr || yr > endyr
            return repmat([NaN], 12)
        else
            time = (yr - startyr) * 1year
            return Array(era.array[(lat - thisstep1/2)..(lat + thisstep1/2),
                  (lon - thisstep2/2)..(lon + thisstep2/2),
                  time .. (time + 11month)][1,1,:])
        end
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
