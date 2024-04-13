using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using NetCDF
using RCall
using HDF5

"""
    processMet(dir::String, param::String, xmin::Float64 = 5513.0, xmax::Float64 = 470513.0, ymin::Float64 = 530301.5, ymax::Float64 = 1220302.0, res::Int64 = 1000)

Function to process Met Office data and transform to same projection / resolution as Scottish demographics.
"""
function processMet(dir::String, param::String, xmin::Float64 = 5513.0,
                    xmax::Float64 = 470513.0, ymin::Float64 = 530301.5,
                    ymax::Float64 = 1220302.0, res::Int64 = 1000)
    # Alter projection of netcdf data and resample to BNG
    @rput xmin xmax ymin ymax res
    @rput dir
    R"library(raster);
    tmpin = raster(dir)
    crs(tmpin) = '+proj=ob_tran +o_proj=longlat +lon_0=357.5 +o_lon_p=0 +o_lat_p=37.5 +a=6371229 +b=6371229 +to_meter=0.0174532925199 +wktext'
    extent(tmpin)[1] = extent(tmpin)[1] - 360
    extent(tmpin)[2] = extent(tmpin)[2] - 360
    ext = extent(xmin, xmax, ymin, ymax)
    proj = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs'
    resol = c(res, res)
    r = raster(ext = ext, resolution = resol, crs = proj)
    reproj = projectRaster(from = tmpin, crs = proj)
    reproj = resample(reproj, r)
    newdat = reproj@data@values
    x = reproj@nrows
    y = reproj@ncols
    "
    @rget newdat
    @rget x y
    # Reshape output and add units
    xs = ((xmin + res / 2):res:xmax) .* m
    ys = (ymin:res:ymax) .* m
    reshaped_dat = reshape(newdat', y, x)
    reshaped_dat = reshaped_dat[:, end:-1:1]
    reshaped_dat = hcat(reshaped_dat, fill(missing, y))
    reshaped_dat = replace(reshaped_dat, missing => NaN)
    units = ncgetatt(dir, param, "units")
    units = unitdict[units]
    uk = AxisArray(reshaped_dat .* units, Axis{:longitude}(xs),
                   Axis{:latitude}(ys))
    return uk
end

"""
    processMet(dir::String, file::String, param::String, dim::Vector{T}) where T <: Unitful.Time

Function to process several files of Met Office data and transform to same projection / resolution as Scottish demographics, and join them together in a 3D AxisArray.
"""
function processMet(dir::String, file::String, param::String,
                    dim::Vector{T}) where {T <: Unitful.Time}
    filenames = searchdir(dir, file)
    newmet = Array{AxisArray, 1}(undef, length(filenames))
    for i in eachindex(filenames)
        newmet[i] = processMet(joinpath(dir, filenames[i]), param)
    end
    catmet = cat(dims = 3, newmet...)
    catmet = AxisArray(catmet, Axis{:longitude}(catmet.axes[1]),
                       Axis{:latitude}(catmet.axes[2]), Axis{:day}(dim))
    return catmet
end

"""
    writeMet(climatearray::AxisArray, param::String, h5fn=pwd())

Function to write 3D AxisArray of Met Office climate data `climatearray` to a HDF5 file at outputfile `h5fn` for a specified parameter name `param`.
"""
function writeMet(climatearray::AxisArray, param::String, h5fn = pwd())
    axes = (x = string.(axisvalues(climatearray)[1]),
            y = string.(axisvalues(climatearray)[2]),
            times = string.(axisvalues(climatearray)[3]),
            units = string.(unit(climatearray[1])))

    # - initialise HDF5 file
    h5open(h5fn, "w") do fid
        # - create group
        group = g_create(fid, "climate")
        attrs(group)["Description"] = string("Contains climate information for each geographic location and time slice.")
        group[param] = ustrip.(climatearray)

        # fill in axes information
        for k in keys(axes)
            group[string(k)] = axes[k]
        end
    end
end
