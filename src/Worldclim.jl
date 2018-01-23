using DataStructures
#using ArrayViews
using DataFrames
using JLD
using RCall
using AxisArrays
addprocs(6)
using JuliaDB
dat = loadtable("/Users/claireh/Documents/PhD/Data/Worldclim/CSV",
usecache=false)

dat = readtable(
  "/Users/claireh/Documents/PhD/Data/Worldclim/CSV/AnnualMeanTemp.csv",
  header = false)
datarray = DataArray{Float64}(dat)[end:-1:1, :]
avMeanTemp = AxisArray(datarray[2:nrow(dat), 1:(ncol(dat)-1)],
    Axis{:x}(datarray[1:999, 1]),
    Axis{:y}(datarray[end,2:43201]))

using Dagger

load(Dagger.chunk,
"/Users/claireh/Documents/PhD/Data/Worldclim/CSV/AnnualMeanTemp.csv")

using Images
using RCall
using Dagger

R"library(raster)"
file = "~/Documents/PhD/Data/Worldclim/wc2.0_5m_bio/wc2.0_bio_5m_01.tif"
img = rcall(:raster, file)
data = rcall(:getValues, img)
jdata = reshape(rcopy(data), 2160, 4320)
Dagger.save()



@rget vals

img =
Images.load("/Users/claireh/Documents/PhD/Data/Worldclim/wc2.0_30s_tmin/wc2.0_30s_tmin_01.tif")

function closest_index(x::DataArray{Float64, 1}, val::Float64)
    return Int(round((val- x[1]) / stepsize))
end


function findval(x::Float64, y::Float64, array::AxisArray{Float64, 2})
    x < maximum(axes(array)[1].val) && x > minimum(axes(array)[1].val) ||
    error("X coordinate is out of bounds")
    y < maximum(axes(array)[2].val) && y > minimum(axes(array)[2].val) ||
    error("Y coordinate is out of bounds")
    xval = closest_index(axes(array)[1].val, x)
    yval = closest_index(axes(array)[2].val, y)
    return array[xval, yval]
end

function extractvalues(x::Array{Float64, 1}, y::Array{Float64, 1},
     array::AxisArray{Float64, 2})
     return map((i, j) -> findval(i, j, array), x, y)
end
