using NetCDF
using Unitful
using Unitful.DefaultSymbols
using AxisArrays

ncinfo("/Users/claireh/Downloads/test.nc")

x = ncread("/Users/claireh/Downloads/test.nc", "t2m")

y = x * 1.0

y[y .≈ -32767] = NaN
temps = y .* 0.0009362609303530759K .+ 270.9726153656286K
