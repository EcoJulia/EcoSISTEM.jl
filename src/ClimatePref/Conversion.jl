using AxisArrays
using JuliaDB
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Base.Iterators
import AxisArrays: axes

function worldclim_to_DB(wc::Worldclim)
    gridsize = step(axes(wc.array, 1).val)
    ref = create_reference(ustrip.(gridsize))
    x = collect(axes(wc.array, 1).val)
    y = collect(axes(wc.array, 2).val)
    months = 1:12
    expandedxy = collect(product(x, y, months))
    newx = vcat(map(x-> x[1], expandedxy)...)
    newy = vcat(map(x-> x[2], expandedxy)...)
    months = vcat(map(x-> x[3], expandedxy)...)
    values = wc.array[1:end]
    worldclim_tab = table(newx, newy, months, values,
        names = [:x, :y, :month, :val])
    coords = hcat(select(worldclim_tab, :x), select(worldclim_tab, :y))
    ids = extractvalues(coords[:, 1], coords[:, 2], ref)
    worldclim_tab = pushcol(worldclim_tab, :refval, ids)
    return worldclim_tab
end

function CHELSA_to_DB(ch::CHELSA_monthly)
    gridsize = step(axes(ch.array, 1).val)
    ref = create_reference(ustrip.(gridsize))
    x = collect(axes(ch.array, 1).val)
    y = collect(axes(ch.array, 2).val)
    months = 1:12
    expandedxy = collect(product(x, y, months))
    newx = vcat(map(x-> x[1], expandedxy)...)
    newy = vcat(map(x-> x[2], expandedxy)...)
    months = vcat(map(x-> x[3], expandedxy)...)
    values = ch.array[1:end]
    ch_tab = table(newx, newy, months, values,
        names = [:x, :y, :month, :val])
    coords = hcat(select(ch_tab, :x), select(ch_tab, :y))
    ids = extractvalues(coords[:, 1], coords[:, 2], ref)
    ch_tab = pushcol(ch_tab, :refval, ids)
    return ch_tab
end

function era_to_DB(era::Union{CERA, ERA})
    gridsize = AxisArrays.axes(era.array, 1).val[2] - AxisArrays.axes(era.array, 1).val[1]
    ref = create_reference(Float64(ustrip.(gridsize)))
    x = collect(AxisArrays.axes(era.array, 1).val)
    y = collect(AxisArrays.axes(era.array, 2).val)
    times = collect(AxisArrays.axes(era.array, 3).val)
    yrs = uconvert.(year, times)
    un_years = unique(floor.(Int64, ustrip.(yrs)))
    mnths = 0:11
    expandedxy = collect(product(x, y, un_years, mnths))
    newx = vcat(map(x-> x[1], expandedxy)...)/1.0°
    newy = vcat(map(x-> x[2], expandedxy)...)/1.0°
    newyr =  vcat(map(x-> x[3], expandedxy)...)
    newmonth = vcat(map(x-> x[4], expandedxy)...)
    vals = era.array[1:end]
    era_tab = table(newx * °, newy * °, newmonth .+ 1, newyr, vals, names = [:x, :y, :month, :year, :val])
    coords = hcat(select(era_tab, :x), select(era_tab, :y))
    ids = extractvalues(coords[:, 1], coords[:, 2], ref)
    era_tab = pushcol(era_tab, :refval, ids)
    return era_tab
end

function cera_to_DB(cera::CERA)
    return era_to_DB(cera)
end
