# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM: LatLong
using AxisArrays
using CSV
using RasterDataSources
const RDS = RasterDataSources
using Rasters
using NCDatasets
import Rasters: X, Y, Ti

import Unitful.°, Unitful.°C, Unitful.mm
import ArchGDAL
const AG = ArchGDAL

import Base.read

# Convert a Rasters raster (2-D `X`/`Y`, or 3-D `X`/`Y` with a band/time dimension) into the
# `AxisArray` EcoSISTEM uses internally: dim 1 = `:latitude` (ascending), dim 2 = `:longitude`
# (ascending), and an optional dim 3 = `:time` (`(1:n)·month`) or `:var` (`1:n`). Coordinates
# come from the raster's `Y` (latitude) / `X` (longitude) lookups read from the file metadata —
# no hardcoded global extent — and missing/nodata cells become `NaN`. `unit` is attached
# (Unitful); pass `NoUnits` for a dimensionless layer. This replaces the hand-rolled ArchGDAL
# band loops + fixed `-180°..180°`/`-90°..90°` axes the readers used to build.
#
# NB the geographically-correct convention (`:latitude` = `Y`) is a fix: the old readers put the
# longitude range on the `:latitude` axis (and vice versa) — see the git history / plan.
function _rastertoaxisarray(ras::Rasters.AbstractRaster; unit = NoUnits,
                            thirdaxis::Symbol = :time)
    r = Rasters.replace_missing(ras, NaN)
    # latitude ascending (GeoTIFF `Y` is usually stored north→south)
    (first(Rasters.lookup(r, Y)) < last(Rasters.lookup(r, Y))) ||
        (r = reverse(r; dims = Y))
    # reorder axes to (latitude, longitude[, band/time])
    others = filter(d -> !(d isa X || d isa Y), Rasters.dims(r))
    r = permutedims(r, (Y, X, others...))
    latv = collect(Rasters.lookup(r, Y)) .* °
    lonv = collect(Rasters.lookup(r, X)) .* °
    data = Float64.(Array(r)) .* unit
    ndims(data) == 2 &&
        return AxisArray(data, Axis{:latitude}(latv), Axis{:longitude}(lonv))
    n = size(data, 3)
    third = thirdaxis === :time ? Axis{:time}((1:n) .* month) : Axis{:var}(1:n)
    return AxisArray(data, Axis{:latitude}(latv), Axis{:longitude}(lonv), third)
end

# Combine per-file 2-D lat×long `AxisArray`s into one layer: a single 2-D array for one file,
# else a 3-D `AxisArray` stacked along `third` (an `Axis{:var}` or `Axis{:time}`).
function _stacklayers(aas::AbstractVector, third)
    length(aas) == 1 && return first(aas)
    return AxisArray(cat(getfield.(aas, :data)...; dims = 3),
                     AxisArrays.axes(first(aas))[1],
                     AxisArrays.axes(first(aas))[2],
                     third)
end

# Apply the optional `cut` (a `LatLong` of `°` intervals `.lat`, `.long`) to a lat×long[×z] world.
# With the corrected axis convention `.lat` selects latitude (dim 1) and `.long` longitude (dim 2).
_applycut(world, ::Nothing) = world
function _applycut(world, cut)
    return ndims(world) == 2 ? world[cut.lat, cut.long] :
           world[cut.lat, cut.long, :]
end

# Physical unit for a WorldClim/CRUTS-style file, from the variable name in its path via `VARDICT`
# (e.g. `tavg` → K, `prec` → mm); dimensionless (`1.0`) if none matches.
function _fileunit(file::AbstractString)
    vars = split(file, "_")
    i = findfirst(v -> v ∈ keys(VARDICT), vars)
    return isnothing(i) ? 1.0 : VARDICT[vars[i]]
end

# Read a raster file lazily (so large files are streamed, not loaded whole), optionally
# coarsening it by an integer `scale` factor aggregating each block with `fn` (replaces the
# old read-time `downresolution!`). `_rastertoaxisarray` then materialises the result. The read
# is run under a `NullLogger` to drop Rasters' benign metadata warnings (e.g. a nodata value a
# file declares that its integer eltype can't hold, as in CHELSA's `-99999` UInt16 bands) —
# proper nodata masking for float layers still happens.
function _readraster(f::AbstractString; scale::Integer = 1, fn = mean)
    r = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        return Raster(f; lazy = true)
    end
    return scale > 1 ? Rasters.aggregate(fn, r, scale) : r
end

# --- per-source read traits -------------------------------------------------
# After the Rasters migration the RasterDataSources-backed readers differ in only four small,
# type-keyed ways; one generic `read` (below) drives them all by dispatching these traits on the
# source type `T`. Each has a safe default and a per-source override where it actually differs.

# Kind of third axis stacked over a multi-file layer: bands (`:var`) or a monthly series (`:time`).
_thirdaxis(::Type{<:RDS.RasterDataSource}) = :var
_thirdaxis(::Type{<:WorldClim{Climate}}) = :time

# Physical unit attached to the data: none by default; WorldClim monthly climate infers it from the
# variable name in the filename (via `_fileunit`/`VARDICT`).
_layerunit(::Type{<:RDS.RasterDataSource}, files) = NoUnits
_layerunit(::Type{<:WorldClim{Climate}}, files) = _fileunit(first(files))

# Default read-time block-aggregation factor + reducer (land cover is coarsened 10× by default).
_defaultscale(::Type{<:RDS.RasterDataSource}) = 1
_defaultscale(::Type{<:EarthEnv{<:LandCover}}) = 10
_defaultfn(::Type{<:RDS.RasterDataSource}) = mean
_defaultfn(::Type{<:EarthEnv{<:LandCover}}) = x -> round(mean(x))

# Default keywords forwarded to `getraster` (WorldClim monthly climate needs a month range).
_getrasterkw(::Type{<:RDS.RasterDataSource}) = (;)
_getrasterkw(::Type{<:WorldClim{Climate}}) = (month = 1:12,)

# Build the stacked third axis of a multi-file layer from its kind and length.
function _mkthirdaxis(kind::Symbol, n::Integer)
    return kind === :time ? Axis{:time}((1:n) .* month) : Axis{:var}(1:n)
end

# Normalise `getraster`'s return (a single path, a vector, or a keyed collection of paths) to a
# plain `Vector{String}`.
_filelist(x::AbstractString) = [String(x)]
_filelist(x) = String[String(f) for f in values(x)]

# Read a resolved set of raster file paths into a `ClimateRaster` of source `T`, applying the
# per-source unit, third axis and (for °C-stored Kelvin data) the °C→K shift. Shared by `read` and
# the deprecated `readworldclim`.
function _readsource(T::Type{<:RDS.RasterDataSource}, files::Vector{String};
                     cut = nothing, scale = _defaultscale(T),
                     fn = _defaultfn(T))
    u = _layerunit(T, files)
    aas = map(f -> _rastertoaxisarray(_readraster(f; scale = scale, fn = fn);
                                      unit = u), files)
    world = _stacklayers(aas, _mkthirdaxis(_thirdaxis(T), length(aas)))
    # WorldClim/CHELSA store temperature in °C; the unit is attached as K, so shift to true Kelvin.
    u == K && (world .+= uconvert(K, 0.0°C))
    return ClimateRaster(T, _applycut(world, cut))
end

"""
    read(T::Type{<:RasterDataSource}, layers = RasterDataSources.layers(T);
         cut = nothing, scale, fn, kw...)

Download (via `getraster`) and read a RasterDataSources layer set into a [`ClimateRaster`](@ref).
`layers` chooses which layers/variables to read (default: all of them). `cut`, if given, restricts
the result to a `(lat = a .. b, long = c .. d)` box of `°` intervals. `scale`/`fn` coarsen each
raster by an integer block-aggregation factor with reducer `fn` (source-specific defaults — e.g.
`EarthEnv{LandCover}` is aggregated 10×). Any remaining keywords (e.g. `month`) pass through to
`getraster`.
"""
function Base.read(T::Type{<:RDS.RasterDataSource}, layers = RDS.layers(T);
                   cut = nothing, scale = _defaultscale(T), fn = _defaultfn(T),
                   kw...)
    files = _filelist(getraster(T, layers; _getrasterkw(T)..., kw...))
    return _readsource(T, files; cut = cut, scale = scale, fn = fn)
end

# Snap a `°` coordinate `x` to the nearest multiple of step `r` (a `°` quantity) in direction
# `dirn` (`floor` = down, `ceil` = up). Applied as floor to the low edge and ceil to the high edge,
# this always rounds *outward* so a rounded box encloses the exact one.
_snapout(dirn, x, r) = dirn(uconvert(NoUnits, x / r)) * r

"""
    boundingbox(region::AbstractString; islands = false, round = false)

Return the geographic bounding box of `region` as a [`LatLong`](@ref) of `°` intervals, ready to
pass as the `cut` keyword to [`read`](@ref) and the other raster readers. Boxes are read
from the shipped `data/bounding_boxes.csv` table. `islands = true` selects the island-inclusive
extent (the table's `Islands` coverage) instead of the mainland one. `round`, when given a degree
step (e.g. `round = 5°`), snaps the box *outwards* to the nearest multiple of that step so the
rounded box fully contains the exact one; the default `false` leaves it unrounded.
"""
function boundingbox(region::AbstractString; islands::Bool = false,
                     round = false)
    coverage = islands ? "Islands" : "Mainland"
    table = CSV.File(pkgdir(@__MODULE__, "data", "bounding_boxes.csv"))
    matches = filter(r -> r.Region == region && r.Coverage == coverage, table)
    isempty(matches) &&
        error("No bounding box for region \"$region\" ($coverage) in bounding_boxes.csv")
    row = only(matches)
    south, north = row.South * °, row.North * °
    west, east = row.West * °, row.East * °
    if round !== false
        south, west = _snapout(floor, south, round),
                      _snapout(floor, west, round)
        north, east = _snapout(ceil, north, round), _snapout(ceil, east, round)
    end
    return LatLong(south .. north, west .. east)
end

const VARDICT = Dict("bio" => NaN, "prec" => mm,
                     "srad" => u"kJ" * u"m"^-2 * day^-1, "tavg" => K,
                     "tmax" => K, "tmin" => K, "vapr" => u"kPa",
                     "wind" => u"m" * u"s"^-1)
const UNITDICT = Dict("K" => K, "m" => m, "J m**-2" => J / m^2,
                      "m**3 m**-3" => m^3)

"""
    searchdir(path,key)

Function to search a directory `path` using a given `key` string.
"""
searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))

"""
    readfile(file::String; cut = nothing)

Import a raster file from a path string into an `AxisArray`. `cut`, if given, restricts the result
to a [`LatLong`](@ref) box of `°` intervals (e.g. from `boundingbox`).
"""
function readfile(file::String; cut = nothing, xmin = nothing, xmax = nothing,
                  ymin = nothing, ymax = nothing)
    n = count(!isnothing, (xmin, xmax, ymin, ymax))
    if n == 4
        isnothing(cut) ||
            error("`readfile`: pass either `cut` or the `xmin`/`xmax`/`ymin`/`ymax` extent, not both.")
        Base.depwarn("The `xmin`/`xmax`/`ymin`/`ymax` keywords are deprecated; they are being used " *
                     "as `cut = LatLong(ymin .. ymax, xmin .. xmax)`. Pass `cut` (e.g. from " *
                     "`boundingbox`) instead.", :readfile)
        cut = LatLong(ymin .. ymax, xmin .. xmax)
    elseif n != 0
        error("`readfile` needs all four of `xmin`/`xmax`/`ymin`/`ymax` (deprecated) or none of " *
              "them; got $n.")
    end
    return _applycut(_rastertoaxisarray(_readraster(file)), cut)
end

# Deprecated positional extent form `readfile(file, xmin, xmax, ymin, ymax)`; convert to `cut` and
# forward to the keyword method above.
function readfile(file::String, xmin, xmax, ymin, ymax)
    Base.depwarn("`readfile(file, xmin, xmax, ymin, ymax)` is deprecated; pass " *
                 "`cut = LatLong(ymin .. ymax, xmin .. xmax)` (e.g. from `boundingbox`) " *
                 "instead.", :readfile)
    return readfile(file; cut = LatLong(ymin .. ymax, xmin .. xmax))
end

function readworldclim(::Type{WorldClim{Climate}}, files::Vector{String};
                       cut = nothing)
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

    return ClimateRaster(RDS.WorldClim{RDS.Climate}, world)
end

"""
    readworldclim(T::Type{WorldClim{Climate}}, files; cut = nothing)

Deprecated — use `read(WorldClim{Climate}, layers; …)` instead, which downloads via `getraster`
and reads through the same machinery. Retained to read an already-resolved set of monthly climate
raster file paths.
"""
function readworldclim(T::Type{WorldClim{Climate}}, files; cut = nothing)
    Base.depwarn("`readworldclim` is deprecated; use `read(WorldClim{Climate}, layers; …)`.",
                 :readworldclim)
    return _readsource(T, _filelist(files); cut = cut)
end

# Roll a lat×long×time ERA `AxisArray` from the ERA5 0–360° longitude convention onto (-180, 180],
# reordering the data columns so longitude stays ascending — for consistency with the tif readers
# (WorldClim/CHELSA/EarthEnv), which are already on -180..180. A no-op for data already in range.
function _wraplong180(aa::AxisArray)
    wrapped = mod.(ustrip.(AxisArrays.axisvalues(aa)[2]) .+ 180, 360) .- 180
    perm = sortperm(wrapped)
    return AxisArray(aa.data[:, perm, :], AxisArrays.axes(aa)[1],
                     Axis{:longitude}(wrapped[perm] .* °),
                     AxisArrays.axes(aa)[3])
end

"""
    readERA(file::String, param::String, dim::Vector{<:Unitful.Time}; cut = nothing)

Read variable `param` from a **single** ERA netCDF `file` and return an [`ERA`](@ref). `dim` is the
time coordinate attached to the file's layers (one entry per monthly layer); ERA longitudes are
wrapped onto `(-180°, 180°]`. To read and time-concatenate a whole *directory* of files, use the
four-argument method below.
"""
function readERA(file::String, param::String,
                 dim::Vector{<:Unitful.Time}; cut = nothing)
    # NCDatasets (via Rasters) reads the CF coordinate variables, `units` attribute, and
    # `scale_factor`/`add_offset`/`_FillValue` automatically — replacing the manual `ncread`/scale/
    # offset/mask code — and `_rastertoaxisarray` yields the standard (latitude, longitude, time)
    # layout. The physical unit comes from the variable's `units` attribute via `UNITDICT`.
    # `source = NCDsource()` is forced because CDS/`retrieve_era5` files carry no `.nc` extension,
    # so Rasters' extension-based backend guess would otherwise fall back to GDAL — which reads the
    # data but drops the CF coordinates (→ integer indices) and the `units` attribute.
    ras = Raster(file; name = Symbol(param), source = Rasters.NCDsource())
    u = get(UNITDICT, string(get(Rasters.metadata(ras), "units", "")), NoUnits)
    aa = _rastertoaxisarray(ras; unit = u)
    # ERA5 longitudes run 0–360°; roll them onto (-180, 180] to match the other (tif) readers
    aa = _wraplong180(aa)
    # replace the file's time coordinate with the caller-supplied `dim`
    world = AxisArray(aa.data, AxisArrays.axes(aa)[1], AxisArrays.axes(aa)[2],
                      Axis{:time}(collect(dim)))
    return ERA(_applycut(world, cut))
end

# Read the files in `dir` matching `file` via the single-file `readERA`, one per time-vector in
# `dims`, and concatenate them along time (dim 3). Shared by the directory `readERA`/`readCERA`.
function _readeradir(dir::String, file::String, param::String, dims;
                     cut = nothing)
    filenames = searchdir(dir, file)
    arrays = [readERA(joinpath(dir, filenames[i]), param, dims[i]; cut = cut).array
              for i in eachindex(dims)]
    return cat(arrays...; dims = 3)
end

"""
    readERA(dir::String, file::String, param::String,
            dim::Vector{<:AbstractVector{<:Unitful.Time}}; cut = nothing)

Read variable `param` from **every** ERA netCDF file in directory `dir` whose name contains `file`,
concatenating them along time into one [`ERA`](@ref). `dim` gives the time coordinate *per file* — a
vector of time-vectors, one per matched file (in `readdir` order). This is the multi-file wrapper
over the single-file three-argument method above.
"""
function readERA(dir::String, file::String, param::String,
                 dim::Vector{<:AbstractVector{<:Unitful.Time}}; cut = nothing)
    return ERA(_readeradir(dir, file, param, dim; cut = cut))
end

"""
    readCERA(dir::String, file::String, param::String; cut = nothing)

Read variable `param` from the CERA-20C netCDF files in directory `dir` whose names contain `file`
— one per decade — and concatenate them along time into a single [`CERA`](@ref). The monthly time
coordinate for each decade is generated internally; ERA/CERA longitudes are wrapped onto
`(-180°, 180°]`. Built on the single-file three-argument `readERA` (once per file).
"""
function readCERA(dir::String, file::String, param::String; cut = nothing)
    # one decade-length monthly time vector per file (CERA-20C is archived a decade per file)
    times = [collect((1901year + 1month):(1month):(1910year))]
    for i in 2:12
        push!(times,
              1900year .+ ifelse(i == 12,
                     collect(((i - 1) * 120month + 1month):(1month):((i - 1) * 120month + 1year)),
                     collect(((i - 1) * 120month + 1month):(1month):(i * 10year))))
    end
    return CERA(_readeradir(dir, file, param, times; cut = cut))
end

# Assemble a monthly time series from every `.tif` in `dir`, attaching the unit for `var_name`
# (via `VARDICT`) and shifting °C→K when that unit is Kelvin. Shared by the directory-based
# `readCRUTS`/`readCHELSA_monthly`, which differ only in their wrapper type.
function _readmonthlydir(dir::String, var_name::String; scale = 1, fn = mean,
                         cut = nothing)
    files = joinpath.(dir, searchdir(dir, ".tif"))
    u = VARDICT[var_name]
    aas = map(f -> _rastertoaxisarray(_readraster(f; scale = scale, fn = fn);
                                      unit = u), files)
    world = _stacklayers(aas, _mkthirdaxis(:time, length(aas)))
    u == K && (world .+= uconvert(K, 0.0°C))
    return _applycut(world, cut)
end

"""
    readCRUTS(dir::String, var_name::String; cut = nothing)

Read every `.tif` in `dir` as a monthly time series of variable `var_name` and return a `CRUTS`.
"""
function readCRUTS(dir::String, var_name::String; cut = nothing)
    return CRUTS(_readmonthlydir(dir, var_name; cut = cut))
end

"""
    readCHELSA_monthly(dir::String, var_name::String; scale = 1, fn = mean, cut = nothing)

Read every `.tif` in `dir` as a monthly time series of variable `var_name` (optionally coarsened by
`scale`/`fn`) and return a `ClimateRaster{CHELSA{Climate}}`.
"""
function readCHELSA_monthly(dir::String, var_name::String; scale = 1, fn = mean,
                            cut = nothing)
    return ClimateRaster(RDS.CHELSA{RDS.Climate},
                         _readmonthlydir(dir, var_name; scale = scale, fn = fn,
                                         cut = cut))
end
