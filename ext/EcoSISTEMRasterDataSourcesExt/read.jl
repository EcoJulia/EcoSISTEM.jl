# SPDX-License-Identifier: LGPL-3.0-or-later

# --- Reading real datasets ------------------------------------------------------------------------
#
# **Only the parts of `src/datasetread.jl` that genuinely need `RasterDataSources` are
# here** - the two entry points that call `getraster` (so the download is theirs), and the two
# directory readers that name a concrete dataset. The other ~870 lines of that file are generic
# raster plumbing (CRS units, extents, arcsecond snapping, block aggregation, stacking) which the
# parent needs for every raster it handles, whatever its provenance, and which stayed put.
#
# `read(::Type{CRUTS}, ...)` is here for a reason worth knowing: **CRU TS has no layer table of its
# own** and borrows `WorldClim{Climate}`'s, so naming that table means naming a `RasterDataSources`
# type. Reading CRU TS therefore needs this extension even though CRU TS is not one of its datasets.

# **Documented on the stub in `src/datasetread.jl`, not here** - and it was briefly
# documented in *both*, which rendered two `sourcecrs` entries with the same HTML anchor. A name that
# can take a parent stub must be documented there and only there.
function EcoSISTEM.sourcecrs(T::Type{<:RDS.RasterDataSource},
                             layers = RDS.layers(T);
                             cut = nothing, scale = nothing, fn = nothing,
                             kw...)
    raw = getraster(T, layers; _getrasterkw(T)..., kw...)
    r = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        return Raster(_firstfile(raw); lazy = true)
    end
    c = Rasters.crs(r)
    return _isblankcrs(c) ? nothing : c
end

"""
    read(T::Type{<:RasterDataSource}, layers = RasterDataSources.layers(T);
         cut = nothing, scale, fn, kw...)

Download (via `getraster`) and read a RasterDataSources layer set into a [`ClimateRaster`](@ref).
`layers` chooses which layers/variables to read (default: all of them); a vector of *distinct
same-unit* layers (e.g. `[:tmin, :tavg, :tmax]`) is read into one `Dim{:layer}`-stacked array
rather than being conflated with a time axis. `cut`, if given, restricts the result to a
`Extents.Extent(Y = (a, b), X = (c, d))` of `°` bounds. `scale`/`fn` coarsen each raster by an
integer block-aggregation factor with reducer `fn` (source-specific defaults - e.g.
`EarthEnv{LandCover}` is aggregated 10×). Any remaining keywords (e.g. `month`) pass through to
`getraster`.
"""
function Base.read(T::Type{<:RDS.RasterDataSource}, layers = RDS.layers(T);
                   cut = nothing, scale = _defaultscale(T), fn = _defaultfn(T),
                   kw...)
    # Merged once and reused, because this is the only point at which *which* months were asked
    # for is known: `getraster` turns them into paths and the months are gone. `_getrasterkw`
    # supplies the source's own default (WorldClim monthly climate is `month = 1:12`), which the
    # caller's `kw` overrides - so this is the effective request either way, and it is what labels
    # the time axis. Without it a partial read is labelled `1:n` and names each slice after the
    # wrong month.
    rasterkw = (; _getrasterkw(T)..., kw...)
    raw = getraster(T, layers; rasterkw...)
    out = _readraw(T, raw; cut = cut, scale = scale, fn = fn,
                   slices = get(rasterkw, :month, nothing))
    return _rescalepublished(T, layers, out)
end

"""
    read(::Type{CRUTS}, dir::String, var_name::String; cut = nothing)

Read every `.tif` in `dir` as a monthly time series of variable `var_name`, returning a
[`ClimateRaster`](@ref)`{CRUTS}` - `CRUTS` names the *source*, and the raster carries the data.
`var_name` uses the same short codes as `WorldClim{Climate}`'s shipped layer table (`tavg`, `wind`,
`prec`, ...), which is where the attached unit comes from - CRU TS has no layer table of its own.
"""
function Base.read(::Type{CRUTS}, dir::AbstractString, var_name::AbstractString;
                   cut = nothing)
    u = layerunit(RDS.WorldClim{RDS.Climate}, var_name)
    return EcoSISTEM._timeseriesraster(CRUTS,
                                       _readmonthlydir(dir, u; cut = cut))
end

"""
    read(T::Type{CHELSA{Climate}}, dir::String, var_name::String; scale = 1, fn = mean, cut = nothing)

Read every `.tif` in `dir` as a monthly time series of variable `var_name` (optionally coarsened by
`scale`/`fn`) and return a `ClimateRaster{CHELSA{Climate}}`.
"""
function Base.read(T::Type{RDS.CHELSA{RDS.Climate}}, dir::AbstractString,
                   var_name::AbstractString; scale = 1, fn = mean,
                   cut = nothing)
    u = layerunit(T, var_name)
    return ClimateRaster(T,
                         _readmonthlydir(dir, u; scale = scale, fn = fn,
                                         cut = cut))
end
