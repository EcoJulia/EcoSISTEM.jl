# SPDX-License-Identifier: LGPL-3.0-or-later

# --- What RasterDataSources answers about a dataset's files ---------------------------------------
#
# The three hooks whose sole methods need `getraster`: where a dataset's files are (`_fetchfiles`,
# which downloads them on first use), the first of them opened lazily (`_lazysource`), and its CRS
# (`sourcecrs`). Reading is the parent's - `read(::RasterSpec)` in `src/rasters.jl` - and never
# names this package; a spec's files reach it through `_fetchfiles`.

# **Documented on the stub in `src/datasetread.jl`, not here** - and it was briefly
# documented in *both*, which rendered two `sourcecrs` entries with the same HTML anchor. A name that
# can take a parent stub must be documented there and only there.
function EcoSISTEM.sourcecrs(T::Type{<:RDS.RasterDataSource},
                             layers = RDS.layers(T);
                             cut = nothing, scale = nothing, fn = nothing,
                             kw...)
    c = Rasters.crs(EcoSISTEM._lazysource(T, layers; kw...))
    return _isblankcrs(c) ? nothing : c
end

# The first file the source resolves to, opened lazily. `nothing` for `layers` means the whole
# dataset, as it does for a spec without a code. The read options are accepted and ignored, so a
# spec's stored keywords can be splatted in unchanged.
function EcoSISTEM._lazysource(T::Type{<:RDS.RasterDataSource},
                               layers = RDS.layers(T);
                               cut = nothing, scale = nothing, fn = nothing,
                               kw...)
    raw = getraster(T, something(layers, RDS.layers(T));
                    _getrasterkw(T)..., kw...)
    return Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        return Raster(_firstfile(raw); lazy = true)
    end
end

# Where a dataset's layers are, through `getraster`, which downloads them on first use. Every file
# it hands back that has no provenance record yet gets one - the cache is this package's own, so
# a file there without a record is one this package fetched before it wrote them - naming the
# dataset, the layer where one file is one layer, and the URL RasterDataSources fetched it from
# where its interface can say.
function EcoSISTEM._fetchfiles(T::Type{<:RDS.RasterDataSource}, code; kw...)
    layers = something(code, RDS.layers(T))
    raw = getraster(T, layers; kw...)
    for (layer, path) in _layerfiles(T, layers, raw)
        isfile(EcoSISTEM._sidecarpath(path)) && continue
        record = EcoSISTEM._baserecord(path, :habitat)
        delete!(record, "fetched")      # not known for a file already in the cache
        url = _rasterurl(T, layer; kw...)
        isnothing(url) || (record["url"] = url)
        merge!(record,
               EcoSISTEM._cataloguerecord(T,
                                          isnothing(layer) ? nothing :
                                          EcoSISTEM.layerinfo(T, layer)))
        EcoSISTEM._writesidecar(path, record)
    end
    return raw
end

# Each path `getraster` returned with the layer it holds, or `nothing` for a file holding several
# or a time series of one: the three shapes `getraster` returns, as `_readraw` reads them.
function _layerfiles(T, layers, raw::Vector{<:NamedTuple})
    return [(name, String(p)) for nt in raw for (name, p) in pairs(nt)]
end

function _layerfiles(T, layers, raw::NamedTuple)
    return [(name, String(p)) for (name, p) in pairs(raw)]
end

function _layerfiles(T, layers, raw::AbstractString)
    return [(layers isa Union{Tuple, AbstractVector} ? nothing : layers,
             String(raw))]
end

function _layerfiles(T, layers, raw)
    return [(nothing, String(p)) for p in EcoSISTEM._filelist(raw)]
end

# The URL RasterDataSources fetches a layer from, as text, or `nothing` where its interface has no
# answer for that dataset and layer (a zip holding every layer answers with the zip's).
function _rasterurl(T, layer; kw...)
    for f in (RDS.rasterurl, RDS.zipurl)
        url = try
            isnothing(layer) ? f(T; kw...) : f(T, layer; kw...)
        catch
            nothing
        end
        isnothing(url) && continue
        url isa AbstractVector && (url = isempty(url) ? nothing : first(url))
        isnothing(url) || return string(url)
    end
    return nothing
end

# A dataset's files for `code` that are already on disk, fetching nothing.
function EcoSISTEM._localfiles(T::Type{<:RDS.RasterDataSource}, code; kw...)
    paths = try
        RDS.rasterpath(T, something(code, RDS.layers(T)); kw...)
    catch
        return String[]
    end
    return filter(isfile, EcoSISTEM._filelist(paths))
end
