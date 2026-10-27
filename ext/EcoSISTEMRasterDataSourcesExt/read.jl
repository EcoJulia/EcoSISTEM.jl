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
    _refusemisnamed(T, something(layers, RDS.layers(T)); kw...)
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
    _refusemisnamed(T, layers; kw...)
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
                                          EcoSISTEM.layerinfo(T, layer), path,
                                          route = :published))
        EcoSISTEM._writesidecar(path, record)
    end
    return raw
end

# A fetch that cannot succeed, refused before it is tried: CHELSA-BIOCLIM+'s four `rsds` files, for
# which RasterDataSources asks CHELSA's server for `CHELSA_rsds_<stat>_1981-2010_V.2.1.tif` where it
# has `CHELSA_rsds_1981-2010_<stat>_V.2.1.tif`, so the download is a 404. A file already in place
# reads as usual. Delete this method once RasterDataSources names them as the server does, which
# `test/ext_EcoSISTEMRasterDataSourcesExt.jl` reports by failing.
_refusemisnamed(T, layers; kw...) = nothing

function _refusemisnamed(T::Type{RDS.CHELSA{RDS.BioClimPlus}}, layers; kw...)
    versioning = (k => v for (k, v) in kw if k in (:version, :patch))
    for layer in (layers isa Union{Tuple, AbstractVector} ? layers : (layers,))
        code = Symbol(string(layer))
        startswith(string(code), "rsds_") || continue
        path = RDS.rasterpath(T, code; versioning...)
        m = match(r"^CHELSA_rsds_([a-z]+)_1981-2010_(V\.2\.\d+)\.tif$",
                  basename(path))
        (isfile(path) || isnothing(m)) && continue
        served = "CHELSA_rsds_1981-2010_$(m[1])_$(m[2]).tif"
        url = dirname(string(RDS.rasterurl(T, code; versioning...))) * "/" *
              served
        error("CHELSA-BIOCLIM+'s `$code` cannot be fetched: RasterDataSources asks CHELSA's " *
              "server for `$(basename(path))`, which it does not have; its file is `$served`. " *
              "Download $url to `$path` to read it.")
    end
    return nothing
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

# Where a dataset's files for `code` are or would be, fetching nothing. `rasterpath` takes only the
# dataset's own code, where `getraster` also takes a layer's name, so each code is resolved as a read
# resolves it; it has no default resolution, so `getraster`'s is supplied; and it names one file a
# call, so a range such as `month = 1:12`, which `getraster` maps over, is taken a value at a time.
function EcoSISTEM._sourcepaths(T::Type{<:RDS.RasterDataSource}, code; kw...)
    return [f.path for f in EcoSISTEM._sourcefiles(T, code; kw...)]
end

# The same files, each with the code it holds and the URL RasterDataSources fetches it from.
function EcoSISTEM._sourcefiles(T::Type{<:RDS.RasterDataSource}, code; kw...)
    codes = EcoSISTEM._preferredcode(T, something(code, collect(RDS.layers(T))))
    return [(path = String(path), code = c, url = _rasterurl(T, c; pathkw...))
            for c in vcat(codes)
            for pathkw in _eachfilekw(_rasterpathkw(T; kw...))
            for path in EcoSISTEM._filelist(RDS.rasterpath(T, c; pathkw...))]
end

# Keywords naming one file each: every combination of the values of those holding several, a range
# or a vector, the rest as they are.
function _eachfilekw(kw)
    ks = collect(keys(kw))
    vs = [v isa AbstractArray ? collect(v) : [v] for v in values(kw)]
    return vec([NamedTuple{Tuple(ks)}(combo)
                for combo in Iterators.product(vs...)])
end

# Read options with `getraster`'s default resolution added, for a dataset that has one and where
# none is given.
function _rasterpathkw(T::Type{<:RDS.RasterDataSource}; kw...)
    (haskey(kw, :res) || !hasmethod(RDS.defres, Tuple{Type{T}})) && return kw
    return (; res = RDS.defres(T), kw...)
end
