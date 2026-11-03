# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The two things this package fetches, each implementing `Provenancer`'s `IsFetchable`: a file
# downloaded once into the asset cache (`CachedAsset`) and a request to the Copernicus Climate Data
# Store (`CDSRequest`); the cache's own directory; and the provenance records their fetches write,
# with the catalogue's facts about the file folded in. The download, the lock, the sidecar and the
# cache root are `Provenancer`'s.

"""
    CachedAsset(owner::Type, url::AbstractString; path = nothing)

An immutable descriptor for a file downloaded from `url` on first use: into
`EcoSISTEM.assetdir(owner = owner)` under the URL's own file name, or to `path` where one is given,
which is how a project keeps a download under its own data directory. Nothing is downloaded at
construction - call [`assetpath`](@ref) to get the local path, downloading it first if it is not
there. The download is resumed from whatever an interrupted one left, and a provenance record is
written beside the file (see [`provenance`](@ref)).

# Fields

  - `owner`: the type the download belongs to, which names its own subdirectory of the cache so that
    one type's downloads cannot collide with another's, and the catalogue row its provenance record
    is filled from where the owner is a data source.
  - `url`: where to fetch it from.
  - `path`: where it lives once fetched, or `nothing` for the owner's cache directory.
"""
struct CachedAsset
    owner::Type
    url::String
    path::Union{Nothing, String}
    function CachedAsset(owner::Type, url::AbstractString; path = nothing)
        return new(owner, String(url),
                   isnothing(path) ? nothing : String(path))
    end
end

@traitimpl IsFetchable{CachedAsset}

Provenancer.route(::CachedAsset) = PublishedByProvider()

"""
    CDSRequest(dataset::AbstractString, request::AbstractDict, path::AbstractString;
               region = nothing)

A file fetched from the Copernicus Climate Data Store on first use and kept at `path`: the CDS
`dataset` to ask (`"reanalysis-era5-single-levels-monthly-means"`), the `request` as the CDS
takes it (product type, variable, years, months, area, format), and where the answer is written.
Nothing is fetched at construction; [`assetpath`](@ref) fetches it through `CDSAPI.jl` when
`path` does not yet exist, so a file already downloaded is never asked for twice, and writes a
provenance record beside it holding the request and the CDS job. An entry of a
[`RasterSpec`](@ref)'s `files`. `CDSRequest(ERA, code; years, path)` builds one from a catalogued
layer instead of the CDS's own names.

# Fields

  - `dataset`: the CDS dataset name.
  - `request`: the request body, keys as the CDS documents them, values strings or vectors of
    strings.
  - `path`: where the file lives once fetched.
  - `region`: the named region the request's area was taken from, as the
    [`NaturalEarthSpec`](@ref) naming it, which the file's provenance record keeps beside the
    box; `nothing` where the area was given as a box, or there is none.
"""
struct CDSRequest
    dataset::String
    request::Dict{String, Any}
    path::String
    region::Union{Nothing, String}
    function CDSRequest(dataset::AbstractString, request::AbstractDict,
                        path::AbstractString; region = nothing)
        return new(String(dataset),
                   Dict{String, Any}(string(k) => v for (k, v) in request),
                   String(path), isnothing(region) ? nothing : String(region))
    end
end

@traitimpl IsFetchable{CDSRequest}

Provenancer.route(::CDSRequest) = ExtractedByService()

# == Functions ======================================================================================

"""
    EcoSISTEM.assetdir(mod::Module = EcoSISTEM; owner::Union{Type, Nothing} = nothing)

Path to a subdirectory of EcoSISTEM's asset cache, for storing downloaded data outside the
repository. Creating the directory is part of answering.

The cache is EcoSISTEM's Scratch.jl space, or the directory [`EcoSISTEM.set_assetdir!`](@ref)
recorded, which is read on every call. One cache with a subdirectory per package, rather than a
space each, keeps every download under one lifecycle: a scratch space is created on first use, and
reclaimed by `Pkg.gc()` when EcoSISTEM is removed. `RasterDataSources` is put under it the same way,
by the extension's `__init__` setting `RASTERDATASOURCES_PATH`.

# Arguments

  - `mod`: whose subdirectory to return, defaulting to EcoSISTEM's own.
  - `owner`: a type to nest a further subdirectory under, so each owner's downloads are separate -
    see [`CachedAsset`](@ref).
"""
function assetdir(mod::Module = EcoSISTEM;
                  owner::Union{Type, Nothing} = nothing)
    return Provenancer.assetdir(mod, owner = owner, cache = EcoSISTEM)
end

"""
    EcoSISTEM.set_assetdir!(dir::AbstractString)
    EcoSISTEM.set_assetdir!(nothing)

Keep EcoSISTEM's asset cache - every download, and the memoised aggregates - in `dir` from now on,
rather than in EcoSISTEM's Scratch.jl space, which lives in the first Julia depot and may be on a
disk too small for climate data. `nothing` returns the cache to the scratch space. Return the
directory the cache is now in, made if need be.

The choice is a preference in the `LocalPreferences.toml` of the first environment on the load path
that has EcoSISTEM, so it holds for every later session taking EcoSISTEM from that environment, and
[`EcoSISTEM.assetdir`](@ref) reads it on every call. Files already downloaded stay where they are,
and are fetched again unless moved across. `RasterDataSources`' downloads move with the cache where
the extension put them in it, and stay where `RASTERDATASOURCES_PATH` was set otherwise. `Pkg.gc()`
never reclaims a directory chosen here.

# Arguments

  - `dir`: the directory to keep the cache in, `~` and a relative path expanded.
"""
function set_assetdir!(dir::AbstractString)
    root = Provenancer.set_assetdir!(EcoSISTEM, dir)
    _repointextensions()
    return root
end

function set_assetdir!(::Nothing)
    root = Provenancer.set_assetdir!(EcoSISTEM, nothing)
    _repointextensions()
    return root
end

"""
    assetpath(asset::CachedAsset; verify = false)
    assetpath(request::CDSRequest; verify = false)
    assetpath(relative_path::AbstractString)

Return a local path to a file, downloading it first if it is not there yet, and writing a
provenance record beside it when it does - see [`provenance`](@ref).

# Arguments

  - `asset`: a [`CachedAsset`](@ref), fetched from its `url` into its owner's subdirectory of
    EcoSISTEM's asset cache, or to its own `path`, if the file is not already there. An interrupted
    download is resumed from where it stopped, by this or a later run; a second process asking for
    the same file waits for the first rather than fetching it again.
  - `request`: a [`CDSRequest`](@ref), fetched from the Climate Data Store into its own `path` if
    that file does not already exist. The fetch needs `CDSAPI` loaded and a CDS key in
    `~/.cdsapirc`; a file already present needs neither. A present file whose provenance record
    says it was fetched for a different request - another box, or other months - is refused,
    naming what differs, and left as it is.
  - `verify`: `true` checks a file that is already present against the checksum its provenance
    record holds, erroring on a mismatch - a truncated or replaced file. Off by default, since
    hashing a large file on every read would cost more than the read.
  - `relative_path`: a path resolved under [`assetdir`](@ref) directly, for a cached file whose
    name is already known and which needs no descriptor.
"""
function assetpath(asset::CachedAsset; verify::Bool = false)
    path = something(asset.path,
                     joinpath(assetdir(owner = asset.owner),
                              basename(asset.url)))
    return Provenancer.fetchfile(asset.url, path, verify = verify) do response
        return _assetrecord(asset, path, response)
    end
end

function assetpath(request::CDSRequest; verify::Bool = false)
    if isfile(request.path)
        _checkrequest(request)
        Provenancer.recordpresent(request.path) do
            return _requestrecord(request, request.path, nothing)
        end
        verify && Provenancer.verifyfile(request.path)
        return request.path
    end
    isnothing(Base.get_extension(EcoSISTEM, :EcoSISTEMERAExt)) &&
        error("`$(request.path)` is not there, and fetching it from the Climate Data Store " *
              "needs `using CDSAPI` (and a CDS key in `~/.cdsapirc`).")
    mkpath(dirname(abspath(request.path)))
    Provenancer.fetchlocked(request.path) do
        status = _fetchcds(request)
        isfile(request.path) ||
            error("the Climate Data Store fetch for `$(request.path)` finished without writing it.")
        return Provenancer.writesidecar(request.path,
                                        _requestrecord(request, request.path,
                                                       status))
    end
    return request.path
end

assetpath(relative_path::AbstractString) = joinpath(assetdir(), relative_path)

# A record of a file this package fetched from its catalogue: it carries the fetch's time, which
# a file the user gives or cuts from a source is recorded without, and names a catalogued dataset
# as its source, as `_cataloguerecord` writes it.
function Provenancer.iscatalogued(record::ProvenanceRecord)
    return !isnothing(record.fetched) &&
           !isnothing(_datasetrecord(record.source))
end

# Tell each loaded extension that keeps files under the asset cache that the cache has moved.
function _repointextensions()
    isnothing(Base.get_extension(EcoSISTEM, :EcoSISTEMRasterDataSourcesExt)) ||
        _repointrasterdatasources()
    return nothing
end

# The catalogue's facts about a file of a catalogued source, as the record of every file of it
# carries them: `Provenancer`'s entries for the row's facts, with the file's own title and the layer
# it holds where it is one. Empty where the owner is no source.
function _cataloguerecord(owner::Type, layer, path::AbstractString;
                          route::AbstractFileRoute)
    rec = _datasetrecord(owner)
    isnothing(rec) && return Dict{String, Any}()
    return Provenancer.sourceentries(SourceFacts(rec), route,
                                     title = _filetitle(rec, layer, path),
                                     code = isnothing(layer) ? nothing :
                                            first(layer.aliases))
end

# The provenance record of a file fetched over https: what `Provenancer` records of any download
# - the URLs, the validators, the provider's time - plus the catalogue's facts where the owner is a
# source and the URL is one of its layers' files. `response` is `nothing` for a file recorded when
# already present, which then carries the URL alone.
function _assetrecord(asset::CachedAsset, path::AbstractString, response)
    record = Provenancer.downloadrecord(path, asset.url, response,
                                        _assetrole(asset.owner))
    rec = _datasetrecord(asset.owner)
    isnothing(rec) && return record
    layer = rec.format === :Shapefile ? nothing :
            _layerbyfile(asset.owner, asset.url)
    merge!(record,
           _cataloguerecord(asset.owner, layer, path,
                            route = Provenancer.route(asset)))
    version = _fileversion(rec, path)
    isnothing(version) || (record["version"] = version)
    return record
end

# The part of the system a download owned by `owner` enters by: an outline for a shape spec or a
# named region, and a layer for everything else.
function _assetrole(owner::Type)
    return owner <: AbstractShapeSpec || owner === NaturalEarthLevel ? :region :
           :habitat
end

# The version Natural Earth writes inside each of its zips as `<name>.VERSION.txt`, or `nothing`
# for a zip holding none - read through GDAL's zip filesystem, which opens that one entry.
function _zipversion(path::AbstractString)
    entry = "/vsizip/" * abspath(path) * "/" *
            chopsuffix(basename(path), ".zip") * ".VERSION.txt"
    handle = ArchGDAL.GDAL.vsifopenl(entry, "rb")
    handle == C_NULL && return nothing
    buffer = Vector{UInt8}(undef, 256)
    n = try
        ArchGDAL.GDAL.vsifreadl(buffer, 1, length(buffer), handle)
    finally
        ArchGDAL.GDAL.vsifclosel(handle)
    end
    version = strip(String(buffer[1:n]))
    return isempty(version) ? nothing : String(version)
end

# The provenance record of a file fetched from the Climate Data Store: the dataset asked, the
# request body as sent, the job the CDS ran, the file's own `history` attribute where the fetch
# replaced it, and the catalogue's facts where the dataset is one the catalogue knows (ERA5 monthly
# means are `ERA`'s) and the variable one of its layers. An extract a service made to a request is
# never its source's own.
function _requestrecord(request::CDSRequest, path::AbstractString, status)
    record = Provenancer.baserecord(path, :habitat)
    record["dataset"] = request.dataset
    record["request"] = request.request
    record["primary"] = false
    isnothing(request.region) || (record["region"] = request.region)
    for (key, field) in (("jobID", "job"), ("history", "history"))
        value = status isa AbstractDict ? get(status, key, nothing) : nothing
        isnothing(value) || (record[field] = string(value))
    end
    if request.dataset == _ERA5_MONTHLY
        variables = get(request.request, "variable", String[])
        layer = length(variables) == 1 ?
                _layerbyrequest(ERA, string(only(variables))) : nothing
        merge!(record,
               _cataloguerecord(ERA, layer, path,
                                route = Provenancer.route(request)))
    end
    return record
end

# What a fetched netCDF file's global `history` attribute is replaced by: the attribute records when
# and in which temporary directory the provider converted the file, so two fetches of the same data
# differ in those bytes alone, and with this in its place they are identical.
const _HISTORYNOTE = "Rewritten by EcoSISTEM so that fetches of the same data are identical; " *
                     "the original is in the provenance record."

# Replace the global `history` attribute of the netCDF file at `path` with `_HISTORYNOTE`, padded or
# cut to the original's length so the file keeps its size and layout, and return the original text;
# `nothing`, touching nothing, for a file with no `history` text or one holding the note already.
# Done to a Climate Data Store file before it is renamed into place and its checksum taken.
function _normalisehistory!(path::AbstractString)
    return NCDatasets.Dataset(path, "a") do ds
        history = get(ds.attrib, "history", nothing)
        history isa AbstractString || return nothing
        n = ncodeunits(history)
        note = rpad(first(_HISTORYNOTE, n), n)
        history == note && return nothing
        ds.attrib["history"] = note
        return String(history)
    end
end

# Refuse a present file whose record of ours says it was fetched for another request - a box or
# months changed under one path - rather than serve it as this one. A file with no record, or with
# another program's, passes, presence being all there is to go on.
function _checkrequest(request::CDSRequest)
    sidecar = Provenancer.sidecarpath(request.path)
    isfile(sidecar) || return nothing
    record = TOML.tryparsefile(sidecar)
    (record isa AbstractDict && Provenancer._isourrecord(record) &&
     haskey(record, "request")) || return nothing
    differ = _requestdifferences(record["request"], request.request)
    recorded = get(record, "dataset", request.dataset)
    recorded == request.dataset ||
        pushfirst!(differ,
                   "dataset: recorded $(repr(recorded)), asked $(repr(request.dataset))")
    isempty(differ) && return nothing
    lines = join(differ, "\n  ")
    return error("`$(request.path)` was fetched for a different request, so it is not served as " *
                 "this one:\n  $lines\nGive this request a path of its own, or delete the file " *
                 "and its `.provenance.toml` to fetch it again.")
end

# What two Climate Data Store request bodies differ in, a line per key, empty where they ask for
# the same data: years and months compare as sets, and a single value as a list of one.
function _requestdifferences(recorded::AbstractDict, asked::AbstractDict)
    lines = String[]
    for key in sort!(collect(union(keys(recorded), keys(asked))))
        r = get(recorded, key, nothing)
        a = get(asked, key, nothing)
        _requestvalue(key, r) == _requestvalue(key, a) ||
            push!(lines, "$key: recorded $(repr(r)), asked $(repr(a))")
    end
    return lines
end

# A request value in the form two requests are compared in.
_requestvalue(::Any, ::Nothing) = nothing

function _requestvalue(key, value)
    list = value isa AbstractVector ? collect(value) : [value]
    return key in ("year", "month") ? sort!(unique!(string.(list))) : list
end

# The body of a Climate Data Store request for ERA5 monthly means on single levels: one
# `variable` by its CDS name, the `years` and `months` in order as strings, each month of two
# digits, an `area` as `[north, west, south, east]` in degrees where given, netCDF unarchived. The
# shape the CDS takes and `CDSAPI.retrieve` sends.
function _era5request(variable::AbstractString, years; months = 1:12,
                      area = nothing)
    ms = sort!(unique(collect(months)))
    !isempty(ms) && all(in(1:12), ms) ||
        error("`months` are months of the year, 1 to 12; got $(repr(months)).")
    request = Dict{String, Any}("product_type" =>
                                    ["monthly_averaged_reanalysis"],
                                "variable" => [String(variable)],
                                "year" =>
                                    string.(sort!(unique(collect(years)))),
                                "month" => lpad.(string.(ms), 2, '0'),
                                "time" => ["00:00"],
                                "data_format" => "netcdf",
                                "download_format" => "unarchived")
    isnothing(area) || (request["area"] = collect(area))
    return request
end

# The CDS dataset ERA5 monthly means on single levels are served from.
const _ERA5_MONTHLY = "reanalysis-era5-single-levels-monthly-means"

# Search a directory `path` for the files whose names contain `key` -- a plain `occursin` test, not a
# glob or a regular expression - leaving out the provenance sidecars, whose names contain their
# file's. It knows nothing about climate data: its callers are listing a cache or data folder,
# which is why it sits beside the other filesystem helpers and is private.
function _searchdir(path, key)
    return filter(x -> occursin(key, x) && !endswith(x, ".provenance.toml"),
                  readdir(path))
end
