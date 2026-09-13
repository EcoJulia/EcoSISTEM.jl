# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A file downloaded once and cached outside the repository, under `EcoSISTEM.assetdir`, or fetched
# once to a path of the caller's; and the provenance record every such fetch writes beside it.

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

"""
    CDSRequest(dataset::AbstractString, request::AbstractDict, path::AbstractString)

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
"""
struct CDSRequest
    dataset::String
    request::Dict{String, Any}
    path::String
    function CDSRequest(dataset::AbstractString, request::AbstractDict,
                        path::AbstractString)
        return new(String(dataset),
                   Dict{String, Any}(string(k) => v for (k, v) in request),
                   String(path))
    end
end

"""
    EcoSISTEM.assetdir(mod::Module = EcoSISTEM; owner::Union{Type, Nothing} = nothing)

Path to a subdirectory of EcoSISTEM's Scratch.jl space, for storing downloaded data outside the
repository. Creating the directory is part of answering.

One EcoSISTEM-owned space with a subdirectory per package, rather than a space each, keeps the whole
cache under one lifecycle: created on first use, and reclaimed by `Pkg.gc()` when EcoSISTEM is
removed. `RasterDataSources` is put under it the same way, by the extension's `__init__` setting
`RASTERDATASOURCES_PATH`.

# Arguments

  - `mod`: whose subdirectory to return, defaulting to EcoSISTEM's own.
  - `owner`: a type to nest a further subdirectory under, so each owner's downloads are separate -
    see [`CachedAsset`](@ref).
"""
function assetdir(mod::Module = EcoSISTEM;
                  owner::Union{Type, Nothing} = nothing)
    dir = joinpath(get_scratch!(EcoSISTEM, "assets"), string(nameof(mod)))
    isnothing(owner) || (dir = joinpath(dir, string(nameof(owner))))
    return mkpath(dir)
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
    `~/.cdsapirc`; a file already present needs neither.
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
    if isfile(path)
        verify && _verifyfile(path)
        return path
    end
    mkpath(dirname(abspath(path)))
    # Downloaded to a sibling part file and renamed in, never written straight to `path`.
    # `isfile(path)` is the only cache check, so a half-written file there would be served as valid
    # forever, failing with an opaque reader error that re-running never clears; the rename, atomic
    # within a directory, closes the window a process killed outright would leave. One part file
    # per destination, under a lock: a second process wanting the same file waits, then finds it
    # present, and the part file of an interrupted run is what the next run resumes from.
    part = path * ".part"
    _fetchlocked(path) do
        response = _download(asset.url, part)
        mv(part, path, force = true)
        return _writesidecar(path, _assetrecord(asset, path, response))
    end
    return path
end

function assetpath(request::CDSRequest; verify::Bool = false)
    if isfile(request.path)
        verify && _verifyfile(request.path)
        return request.path
    end
    isnothing(Base.get_extension(EcoSISTEM, :EcoSISTEMERAExt)) &&
        error("`$(request.path)` is not there, and fetching it from the Climate Data Store " *
              "needs `using CDSAPI` (and a CDS key in `~/.cdsapirc`).")
    mkpath(dirname(abspath(request.path)))
    _fetchlocked(request.path) do
        status = _fetchcds(request)
        isfile(request.path) ||
            error("the Climate Data Store fetch for `$(request.path)` finished without writing it.")
        return _writesidecar(request.path,
                             _requestrecord(request, request.path, status))
    end
    return request.path
end

assetpath(relative_path::AbstractString) = joinpath(assetdir(), relative_path)

# Run `fetch` under the lock for `path`, unless the file has appeared by the time the lock is held,
# which is what a process that waited on another's download finds. The lock is a pid file beside
# the target, released on return however `fetch` ends; a lock older than `_STALE_LOCK_SECONDS`
# whose process is gone is taken over.
function _fetchlocked(fetch::Function, path::AbstractString)
    lock = mkpidlock(path * ".lock", stale_age = _STALE_LOCK_SECONDS)
    try
        isfile(path) || fetch()
    finally
        close(lock)
    end
    return path
end

# How old a download lock may be before a process that stopped holding it is presumed dead.
const _STALE_LOCK_SECONDS = 4 * 3600.0

# A downloader held to HTTP/1.1. Some networks cut a long HTTP/2 stream off part way through with
# an `INTERNAL_ERROR`, and the reanalysis archives are hundreds of megabytes a file.
function _http11downloader()
    downloader = Downloads.Downloader()
    downloader.easy_hook = (easy, _) -> Downloads.Curl.setopt(easy,
                                                              Downloads.Curl.CURLOPT_HTTP_VERSION,
                                                              Downloads.Curl.CURL_HTTP_VERSION_1_1)
    return downloader
end

# The request headers that resume a transfer from the `have` bytes already held: a byte range
# from there, or none at all from the start.
function _resumeheaders(have::Integer)
    return have > 0 ? ["Range" => "bytes=$(have)-"] : Pair{String, String}[]
end

# Fetch `url` to `path` over HTTP/1.1, appending to whatever `path` already holds through a byte
# range on each retry, up to `attempts` transfers, the last failure rethrown; the part file is left
# for a later call to resume from. A server that ignores the range and sends the whole file
# restarts the file rather than appending a second copy. Returns the last response, for its final
# URL and headers.
function _download(url::AbstractString, path::AbstractString;
                   attempts::Integer = 5)
    downloader = _http11downloader()
    for attempt in 1:attempts
        have = isfile(path) ? filesize(path) : 0
        try
            return open(path, "a") do io
                response = Downloads.request(url, output = io,
                                             headers = _resumeheaders(have),
                                             downloader = downloader,
                                             throw = true)
                # A file URL reports no status at all; a server answers 200 for the whole file and
                # 206 for the range asked.
                if have > 0 && response.status == 200
                    truncate(io, 0)
                    error("the server sent the whole file in answer to a byte range")
                end
                response.status in (0, 200, 206) ||
                    error("HTTP $(response.status) for $url")
                return response
            end
        catch e
            if attempt == attempts
                # Nothing received at all leaves nothing worth resuming from, and a stray empty
                # part file would read as a download in progress.
                isfile(path) && filesize(path) == 0 && rm(path)
                rethrow()
            end
            @warn "download attempt $attempt of $attempts failed; resuming" url exception=e
        end
    end
end

# The sidecar a fetched file's provenance record is written to: the file's own name with
# `.provenance.toml` after it, beside it.
_sidecarpath(path::AbstractString) = path * ".provenance.toml"

# Write `record` as the provenance sidecar of `path`, in TOML, keys sorted so two records of one
# file diff cleanly.
function _writesidecar(path::AbstractString, record::AbstractDict)
    open(_sidecarpath(path), "w") do io
        return TOML.print(io, record, sorted = true)
    end
    return _sidecarpath(path)
end

# The SHA-256 digest of a file, as hex.
_sha256(path::AbstractString) = bytes2hex(open(SHA.sha256, path))

# Check a present file against the checksum its sidecar records, erroring on a mismatch; a file
# with no sidecar, or a record with no checksum, has nothing to check against and passes.
function _verifyfile(path::AbstractString)
    isfile(_sidecarpath(path)) || return nothing
    record = TOML.parsefile(_sidecarpath(path))
    haskey(record, "sha256") || return nothing
    actual = _sha256(path)
    actual == record["sha256"] ||
        error("`$(basename(path))` does not match the checksum its provenance record holds: " *
              "recorded $(record["sha256"]), found $actual. The file has been truncated or " *
              "replaced since it was fetched; delete it to fetch it again.")
    return nothing
end

# The fields every fetched file's record carries: who wrote it, when, the file's name (never its
# directory, since the record travels), its size and checksum. The name is relative by
# construction and no field ever holds a home directory, user name, host name or key.
function _baserecord(path::AbstractString, role::Symbol)
    return Dict{String, Any}("writer" => "EcoSISTEM $(pkgversion(EcoSISTEM))",
                             "role" => string(role),
                             "file" => basename(path),
                             "fetched" => _utcnow(),
                             "bytes" => filesize(path),
                             "sha256" => _sha256(path))
end

# Now, in UTC, written as ISO 8601 with its `Z`.
_utcnow() = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SS") * "Z"

# The catalogue's facts about a fetched file, where the file is a layer of a catalogued source:
# the source, the layer, and the row's DOI, licence, version and citation. Empty where the owner
# is no source, or the file is no layer of it.
function _cataloguerecord(owner::Type, layer)
    rec = _datasetrecord(owner)
    isnothing(rec) && return Dict{String, Any}()
    record = Dict{String, Any}("source" => rec.dataset, "doi" => rec.doi,
                               "licence" => rec.licence,
                               "version" => rec.version,
                               "citation" => rec.citation)
    isnothing(layer) || (record["code"] = first(layer.aliases))
    return record
end

# The provenance record of a file fetched over https: the URL asked for and the one the server
# answered from, and the validators it sent, plus the catalogue's facts where the owner is a
# source and the URL is one of its layers' files.
function _assetrecord(asset::CachedAsset, path::AbstractString, response)
    rec = _datasetrecord(asset.owner)
    record = _baserecord(path, isnothing(rec) ? :region : :habitat)
    record["url"] = asset.url
    response.url == asset.url || (record["final_url"] = response.url)
    for (name, key) in (("etag", "etag"), ("last-modified", "last_modified"))
        i = findfirst(h -> lowercase(first(h)) == name, response.headers)
        isnothing(i) || (record[key] = String(last(response.headers[i])))
    end
    isnothing(rec) ||
        merge!(record,
               _cataloguerecord(asset.owner,
                                _layerbyfile(asset.owner, asset.url)))
    return record
end

# The provenance record of a file fetched from the Climate Data Store: the dataset asked, the
# request body as sent, the job the CDS ran, and the catalogue's facts where the dataset is one
# the catalogue knows (ERA5 monthly means are `ERA`'s) and the variable one of its layers.
function _requestrecord(request::CDSRequest, path::AbstractString, status)
    record = _baserecord(path, :habitat)
    record["dataset"] = request.dataset
    record["request"] = request.request
    job = status isa AbstractDict ? get(status, "jobID", nothing) : nothing
    isnothing(job) || (record["job"] = string(job))
    if request.dataset == _ERA5_MONTHLY
        variables = get(request.request, "variable", String[])
        layer = length(variables) == 1 ?
                _layerbyrequest(ERA, string(only(variables))) : nothing
        merge!(record, _cataloguerecord(ERA, layer))
    end
    return record
end

# The body of a Climate Data Store request for ERA5 monthly means on single levels: one
# `variable` by its CDS name, the `years` and `months` as two-digit strings, an `area` as
# `[north, west, south, east]` in degrees where given, netCDF unarchived. The shape the CDS takes
# and `CDSAPI.retrieve` sends.
function _era5request(variable::AbstractString, years; months = 1:12,
                      area = nothing)
    request = Dict{String, Any}("product_type" =>
                                    ["monthly_averaged_reanalysis"],
                                "variable" => [String(variable)],
                                "year" => string.(collect(years)),
                                "month" => lpad.(string.(collect(months)), 2,
                                      '0'),
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
