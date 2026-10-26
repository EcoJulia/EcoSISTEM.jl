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
    dir = joinpath(_assetroot(), string(nameof(mod)))
    isnothing(owner) || (dir = joinpath(dir, string(nameof(owner))))
    return mkpath(dir)
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
    path = mkpath(abspath(expanduser(dir)))
    # Into the environment that has EcoSISTEM rather than the active one, which may not name it and
    # would then be given it in `[extras]`, its Project.toml rewritten
    Preferences.set_preferences!(EcoSISTEM, "assetdir" => path,
                                 active_project_only = false, force = true)
    _repointextensions()
    return _assetroot()
end

function set_assetdir!(::Nothing)
    Preferences.has_preference(EcoSISTEM, "assetdir") &&
        Preferences.delete_preferences!(EcoSISTEM, "assetdir",
                                        active_project_only = false,
                                        force = true)
    _repointextensions()
    return _assetroot()
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
    if isfile(path)
        _recordpresent(path) do
            record = _assetrecord(asset, path, nothing)
            return record
        end
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
        _checkrequest(request)
        _recordpresent(request.path) do
            return _requestrecord(request, request.path, nothing)
        end
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

# The asset cache's root: the directory `set_assetdir!` recorded, else EcoSISTEM's scratch space.
function _assetroot()
    dir = Preferences.load_preference(EcoSISTEM, "assetdir", nothing)
    return isnothing(dir) ? get_scratch!(EcoSISTEM, "assets") : String(dir)
end

# Tell each loaded extension that keeps files under the asset cache that the cache has moved.
function _repointextensions()
    isnothing(Base.get_extension(EcoSISTEM, :EcoSISTEMRasterDataSourcesExt)) ||
        _repointrasterdatasources()
    return nothing
end

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

# Record a file that is already present when it is first used and has no sidecar - one fetched
# before records were written, or by hand - with everything but a fetch time, which is not known.
# A sidecar that is there already, ours or another program's, is left as it is.
function _recordpresent(record::Function, path::AbstractString)
    isfile(_sidecarpath(path)) && return nothing
    r = record()
    delete!(r, "fetched")
    _writesidecar(path, r)
    return nothing
end

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

# Check a present file against the checksum its sidecar records, erroring on a mismatch, and say
# whether there was one to check: a file with no sidecar, one this package did not write, or a
# record with no checksum has nothing to check against and passes unchecked.
function _verifyfile(path::AbstractString)
    isfile(_sidecarpath(path)) || return false
    record = TOML.parsefile(_sidecarpath(path))
    (_ourrecord(record) && haskey(record, "sha256")) || return false
    actual = _sha256(path)
    actual == record["sha256"] ||
        error("`$(basename(path))` does not match the checksum its provenance record holds: " *
              "recorded $(record["sha256"]), found $actual. The file has been truncated or " *
              "replaced since it was fetched; delete it to fetch it again.")
    return true
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

# The catalogue's facts about a file of a catalogued source, the one place they are computed for
# every record written and every registration made: the source, the layer where the file is one,
# the file's title, and the row's DOI, licence, version, release date and citation, with what its
# provider asks - the references, the notices and the documentation's section quoting them - left
# out where blank. Empty where the owner is no source.
function _cataloguerecord(owner::Type, layer, path::AbstractString)
    rec = _datasetrecord(owner)
    isnothing(rec) && return Dict{String, Any}()
    record = Dict{String, Any}("source" => rec.dataset,
                               "title" => _filetitle(rec, layer, path))
    for (key, value) in (("doi", rec.doi), ("licence", rec.licence),
        ("version", rec.version), ("release_date", rec.releasedate),
        ("citation", rec.citation), ("request_page", _requestpage(rec)))
        isempty(value) || (record[key] = value)
    end
    isempty(rec.references) ||
        (record["references"] = [Dict("doi" => r.doi, "text" => r.text)
                                 for r in rec.references])
    isempty(rec.notices) ||
        (record["notices"] = [Dict("strength" => string(n.strength),
                                   "text" => n.text)
                              for n in rec.notices])
    isnothing(layer) || (record["code"] = first(layer.aliases))
    return record
end

# The provenance record of a file fetched over https: the URL asked for and the one the server
# answered from, and the validators it sent, plus the catalogue's facts where the owner is a
# source and the URL is one of its layers' files. `response` is `nothing` for a file recorded
# when already present, which then carries the URL alone.
function _assetrecord(asset::CachedAsset, path::AbstractString, response)
    rec = _datasetrecord(asset.owner)
    record = _baserecord(path, _assetrole(asset.owner))
    record["url"] = asset.url
    if !isnothing(response)
        response.url == asset.url || (record["final_url"] = response.url)
        for (name, key) in (("etag", "etag"), ("last-modified", "last_modified"))
            i = findfirst(h -> lowercase(first(h)) == name, response.headers)
            isnothing(i) || (record[key] = String(last(response.headers[i])))
        end
    end
    isnothing(rec) && return record
    record["primary"] = true
    layer = rec.format === :Shapefile ? nothing :
            _layerbyfile(asset.owner, asset.url)
    merge!(record, _cataloguerecord(asset.owner, layer, path))
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
# request body as sent, the job the CDS ran, and the catalogue's facts where the dataset is one
# the catalogue knows (ERA5 monthly means are `ERA`'s) and the variable one of its layers.
function _requestrecord(request::CDSRequest, path::AbstractString, status)
    record = _baserecord(path, :habitat)
    record["dataset"] = request.dataset
    record["request"] = request.request
    record["primary"] = false
    isnothing(request.region) || (record["region"] = request.region)
    job = status isa AbstractDict ? get(status, "jobID", nothing) : nothing
    isnothing(job) || (record["job"] = string(job))
    if request.dataset == _ERA5_MONTHLY
        variables = get(request.request, "variable", String[])
        layer = length(variables) == 1 ?
                _layerbyrequest(ERA, string(only(variables))) : nothing
        merge!(record, _cataloguerecord(ERA, layer, path))
    end
    return record
end

# Refuse a present file whose record of ours says it was fetched for another request - a box or
# months changed under one path - rather than serve it as this one. A file with no record, or with
# another program's, passes, presence being all there is to go on.
function _checkrequest(request::CDSRequest)
    sidecar = _sidecarpath(request.path)
    isfile(sidecar) || return nothing
    record = TOML.tryparsefile(sidecar)
    (record isa AbstractDict && _ourrecord(record) &&
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
