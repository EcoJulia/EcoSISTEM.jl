# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A file downloaded once and cached outside the repository, under `EcoSISTEM.assetdir`.

"""
    CachedAsset(owner::Type, url::AbstractString)

An immutable descriptor for a file downloaded from `url` and cached under
`EcoSISTEM.assetdir(owner = owner)`. Nothing is downloaded at construction - call
[`assetpath`](@ref) to get the local path, downloading it into the cache first if it is not there.

# Fields

  - `owner`: the type the download belongs to, which names its own subdirectory of the cache so that
    one type's downloads cannot collide with another's.
  - `url`: where to fetch it from.
"""
struct CachedAsset
    owner::Type
    url::String
    CachedAsset(owner::Type, url::AbstractString) = new(owner, String(url))
end

"""
    CDSRequest(dataset::AbstractString, request::AbstractDict, path::AbstractString)

A file fetched from the Copernicus Climate Data Store on first use and kept at `path`: the CDS
`dataset` to ask (`"reanalysis-era5-single-levels-monthly-means"`), the `request` as the CDS
takes it (product type, variable, years, months, area, format), and where the answer is written.
Nothing is fetched at construction; [`assetpath`](@ref) fetches it through `CDSAPI.jl` when
`path` does not yet exist, so a file already downloaded is never asked for twice. An entry of a
[`RasterSpec`](@ref)'s `files`.

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
    assetpath(asset::CachedAsset)
    assetpath(request::CDSRequest)
    assetpath(relative_path::AbstractString)

Return a local path to a file, downloading it first if it is not there yet.

# Arguments

  - `asset`: a [`CachedAsset`](@ref), fetched from its `url` into its owner's subdirectory of
    EcoSISTEM's asset cache if the cache does not already hold it.
  - `request`: a [`CDSRequest`](@ref), fetched from the Climate Data Store into its own `path` if
    that file does not already exist. The fetch needs `CDSAPI` loaded and a CDS key in
    `~/.cdsapirc`; a file already present needs neither.
  - `relative_path`: a path resolved under [`assetdir`](@ref) directly, for a cached file whose
    name is already known and which needs no descriptor.
"""
function assetpath(asset::CachedAsset)
    path = joinpath(assetdir(owner = asset.owner), basename(asset.url))
    isfile(path) && return path
    # Downloaded to a sibling part file and renamed in, never written straight to `path`.
    # `isfile(path)` is the only cache check, so a half-written file there would be served as valid
    # forever, failing with an opaque reader error that re-running never clears. The rename closes
    # the window a process killed outright part-way through would leave. A rename within a
    # directory is atomic, and concurrent downloads get their own part file, so `force` means
    # last-writer-wins on identical content.
    temp = path * ".part-" * string(getpid())
    try
        _download(asset.url, temp)
        mv(temp, path, force = true)
    finally
        rm(temp, force = true)          # a no-op once the rename has succeeded
    end
    return path
end

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
# range on each retry, up to `attempts` transfers, the last failure rethrown. A server that ignores
# the range and sends the whole file restarts the file rather than appending a second copy.
function _download(url::AbstractString, path::AbstractString;
                   attempts::Integer = 5)
    downloader = _http11downloader()
    for attempt in 1:attempts
        have = isfile(path) ? filesize(path) : 0
        try
            open(path, "a") do io
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
                return response.status in (0, 200, 206) ||
                       error("HTTP $(response.status) for $url")
            end
            return path
        catch e
            attempt == attempts && rethrow()
            @warn "download attempt $attempt of $attempts failed; resuming" url exception=e
        end
    end
end
assetpath(relative_path::AbstractString) = joinpath(assetdir(), relative_path)

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

function assetpath(request::CDSRequest)
    isfile(request.path) && return request.path
    isnothing(Base.get_extension(EcoSISTEM, :EcoSISTEMERAExt)) &&
        error("`$(request.path)` is not there, and fetching it from the Climate Data Store " *
              "needs `using CDSAPI` (and a CDS key in `~/.cdsapirc`).")
    _fetchcds(request)
    isfile(request.path) ||
        error("the Climate Data Store fetch for `$(request.path)` finished without writing it.")
    return request.path
end

# Search a directory `path` for the files whose names contain `key` -- a plain `occursin` test, not a
# glob or a regular expression. It knows nothing about climate data: both callers, `clearcache!` and
# `gettimes`, are listing a cache folder, which is why it sits beside the other filesystem helpers
# and is private.
_searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))
