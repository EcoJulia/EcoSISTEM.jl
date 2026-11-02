# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The asset cache a package keeps its downloads in, the download that resumes, locks and records
# itself, the provenance sidecar every fetched or first-used file gets, and the recording of a
# file the user has from what is known of its source.

"""
    record_input(path::AbstractString, facts::SourceFacts; role, route = PublishedByProvider(),
                 url = nothing, title = nothing, code = nothing, sidecar = sidecarpath(path))
    record_input(path::AbstractString, doi::AbstractString; role, route = PublishedByProvider(),
                 url = nothing, title = nothing, code = nothing, sidecar = sidecarpath(path))

Record where the file at `path` came from, once, and return its [`InputRecord`](@ref): the
provenance sidecar is written beside the file as `<file>.provenance.toml`, or at `sidecar`, from
what the file itself says - its name, size and checksum, under `role` - and what `facts` say of
its source ([`sourceentries`](@ref Provenancer.sourceentries)); given a DOI, the facts are what
its metadata gives ([`lookup_doi`](@ref)). A file recorded already is left as it is and its record
returned, with nothing looked up, so a call repeated on every run costs a read; delete the sidecar
to record the file afresh. A sidecar another program wrote under that name is neither read nor
overwritten: name another place with `sidecar`.

No fetch time is recorded, since this call fetched nothing: a file fetched by
[`fetchfile`](@ref Provenancer.fetchfile) carries one. Offline, a DOI cannot be looked up, and the
error says so: write the facts by hand as [`SourceFacts`](@ref) and give those.

# Arguments

  - `path`: the file, which must be there.
  - `facts`: what is known of its source; or `doi`, the DOI whose metadata gives them.
  - `role`: which part of the system the file enters by, a word of the package using the record.
  - `route`: how the file came to be: as its provider published it
    ([`PublishedByProvider`](@ref)), the default, so that with an identifier naming the data it is
    the cited item or one of its files; made by a service to a request, carrying the identifier
    the service gave it ([`ExtractedByService`](@ref)); or altered or made by you
    ([`SuppliedByUser`](@ref)). With the facts, it decides whether the record calls the file
    `primary` ([`isprimary`](@ref Provenancer.isprimary)).
  - `url`: where the file was fetched from, where known.
  - `title`: the file's own title, the source's by default.
  - `code`: which layer of the source the file holds, where it is one.
  - `sidecar`: where the record is written, the file's `.provenance.toml` beside it by default;
    somewhere else for a file in a place that cannot be written to. [`provenance`](@ref) takes
    the same to read it back.
"""
function record_input(path::AbstractString, facts::SourceFacts; role::Symbol,
                      route::AbstractFileRoute = PublishedByProvider(),
                      url::Union{Nothing, AbstractString} = nothing,
                      title::Union{Nothing, AbstractString} = nothing,
                      code = nothing,
                      sidecar::AbstractString = sidecarpath(path))
    recorded = _priorrecord(path, sidecar)
    isnothing(recorded) || return recorded
    recordpresent(path, sidecar = sidecar) do
        record = baserecord(path, role)
        merge!(record,
               sourceentries(facts, route,
                             title = something(title, facts.title),
                             code = code))
        isnothing(url) || (record["url"] = String(url))
        return record
    end
    return recordat(path, sidecar = sidecar)
end

function record_input(path::AbstractString, doi::AbstractString; role::Symbol,
                      sidecar::AbstractString = sidecarpath(path), kw...)
    recorded = _priorrecord(path, sidecar)
    isnothing(recorded) || return recorded
    return record_input(path, lookup_doi(doi); role = role, sidecar = sidecar,
                        kw...)
end

"""
    assetdir(mod::Module; owner::Union{Type, Nothing} = nothing, cache::Module = mod)

Return the path of `mod`'s subdirectory of an asset cache, for storing downloaded data outside the
repository, making the directory if need be: the cache's root - the directory
[`set_assetdir!`](@ref Provenancer.set_assetdir!) recorded for the package, read on every call,
else the package's Scratch.jl space, created on first use and reclaimed by `Pkg.gc()` when the
package is removed - then `mod`'s name, and the owner's where one is given.

# Arguments

  - `mod`: the module whose subdirectory it is, named by it.
  - `owner`: a type to nest a further subdirectory under, so each owner's downloads are separate
    and one type's cannot collide with another's.
  - `cache`: the package whose cache root the directory sits in, `mod`'s own by default. A package
    keeping another's downloads under its own root, so that one lifecycle covers every download,
    names itself here.
"""
function assetdir(mod::Module; owner::Union{Type, Nothing} = nothing,
                  cache::Module = mod)
    dir = joinpath(_assetroot(cache), string(nameof(mod)))
    isnothing(owner) || (dir = joinpath(dir, string(nameof(owner))))
    return mkpath(dir)
end

"""
    set_assetdir!(mod::Module, dir::AbstractString)
    set_assetdir!(mod::Module, nothing)

Keep the asset cache of the package `mod` belongs to in `dir` from now on, rather than in its
Scratch.jl space, which lives in the first Julia depot and may be on a disk too small for the
data; `nothing` returns the cache to the scratch space. Return the directory the cache is now in,
made if need be.

The choice is a preference in the `LocalPreferences.toml` of the first environment on the load
path that has the package, so it holds for every later session taking the package from that
environment, and [`assetdir`](@ref Provenancer.assetdir) reads it on every call. Files already
downloaded stay where they are, and are fetched again unless moved across. `Pkg.gc()` never
reclaims a directory chosen here.

# Arguments

  - `mod`: the package whose cache it is, or any module within it.
  - `dir`: the directory to keep the cache in, `~` and a relative path expanded.
"""
function set_assetdir!(mod::Module, dir::AbstractString)
    root = Base.moduleroot(mod)
    path = mkpath(abspath(expanduser(dir)))
    # Into the environment that has the package rather than the active one, which may not name it
    # and would then be given it in `[extras]`, its Project.toml rewritten.
    Preferences.set_preferences!(root, "assetdir" => path,
                                 active_project_only = false, force = true)
    return _assetroot(root)
end

function set_assetdir!(mod::Module, ::Nothing)
    root = Base.moduleroot(mod)
    Preferences.has_preference(root, "assetdir") &&
        Preferences.delete_preferences!(root, "assetdir",
                                        active_project_only = false,
                                        force = true)
    return _assetroot(root)
end

"""
    fetchlocked(fetch::Function, path::AbstractString)

Run `fetch` under the lock for the file at `path`, unless the file has appeared by the time the
lock is held, which is what a process that waited on another's download finds; return `path`. The
lock is a pid file beside the target, released on return however `fetch` ends; a lock older than
four hours whose process is gone is taken over.

# Arguments

  - `fetch`: what fetches the file, called with no arguments.
  - `path`: where the file will be.
"""
function fetchlocked(fetch::Function, path::AbstractString)
    lock = mkpidlock(path * ".lock", stale_age = _STALE_LOCK_SECONDS)
    try
        isfile(path) || fetch()
    finally
        close(lock)
    end
    return path
end

"""
    download!(url::AbstractString, path::AbstractString; attempts = 5)

Fetch `url` to `path` over HTTP/1.1, appending to whatever `path` already holds through a byte
range on each retry, and return the last response, for its final URL and headers. The last of
`attempts` failures is rethrown, the part file left for a later call to resume from. A server that
ignores the range and sends the whole file restarts the file rather than appending a second copy.

# Arguments

  - `url`: where to fetch from.
  - `path`: the file to write, resumed where it exists.
  - `attempts`: how many transfers to try before giving up.
"""
function download!(url::AbstractString, path::AbstractString;
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

"""
    fetchfile(record::Function, url::AbstractString, path::AbstractString; verify = false)
    fetchfile(file::URLFile; role::Symbol, verify = false)

Return `path`, the file fetched from `url` into it first where it is not there, with its
provenance record written beside it when it is: `record(response)` gives the record, as a
dictionary of the sidecar's keys, from the last response of the download, or from `nothing` for a
file found present with no record, which is recorded without a fetch time. A file present with a
record is left as it is, and checked against the record's checksum where `verify` is `true`.

The download goes to a sibling part file and is renamed in, never written straight to `path`: a
present file is taken as whole, so a half-written one there would be served as valid forever. One
part file per destination, under a lock, so a second process wanting the same file waits and then
finds it present, and the part file of an interrupted run is what the next run resumes from.

A [`URLFile`](@ref) is recorded with what its download says ([`downloadrecord`](@ref
Provenancer.downloadrecord)) and what its facts say of the source ([`sourceentries`](@ref
Provenancer.sourceentries)), under `role`.

# Arguments

  - `record`: what builds the record, called with the download's last response, or `nothing`.
  - `url`: where to fetch from.
  - `path`: where the file lives.
  - `file`: the file to fetch, as a `URLFile`.
  - `role`: which part of the system the file enters by, a word of the package using the record.
  - `verify`: whether to check a file already present against its record's checksum.
"""
function fetchfile(record::Function, url::AbstractString, path::AbstractString;
                   verify::Bool = false)
    if isfile(path)
        recordpresent(path) do
            return record(nothing)
        end
        verify && verifyfile(path)
        return path
    end
    mkpath(dirname(abspath(path)))
    part = path * ".part"
    fetchlocked(path) do
        response = download!(url, part)
        mv(part, path, force = true)
        return writesidecar(path, record(response))
    end
    return path
end

function fetchfile(file::URLFile; role::Symbol, verify::Bool = false)
    return fetchfile(file.url, file.path, verify = verify) do response
        record = downloadrecord(file.path, file.url, response, role)
        isnothing(file.facts) ||
            merge!(record, sourceentries(file.facts, route(file)))
        return record
    end
end

"""
    recordpresent(record::Function, path::AbstractString; sidecar = sidecarpath(path))

Record a file that is already present when it is first used and has no sidecar - one fetched
before records were written, or by another program, or given by the user - writing the record
`record()` returns beside it with everything but a fetch time, which is not known. A sidecar that
is there already, ours or another program's, is left as it is.

# Arguments

  - `record`: what builds the record, called with no arguments; a dictionary of the sidecar's keys.
  - `path`: the file.
  - `sidecar`: where its record goes, the file's `.provenance.toml` beside it by default.
"""
function recordpresent(record::Function, path::AbstractString;
                       sidecar::AbstractString = sidecarpath(path))
    isfile(sidecar) && return nothing
    r = record()
    delete!(r, "fetched")
    writesidecar(path, r, sidecar = sidecar)
    return nothing
end

"""
    writesidecar(path::AbstractString, record::AbstractDict; sidecar = sidecarpath(path))

Write `record` as the provenance sidecar of the file at `path`, in TOML with its keys sorted so two
records of one file diff cleanly, and return the sidecar's path. Written to a sibling and renamed
in, never into the file already there: a sidecar a data pipeline run put in place is a hard link
to the run's store, and writing into it would change the store's copy. The sibling is each
writer's own, since several processes may record one file at once.

# Arguments

  - `path`: the file, not its sidecar.
  - `record`: the sidecar's keys and values.
  - `sidecar`: where to write it, the file's `.provenance.toml` beside it by default; its
    directory is made if need be.
"""
function writesidecar(path::AbstractString, record::AbstractDict;
                      sidecar::AbstractString = sidecarpath(path))
    part, io = mktemp(mkpath(dirname(abspath(sidecar))))
    try
        TOML.print(io, record, sorted = true)
    finally
        close(io)
    end
    chmod(part, 0o644)
    mv(part, sidecar, force = true)
    return sidecar
end

"""
    verifyfile(path::AbstractString; sidecar = sidecarpath(path))

Check a present file against the checksum its sidecar records, erroring on a mismatch, and return
whether there was one to check: a file with no sidecar, one this machinery did not write, or a
record with no checksum has nothing to check against and passes unchecked.

# Arguments

  - `path`: the file, not its sidecar.
  - `sidecar`: where its record is, the file's `.provenance.toml` beside it by default.
"""
function verifyfile(path::AbstractString;
                    sidecar::AbstractString = sidecarpath(path))
    isfile(sidecar) || return false
    record = TOML.parsefile(sidecar)
    (_isourrecord(record) && haskey(record, "sha256")) || return false
    actual = _filesha256(path)
    actual == record["sha256"] ||
        error("`$(basename(path))` does not match the checksum its provenance record holds: " *
              "recorded $(record["sha256"]), found $actual. The file has been truncated or " *
              "replaced since it was fetched; delete it to fetch it again.")
    return true
end

"""
    baserecord(path::AbstractString, role::Symbol)

Return the fields every fetched file's record carries, as a dictionary of the sidecar's keys: who
wrote it (`writer`), its `role`, the file's name (`file`; never its directory, since the record
travels), when it was fetched (`fetched`, now), its size (`bytes`) and checksum (`sha256`). The
name is relative by construction and no field ever holds a home directory, user name, host name
or key.

# Arguments

  - `path`: the file.
  - `role`: which part of the system the file enters by, a word of the package using the record.
"""
function baserecord(path::AbstractString, role::Symbol)
    return Dict{String, Any}("writer" => _writer(),
                             "role" => string(role),
                             "file" => basename(path),
                             "fetched" => _utcnow(),
                             "bytes" => filesize(path),
                             "sha256" => _filesha256(path))
end

"""
    downloadrecord(path::AbstractString, url::AbstractString, response, role::Symbol)

Return the provenance record of a file fetched over https, as a dictionary of the sidecar's keys:
the [`baserecord`](@ref Provenancer.baserecord), the URL asked for (`url`) and the one the server
answered from where it differs (`final_url`), the validators it sent (`etag`, `last_modified`),
and when the provider last changed the file (`modified`: its own time where the server keeps it as
metadata, else its `Last-Modified`). `response` is what
[`download!`](@ref Provenancer.download!) returned, or `nothing` for a file recorded when already
present, which then carries the URL alone.

# Arguments

  - `path`: the file.
  - `url`: where it was fetched from.
  - `response`: the download's last response, or `nothing`.
  - `role`: which part of the system the file enters by.
"""
function downloadrecord(path::AbstractString, url::AbstractString, response,
                        role::Symbol)
    record = baserecord(path, role)
    record["url"] = url
    isnothing(response) && return record
    response.url == url || (record["final_url"] = response.url)
    for (name, key) in (("etag", "etag"), ("last-modified", "last_modified"))
        value = _headervalue(response.headers, name)
        isnothing(value) || (record[key] = value)
    end
    modified = _providertime(response.headers)
    isnothing(modified) || (record["modified"] = modified)
    return record
end

# What a call recording the file at `path` finds before it writes: the file must be there; a record
# of ours at `sidecar` is returned as it is, no sidecar at all is `nothing`, and a sidecar there
# that is not a record of ours - another program's, or one that does not read - is refused, since
# it is neither read nor overwritten.
function _priorrecord(path::AbstractString, sidecar::AbstractString)
    isfile(path) ||
        throw(ArgumentError("`$path` is not there, so it cannot be recorded: a record carries " *
                            "the file's size and checksum."))
    isfile(sidecar) || return nothing
    record = recordat(path, sidecar = sidecar)
    isnothing(record) &&
        error("`$sidecar` is there but is not a record this machinery wrote and can read, so " *
              "`$(basename(path))` is neither read from it nor recorded over it. Move it, or " *
              "give `sidecar` another path to write the record to.")
    return record
end

# The root of the asset cache of the package `mod` belongs to: the directory `set_assetdir!`
# recorded for it as a preference, read on every call, else the package's Scratch.jl space.
function _assetroot(mod::Module)
    root = Base.moduleroot(mod)
    dir = Preferences.load_preference(root, "assetdir", nothing)
    return isnothing(dir) ? get_scratch!(root, "assets") : String(dir)
end

# How old a download lock may be before a process that stopped holding it is presumed dead.
const _STALE_LOCK_SECONDS = 4 * 3600.0

# A `Downloads.Downloader` held to HTTP/1.1: some networks cut a long HTTP/2 stream off part way
# through with an `INTERNAL_ERROR`, and the files worth caching are hundreds of megabytes each.
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

# The SHA-256 digest of a file, as hex.
_filesha256(path::AbstractString) = bytes2hex(open(SHA.sha256, path))

# Now, in UTC, written as ISO 8601 with its `Z`.
_utcnow() = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SS") * "Z"

# A response header's value by its lower-case name, matched in any case, or `nothing` where it was
# not sent.
function _headervalue(headers, name::AbstractString)
    i = findfirst(h -> lowercase(first(h)) == name, headers)
    return isnothing(i) ? nothing : String(last(headers[i]))
end

# When the provider last changed the file a response describes, as ISO 8601 text in UTC with its
# `Z`: the file's own time where the server keeps it as metadata (an object store sends it as
# `x-amz-meta-mtime`, its `Last-Modified` being when it was uploaded there), else the
# `Last-Modified` it sent; `nothing` where it sent neither in a form that reads.
function _providertime(headers)
    mtime = _headervalue(headers, "x-amz-meta-mtime")
    seconds = isnothing(mtime) ? nothing : tryparse(Float64, mtime)
    isnothing(seconds) ||
        return Dates.format(Dates.unix2datetime(floor(seconds)),
                            "yyyy-mm-ddTHH:MM:SS") * "Z"
    modified = _headervalue(headers, "last-modified")
    isnothing(modified) && return nothing
    t = tryparse(Dates.DateTime, replace(modified, r"\s+GMT$" => ""),
                 Dates.dateformat"e, d u y H:M:S")
    return isnothing(t) ? nothing :
           Dates.format(t, "yyyy-mm-ddTHH:MM:SS") * "Z"
end

# The provider's time for the file at `url`, asked by a `HEAD` and read as `_providertime` reads a
# response, where the answer can be tied to the bytes held at `held`: the size the server states
# must be that file's, since a record may be written long after the file was fetched and the
# server's file changed since. `nothing` where the server cannot be asked, states no time or no
# size, or describes another file.
function _headprovidertime(url::AbstractString, held::AbstractString)
    response = try
        Downloads.request(url, method = "HEAD", timeout = 30,
                          downloader = _http11downloader(), throw = false)
    catch
        return nothing
    end
    response.status == 200 || return nothing
    size = _headervalue(response.headers, "content-length")
    (isnothing(size) || tryparse(Int, size) != filesize(held)) && return nothing
    return _providertime(response.headers)
end
