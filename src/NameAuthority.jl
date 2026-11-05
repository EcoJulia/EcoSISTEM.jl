# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Name authorities: what a species name means, by a checklist's own identifiers. The abstract type
# and its interface (`resolve_name`, `lineage`, `provenance`), and the one authority shipped, the
# World Flora Online Plant List, read from its Darwin Core release into a table cached beside the
# asset cache.

using CSV: CSV
using SHA: sha256

"""
    AbstractNameAuthority

A checklist of names that says which species a name means, by the checklist's own identifiers: a
recognised synonym is followed to its accepted species, so tables naming species differently
still meet on one identifier. A [`SpeciesPool`](@ref) built with one standardises every name it is
given through it. A subtype answers [`EcoSISTEM.resolve_name`](@ref) and
[`EcoSISTEM.lineage`](@ref), and its release is an input recorded in the provenance
([`provenance`](@ref)). [`WFOPlantList`](@ref) is the one for plants.
"""
abstract type AbstractNameAuthority end

"""
    WFOPlantList(; release, allow_unchecked = false, prefer_accepted = false, path = nothing,
                 track = nothing)

The World Flora Online Plant List as a [name authority](@ref EcoSISTEM.AbstractNameAuthority): the
consensus checklist of vascular plants and bryophytes, issued as a release every June and
December, each published on Zenodo under its own DOI with CC0. A release is read from its Darwin
Core file, fetched into the asset cache on first use with the catalogue's provenance record, and
its names are matched exactly, as written in its `scientificName` column, without authors.

Each of the file's rows is a name with a WFO identifier and a status: **Accepted** (a species of
its own), **Synonym** (another name for an accepted species, which it names) or **Unchecked** (a
name not yet placed). A name resolves to a species by these rules:

 1. Its Accepted and Synonym rows are used where it has any. Otherwise its Unchecked rows are used
    with `allow_unchecked`, each one standing for itself; a name with none of these is not
    recognised.
 2. Where those rows all reach one accepted species, the name is that species, under its accepted
    name and that name's identifier: a synonym is translated to the name it is a synonym of.
 3. Otherwise the name is a homonym - several species written alike - and is kept under its own
    name with its rows' identifiers joined by `+`, so that every table naming it still finds one
    species: `"wfo-0000185290+wfo-0000202067"`. With `prefer_accepted`, only its Accepted rows'
    identifiers are joined, where it has any, so a homonym that is accepted once takes that
    identifier alone.

An authorship given with a name picks among a homonym's rows first, by the file's
`scientificNameAuthorship` (spaces trimmed and runs of them collapsed).

# Arguments

  - `release`: the release, as WFO names it - `"2026-06"`. Required, so that a script always says
    which release it was resolved against; the releases this package catalogues are the ones
    accepted, and the error naming none lists them.
  - `allow_unchecked`: `true` keeps names WFO has not yet checked, each under its own identifier,
    where a name has no Accepted or Synonym row.
  - `prefer_accepted`: `true` resolves a homonym with one Accepted row to that row alone.
  - `path`: where the release's file lives once fetched, or `nothing` for its directory in the asset
    cache; how a project keeps its own copy.
  - `track`: a [`EcoSISTEM.Track`](@ref) naming the data pipeline run the release's file goes
    through, or `nothing`, the default, for none.
"""
struct WFOPlantList <: AbstractNameAuthority
    release::String
    allow_unchecked::Bool
    prefer_accepted::Bool
    path::Union{Nothing, String}
    track::Union{Nothing, Track}

    function WFOPlantList(; release = _norelease(),
                          allow_unchecked::Bool = false,
                          prefer_accepted::Bool = false,
                          path::Union{Nothing, AbstractString} = nothing,
                          track::Union{Nothing, Track} = nothing)
        _wforelease(release)
        _refusefetchedrelease(track, nothing)
        return new(String(release), allow_unchecked, prefer_accepted,
                   isnothing(path) ? nothing : String(path), track)
    end
end

# Its file is the release as WFO deposited it.
@traitimpl IsFetchable{WFOPlantList}

Provenancer.route(::WFOPlantList) = VersionOfRecord()

function Base.show(io::IO, a::WFOPlantList)
    print(io, "WFOPlantList(release = ", repr(a.release))
    a.allow_unchecked && print(io, ", allow_unchecked = true")
    a.prefer_accepted && print(io, ", prefer_accepted = true")
    isnothing(a.path) || print(io, ", path = ", repr(a.path))
    isnothing(a.track) || print(io, ", track = ", a.track)
    return print(io, ")")
end

# A release's names as the reader keeps them in memory: one entry per row of the file, with the
# row each name's species is (a synonym's accepted row, any other row itself) and its parent in
# the classification, and the rows of each name and of each identifier. Read once per session.
struct _WFONames
    id::Vector{String}
    name::Vector{String}
    authorship::Vector{String}
    rank::Vector{Symbol}
    status::Vector{Symbol}
    accepted::Vector{Int32}
    parent::Vector{Int32}
    byname::Dict{String, Vector{Int32}}
    byid::Dict{String, Int32}
end

# == Functions ==================================================================================

"""
    provenance(authority::WFOPlantList)

Return the provenance record of the release's file, as a one-element vector, `nothing` in it where
the file has not been fetched yet - as [`provenance`](@ref) answers for a spec.

# Arguments

  - `authority`: the release.
"""
function provenance(authority::WFOPlantList)
    return Union{Nothing, ProvenanceRecord}[provenance(_localpath(_wfoasset(authority)))]
end

"""
    resolve_name(authority::AbstractNameAuthority, name::AbstractString; authorship = nothing)

Return which species `name` means to `authority`, as a named tuple: the species' `name` (an accepted
name, or a homonym's own), its `id` (the authority's identifier, or several joined by `+` for a
homonym), and how it `match`ed - `:accepted`, `:synonym` (translated to its accepted name),
`:homonym`, or `:unchecked` (a name kept under its own identifier). A name not recognised has
`nothing` for both, its `match` saying why: `:not_found`, or `:unchecked_only` for a name the
authority holds but has not checked. [`WFOPlantList`](@ref) states its rules.

# Arguments

  - `authority`: the checklist to ask.
  - `name`: the name, as the authority writes it, without authors.
  - `authorship`: the name's authors, as the authority writes them, which pick among a homonym's
    rows; `nothing`, the default, for none.
"""
function resolve_name end

function resolve_name(authority::WFOPlantList, name::AbstractString;
                      authorship::Union{Nothing, AbstractString} = nothing)
    return _resolvewfo(_wfonames(authority), String(name), authorship,
                       authority.allow_unchecked, authority.prefer_accepted)
end

"""
    lineage(authority::AbstractNameAuthority, id::AbstractString)
    lineage(authority::AbstractNameAuthority, id::AbstractString, rank::AbstractString)

Return the classification of the species `id` names, from the species itself up to its kingdom: a
vector of `(rank, name, id)` named tuples, one per level in order upwards - species, genus,
subfamily, family, order and so on, as the authority classifies it. A rank may appear more than
once (`"unranked"`), and a level the authority does not hold is absent. For a synonym's identifier
the path starts from its accepted species; for several identifiers joined by `+`, it is the levels
they share, from the lowest shared one up. Given a `rank`, return the name at that rank alone - the
lowest, where the rank appears more than once - or `missing` where the path has none, as for
species written alike in different genera:

```julia
add_trait!(pool, :genus) do row
    return lineage(wfo, coalesce(row.id, row.copyofid), "genus")
end
```

# Arguments

  - `authority`: the checklist to ask.
  - `id`: an identifier [`EcoSISTEM.resolve_name`](@ref) gave, or one of the authority's own.
  - `rank`: a rank as the authority writes it, such as `"genus"` or `"family"`.
"""
function lineage end

function lineage(authority::WFOPlantList, id::AbstractString)
    t = _wfonames(authority)
    paths = [_wfopath(t, _wforow(t, part, authority))
             for part in split(id, '+')]
    shared = _sharedtop(paths)
    return [(rank = String(t.rank[i]), name = t.name[i], id = t.id[i])
            for i in shared]
end

function lineage(authority::AbstractNameAuthority, id::AbstractString,
                 rank::AbstractString)
    path = lineage(authority, id)
    i = findfirst(l -> l.rank == rank, path)
    return isnothing(i) ? missing : path[i].name
end

# An authority as a report names it.
_authorityphrase(a::WFOPlantList) = "the WFO Plant List $(a.release)"

# The error a `WFOPlantList` given no `release` raises, naming the releases there are.
function _norelease()
    return error("`WFOPlantList` needs the `release` to resolve names against, so that a script " *
                 "always says which; this package catalogues $(_wforeleasephrase()).")
end

# The catalogue's row for a release, refusing one it does not catalogue.
function _wforelease(release::AbstractString)
    rec = datasetinfo(WFOPlantList)
    rec.version == release ||
        error("the WFO Plant List release `$release` is not catalogued; this package catalogues " *
              "$(_wforeleasephrase()).")
    return rec
end

# The releases the catalogue holds, as a phrase for a message.
_wforeleasephrase() = "release `" * datasetinfo(WFOPlantList).version * "`"

# The release's Darwin Core file, as WFO deposits it in the release's Zenodo record, which the
# catalogued DOI names, and where it is kept.
function _wfoasset(a::WFOPlantList)
    rec = _wforelease(a.release)
    record = last(split(rec.doi, "zenodo."))
    url = "https://zenodo.org/records/$record/files/" * _WFOFILE
    path = something(a.path,
                     joinpath(assetdir(owner = WFOPlantList), a.release,
                              _WFOFILE))
    return CachedAsset(WFOPlantList, url, path = path)
end

# The Darwin Core file of a release, holding one table, `classification.csv`, of every name.
const _WFOFILE = "_DwC_backbone_R.zip"

# The columns of `classification.csv` the reader keeps.
const _WFOCOLUMNS = [:taxonID, :scientificName, :scientificNameAuthorship,
    :taxonRank, :taxonomicStatus, :acceptedNameUsageID,
    :parentNameUsageID]

# The reader's own text, which keys its cache, so a change to how a release is read invalidates the
# tables read before.
const _WFOCODEHASH = bytes2hex(sha256(read(@__FILE__)))

# The names read so far this session, by the cache file they came from.
const _WFOLOADED = Dict{String, _WFONames}()

const _WFOLOCK = ReentrantLock()

# A release's names: from memory while its file is the one read, else from its cached table, else
# read from its file - fetched, or provided by the release's pipeline run, first - and cached.
function _wfonames(a::WFOPlantList)
    zip = _localpath(_wfoasset(a))
    return lock(_WFOLOCK) do
        if isfile(zip)
            names = get(_WFOLOADED, _wfocachepath(zip), nothing)
            isnothing(names) || return names
        end
        _providespec!(a)
        assetpath(_wfoasset(a))
        cache = _wfocachepath(zip)
        return get!(_WFOLOADED, cache) do
            isfile(cache) || _writewfocache(zip, cache)
            return _wfotable(cache)
        end
    end
end

# How a pool resolves many names against one authority: through `resolve_name`, or for a release
# read into memory, straight against its names.
function _resolver(authority::AbstractNameAuthority)
    return (name, authorship) -> resolve_name(authority, name,
                                              authorship = authorship)
end

function _resolver(a::WFOPlantList)
    t = _wfonames(a)
    return (name, authorship) -> _resolvewfo(t, String(name), authorship,
                                             a.allow_unchecked,
                                             a.prefer_accepted)
end

# Where a release file's reduced table is cached: named for the file's path, then for its size and
# time and the reader's text, each digested with SHA-256 (never `Base.hash`, whose values change
# between Julia versions), so a replaced file or a changed reader is read afresh.
function _wfocachepath(zip::AbstractString)
    digest(parts...) = bytes2hex(sha256(join(string.(parts), '\0')))[1:16]
    name = digest(abspath(zip)) * "-" *
           digest(filesize(zip), mtime(zip), _WFOCODEHASH) * ".tsv"
    return joinpath(assetdir(owner = WFOPlantList), "names", name)
end

# Read the columns the reader keeps out of a release's `classification.csv`, extracted from the zip
# through GDAL's zip file system into a temporary file beside the cache, and write them to `cache`
# as a tab-separated table, renamed into place so that a reader never sees part of one; the tables
# cached before from the same file, by an older reader or of an older copy, are removed, so a file
# keeps one.
function _writewfocache(zip::AbstractString, cache::AbstractString)
    mkpath(dirname(cache))
    extracted, io = mktemp(dirname(cache))
    try
        _copyzipentry(io, zip, "classification.csv")
        close(io)
        table = CSV.File(extracted, delim = '\t', quotechar = '"',
                         select = _WFOCOLUMNS, types = String)
        partial = cache * ".part"
        CSV.write(partial, table, delim = '\t')
        mv(partial, cache, force = true)
        prefix = first(basename(cache), 17)
        for f in readdir(dirname(cache))
            (startswith(f, prefix) && endswith(f, ".tsv") &&
             f != basename(cache)) &&
                rm(joinpath(dirname(cache), f), force = true)
        end
    finally
        isopen(io) && close(io)
        rm(extracted, force = true)
    end
    return cache
end

# Copy one entry of a zip file into `io`, a block at a time, through GDAL's zip file system.
function _copyzipentry(io::IO, zip::AbstractString, entry::AbstractString)
    handle = ArchGDAL.GDAL.vsifopenl("/vsizip/" * abspath(zip) * "/" * entry,
                                     "rb")
    handle == C_NULL &&
        error("`$zip` holds no `$entry`, so it is not a WFO Plant List Darwin Core file.")
    buffer = Vector{UInt8}(undef, 1 << 24)
    try
        while true
            n = ArchGDAL.GDAL.vsifreadl(buffer, 1, length(buffer), handle)
            n == 0 && break
            write(io, view(buffer, 1:n))
        end
    finally
        ArchGDAL.GDAL.vsifclosel(handle)
    end
    return io
end

# A cached reduced table read into memory: each row's species and parent found by identifier, a
# synonym naming no row of the file standing for itself.
function _wfotable(cache::AbstractString)
    table = CSV.File(cache, delim = '\t', quotechar = '"', types = String)
    text(v) = String[ismissing(x) ? "" : x for x in v]
    id = text(table.taxonID)
    name = text(table.scientificName)
    authorship = _normauthorship.(text(table.scientificNameAuthorship))
    rank = Symbol.(text(table.taxonRank))
    status = Symbol.(text(table.taxonomicStatus))
    byid = Dict{String, Int32}(x => Int32(i) for (i, x) in enumerate(id))
    acceptedid = text(table.acceptedNameUsageID)
    accepted = Int32[s === :Synonym ? get(byid, a, Int32(i)) : Int32(i)
                     for (i, (s, a)) in enumerate(zip(status, acceptedid))]
    parent = Int32[get(byid, p, Int32(0))
                   for p in text(table.parentNameUsageID)]
    byname = Dict{String, Vector{Int32}}()
    for (i, n) in enumerate(name)
        push!(get!(Vector{Int32}, byname, n), Int32(i))
    end
    return _WFONames(id, name, authorship, rank, status, accepted, parent,
                     byname, byid)
end

# An authorship as it is compared: spaces trimmed and runs of them collapsed.
_normauthorship(s::AbstractString) = join(split(s), " ")

# The rows a name matches: its Accepted and Synonym rows where it has any, else its Unchecked rows
# where those are allowed, else none.
function _wforows(t::_WFONames, name::AbstractString, allow_unchecked::Bool)
    rows = get(t.byname, name, Int32[])
    checked = filter(i -> t.status[i] !== :Unchecked, rows)
    isempty(checked) || return checked
    return allow_unchecked ? copy(rows) : Int32[]
end

# Which species a name means, by the rules `WFOPlantList`'s docstring states.
function _resolvewfo(t::_WFONames, name::String, authorship,
                     allow_unchecked::Bool, prefer_accepted::Bool)
    rows = _wforows(t, name, allow_unchecked)
    if isempty(rows)
        return (name = nothing, id = nothing,
                match = haskey(t.byname, name) ? :unchecked_only : :not_found)
    end
    if !isnothing(authorship) && length(unique(t.accepted[rows])) > 1
        wanted = _normauthorship(authorship)
        picked = filter(i -> t.authorship[i] == wanted, rows)
        isempty(picked) || (rows = picked)
    end
    targets = unique(t.accepted[rows])
    length(targets) == 1 ||
        return (name = name, id = _homonymkey(t, rows, prefer_accepted),
                match = :homonym)
    target = only(targets)
    accepted = t.name[target]
    accepted == name &&
        return (name = name, id = t.id[target],
                match = t.status[target] === :Unchecked ? :unchecked :
                        :accepted)
    # A synonym takes the key its accepted name has on its own, which a homonym joins.
    return (name = accepted,
            id = _wfokey(t, accepted, allow_unchecked, prefer_accepted),
            match = :synonym)
end

# The key a name has on its own: its one species' identifier, or a homonym's joined.
function _wfokey(t::_WFONames, name::AbstractString, allow_unchecked::Bool,
                 prefer_accepted::Bool)
    rows = _wforows(t, name, allow_unchecked)
    targets = unique(t.accepted[rows])
    length(targets) == 1 && return t.id[only(targets)]
    return _homonymkey(t, rows, prefer_accepted)
end

# A homonym's key: its rows' own identifiers, sorted and joined by `+` - only its Accepted rows'
# where those are preferred and it has any.
function _homonymkey(t::_WFONames, rows::AbstractVector{<:Integer},
                     prefer_accepted::Bool)
    if prefer_accepted
        accepted = filter(i -> t.status[i] === :Accepted, rows)
        isempty(accepted) || (rows = accepted)
    end
    return join(sort!(t.id[rows]), "+")
end

# The row of one identifier, refusing one the release does not hold.
function _wforow(t::_WFONames, id::AbstractString, a::WFOPlantList)
    row = get(t.byid, id, nothing)
    isnothing(row) &&
        error("`$id` is no identifier of the WFO Plant List release `$(a.release)`.")
    return row
end

# The rows from a species up to the top of the classification, starting from a synonym's accepted
# species, stopping where a row names no parent (or, defensively, at a cycle).
function _wfopath(t::_WFONames, row::Integer)
    i = t.status[row] === :Synonym ? t.accepted[row] : Int32(row)
    path = Int32[]
    while i != 0 && !(i in path)
        push!(path, i)
        i = t.parent[i]
    end
    return path
end

# The rows several paths share, from the lowest shared one up: their common tail.
function _sharedtop(paths::AbstractVector{<:AbstractVector})
    length(paths) == 1 && return only(paths)
    shortest = minimum(length, paths)
    n = 0
    while n < shortest && allequal(p[end - n] for p in paths)
        n += 1
    end
    return first(paths)[(end - n + 1):end]
end

# The run a release's file goes through, as a spec's track is read.
_spectrack(a::WFOPlantList) = a.track

# A release's file enters a run as the species' identities.
_assetrole(::Type{WFOPlantList}) = :species
