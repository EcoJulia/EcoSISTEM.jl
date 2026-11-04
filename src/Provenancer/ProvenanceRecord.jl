# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The record of where a file came from, and the reader of the sidecars the fetches write.

"""
    ProvenanceRecord(; role, source, code = nothing, path = nothing, url = nothing,
                     request = nothing, region = nothing, job = nothing, fetched = nothing,
                     bytes = nothing, sha256 = nothing, doi = nothing, licence = nothing,
                     version = nothing, citation = nothing, title = nothing, primary = nothing,
                     releasedate = nothing, references = [], notices = [],
                     requestpage = nothing, modified = nothing, history = nothing,
                     identifies = nothing, key = nothing, namespace = nothing,
                     catalogue = nothing, deposit = DepositSharesIdentifier(),
                     route = nothing, location = nothing, sidecar = nothing)

The provenance of one file: a file the package fetched, a dataset the user handed it that it never
fetched - an occurrence download, a trait table, a phylogeny - or a file a run of theirs made. One
shape for all, so that everything a result was built from can be listed together and a paper's
references written from the list. [`provenance`](@ref) reads the records the fetches wrote;
[`record_provenance`](@ref) writes one for anything the package could not have known about.

# Fields

  - `role`: which part of the system the file entered by, a word of the package using the record
    (EcoSISTEM's are `:habitat` for a layer, `:region` for a shape or a named region's outline,
    `:species`, `:phylogeny`, `:dispersal`, `:demography` or `:intervention` for what the user
    attached, `:abundance` for records a starting population was placed from, `:state` for a
    saved run a forward run starts from, `:software` for the package itself). Any symbol is
    accepted here; the package using the record says which it takes.
  - `source`: the source's name - a catalogue key such as `TwentyCR` or `WorldClim{BioClim}`, or
    the user's own name for a dataset of theirs.
  - `code`: which layer of the source, where the input is one.
  - `path`: the file's own name, or a path relative to wherever its record sits - never an absolute
    one, since a record is meant to be committed and shared.
  - `url`: where it was fetched from, where it was.
  - `request`: the request body sent to a service, as it was sent.
  - `region`: the named region a request's area was taken from, where it was given as one.
  - `job`: the identifier a service gave the request.
  - `fetched`: when it was fetched, in UTC.
  - `bytes`, `sha256`: the file's size and checksum at fetch, which is what tells a truncated or
    replaced copy from the original.
  - `doi`, `licence`, `version`, `citation`: the source's own, `nothing` where unknown, as every
    text field is, and blank text given is taken as `nothing`; `citation` is the text a paper
    prints.
  - `title`: the file's own title, the release's with the layer and whatever else tells the file
    apart from its siblings, `nothing` where unknown.
  - `primary`: whether the file is the cited item or one of its files (`true`) - the version of
    record of a source whose identifier names the data - or not (`false`): a file of a source cited
    by its paper, an extract a service cut to a request, a file of any other origin; `nothing`
    where unknown.
  - `releasedate`: when the release was published, as ISO 8601 text, `nothing` where unknown.
  - `modified`: when the provider last changed the file, as its server said, in UTC; `nothing`
    where unknown.
  - `history`: the file's own global `history` attribute where the fetch replaced it - a netCDF
    file the Climate Data Store converts says there when and in which temporary directory it did
    so, which would make two fetches of the same data differ; `nothing` where nothing was
    replaced.
  - `references`: the works to cite with the source - its provider's request, or works cited
    alongside it - each its `doi` and `text`, and `notices`: the statements it asks to accompany the data, each its `strength` (`:required`,
    `:requested` or `:optional`) and `text`. `DD-MMM-YYYY` in a text stands for the date the file
    was fetched and `YYYY` for its year, filled in wherever the record is shown or written.
  - `requestpage`: the address of the documentation's section quoting what the provider asks.
  - `identifies`: what the source's identifier names - `:data`, `:paper` or `:other`, the values
    of [`IDENTIFIES`](@ref Provenancer.IDENTIFIES) - and `key`, the key of the citation naming a
    source with no DOI, and `namespace`, the data pipeline namespace its files are registered in,
    `nothing` where unknown: what a registration needs of the source beyond its entries. Where
    `identifies` is not given, a file called `primary` identifies the data and any other is
    `:other`.
  - `catalogue`: the name of the catalogue the file was fetched from, `nothing` for none: with a
    fetch time, the mark of a file the package fetched from its catalogue
    ([`iscatalogued`](@ref Provenancer.iscatalogued)).
  - `deposit`: what the source's identifier stands for in a data pipeline's registry, an
    [`AbstractDeposit`](@ref Provenancer.AbstractDeposit); [`DepositSharesIdentifier`](@ref) by
    default, and for a record written before deposits were recorded.
  - `route`: how the file came to be, an [`AbstractFileRoute`](@ref Provenancer.AbstractFileRoute),
    or `nothing` for a record written before routes were recorded.
  - `location`: the file's path as it was when the record was read, and `sidecar`, the path the
    record was read from, or `nothing` for a record built in memory. Neither is written: the
    sidecar and the account a result writes carry the file's name alone.
"""
struct ProvenanceRecord
    role::Symbol
    source::String
    code::Union{Nothing, String}
    path::Union{Nothing, String}
    url::Union{Nothing, String}
    request::Union{Nothing, Dict{String, Any}}
    job::Union{Nothing, String}
    fetched::Union{Nothing, Dates.DateTime}
    bytes::Union{Nothing, Int}
    sha256::Union{Nothing, String}
    doi::Union{Nothing, String}
    licence::Union{Nothing, String}
    version::Union{Nothing, String}
    citation::Union{Nothing, String}
    region::Union{Nothing, String}
    title::Union{Nothing, String}
    primary::Union{Nothing, Bool}
    releasedate::Union{Nothing, String}
    references::Vector{@NamedTuple{doi::String, text::String}}
    notices::Vector{@NamedTuple{strength::Symbol, text::String}}
    requestpage::Union{Nothing, String}
    modified::Union{Nothing, Dates.DateTime}
    history::Union{Nothing, String}
    identifies::Symbol
    key::Union{Nothing, String}
    namespace::Union{Nothing, String}
    catalogue::Union{Nothing, String}
    deposit::AbstractDeposit
    route::Union{Nothing, AbstractFileRoute}
    location::Union{Nothing, String}
    sidecar::Union{Nothing, String}

    function ProvenanceRecord(; role::Symbol, source::AbstractString,
                              code = nothing, path = nothing, url = nothing,
                              request = nothing, region = nothing,
                              job = nothing, fetched = nothing, bytes = nothing,
                              sha256 = nothing, doi = nothing,
                              licence = nothing,
                              version = nothing, citation = nothing,
                              title = nothing,
                              primary::Union{Nothing, Bool} = nothing,
                              releasedate = nothing, references = [],
                              notices = [], requestpage = nothing,
                              modified = nothing, history = nothing,
                              identifies::Union{Nothing, Symbol} = nothing,
                              key = nothing, namespace = nothing,
                              catalogue = nothing,
                              deposit::AbstractDeposit = DepositSharesIdentifier(),
                              route::Union{Nothing, AbstractFileRoute} = nothing,
                              location = nothing, sidecar = nothing)
        isnothing(path) || !isabspath(path) ||
            throw(ArgumentError("a `ProvenanceRecord` holds a file's name or a relative path, " *
                                "never an absolute one: it is meant to be committed and " *
                                "shared, and `$path` names a directory on one machine."))
        names = something(identifies, primary === true ? :data : :other)
        names in IDENTIFIES ||
            throw(ArgumentError("`identifies` must be one of $(join(IDENTIFIES, ", ")); got " *
                                "`$names`."))
        return new(role, String(source), _optstring(code), _optstring(path),
                   _optstring(url),
                   isnothing(request) ? nothing :
                   Dict{String, Any}(string(k) => v for (k, v) in request),
                   _optstring(job), _optdate(fetched),
                   isnothing(bytes) ? nothing : Int(bytes), _optstring(sha256),
                   _opttext(doi), _opttext(licence), _opttext(version),
                   _opttext(citation), _optstring(region), _opttext(title),
                   primary, _opttext(releasedate), _references(references),
                   _notices(notices), _opttext(requestpage), _optdate(modified),
                   _optstring(history), names, _opttext(key),
                   _opttext(namespace),
                   _opttext(catalogue), deposit, route, _optstring(location),
                   _optstring(sidecar))
    end
end

# Value equality and hashing, not the `===` fallback: a record's references, notices and request
# body are containers, which the fallback compares by identity, so two records of one file read on
# two ranks or twice would differ.
function Base.:(==)(a::ProvenanceRecord, b::ProvenanceRecord)
    return all(getfield(a, f) == getfield(b, f)
               for f in fieldnames(ProvenanceRecord))
end

function Base.isequal(a::ProvenanceRecord, b::ProvenanceRecord)
    return all(isequal(getfield(a, f), getfield(b, f))
               for f in fieldnames(ProvenanceRecord))
end

function Base.hash(r::ProvenanceRecord, h::UInt)
    return hash(Tuple(getfield(r, f) for f in fieldnames(ProvenanceRecord)),
                hash(ProvenanceRecord, h))
end

function Base.show(io::IO, r::ProvenanceRecord)
    print(io, "ProvenanceRecord(", r.role, ", ", repr(r.source))
    isnothing(r.code) || print(io, ", ", repr(r.code))
    isnothing(r.path) || print(io, ", ", repr(r.path))
    isnothing(r.doi) || print(io, ", doi = ", repr(r.doi))
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", r::ProvenanceRecord)
    println(io, "ProvenanceRecord (", r.role, ")")
    for (label, value) in (("source", r.source), ("title", r.title),
        ("code", r.code), ("file", r.path), ("location", r.location),
        ("sidecar", r.sidecar), ("route", r.route), ("url", r.url),
        ("region", r.region), ("job", r.job), ("fetched", r.fetched),
        ("bytes", r.bytes), ("sha256", r.sha256), ("primary", r.primary),
        ("doi", r.doi), ("licence", r.licence), ("version", r.version),
        ("released", r.releasedate), ("modified", r.modified),
        ("namespace", r.namespace), ("catalogue", r.catalogue),
        ("deposit",
        r.deposit isa DepositSharesIdentifier ? nothing : r.deposit),
        ("citation", _citesown(r) ? r.citation : nothing))
        (isnothing(value) || (value isa AbstractString && isempty(value))) &&
            continue
        println(io, "  ", rpad(label, 9), value)
    end
    for line in _requestlines(r)
        println(io, "  ", line)
    end
    return nothing
end

# == Functions ======================================================================================

"""
    provenance(path::AbstractString; sidecar = sidecarpath(path))

Return what is known about where a file came from: the [`ProvenanceRecord`](@ref) its provenance
sidecar holds - written beside every file the package fetched or first used, as
`<file>.provenance.toml`, or where [`record_provenance`](@ref) was told to put it - with `path`
and `sidecar` as its `location` and `sidecar`, or `nothing` for a file with none. A sidecar of
that name that this machinery did not write is left alone and reports `nothing` too, with a
warning saying so, since its fields mean whatever its writer meant.

# Arguments

  - `path`: the file, not its sidecar.
  - `sidecar`: where its record is, the file's `.provenance.toml` beside it by default.
"""
function provenance(path::AbstractString;
                    sidecar::AbstractString = sidecarpath(path))
    isfile(sidecar) || return nothing
    return _readsidecar(path, sidecar)
end

"""
    NOTICE_STRENGTHS

The strengths a provider's notice may carry in a record: `:required`, `:requested` or `:optional`.
"""
const NOTICE_STRENGTHS = (:required, :requested, :optional)

"""
    sidecarpath(path::AbstractString)

Return the path of the provenance sidecar of the file at `path`: its own name with
`.provenance.toml` after it, beside it.

# Arguments

  - `path`: the file, not its sidecar.
"""
sidecarpath(path::AbstractString) = path * ".provenance.toml"

"""
    recordat(path::AbstractString; sidecar = sidecarpath(path))

Return the record of ours beside the file at `path`, or `nothing` - silently - where there is no
sidecar, or one another program wrote, or one that does not parse. For a read recording its
inputs, where a warning on every build would be noise and a malformed file of someone else's must
not stop it; [`provenance`](@ref) is the asking form, which warns.

# Arguments

  - `path`: the file, not its sidecar.
  - `sidecar`: where its record is, the file's `.provenance.toml` beside it by default.
"""
function recordat(path::AbstractString;
                  sidecar::AbstractString = sidecarpath(path))
    isfile(sidecar) || return nothing
    t = TOML.tryparsefile(sidecar)
    (t isa AbstractDict && _isourrecord(t)) || return nothing
    return _recordfromtoml(t, location = path, sidecar = sidecar)
end

# What this module writes as a record's `writer`: its full name and the version of the package it
# is part of, which is what tells its records from a same-named TOML another program wrote.
function _writer()
    return join(string.(Base.fullname(@__MODULE__)), ".") * " " *
           string(something(pkgversion(@__MODULE__), ""))
end

# Whether a parsed sidecar is one this machinery wrote, by its `writer`: this module under any
# package it is part of, or `EcoSISTEM` itself, which wrote the records before the machinery had a
# name of its own and whose records must go on being read.
function _isourrecord(t::AbstractDict)
    writer = string(get(t, "writer", ""))
    name = first(split(writer, ' ', limit = 2))
    return endswith(name, "Provenancer") || startswith(writer, "EcoSISTEM")
end

# The sidecar at `sidecar` of the file at `path` as a `ProvenanceRecord`. A sidecar another program
# wrote under the same name is not read: its `role` may be prose and its `source` a sentence, so a
# record built from it would be wrong in ways no reader could see.
function _readsidecar(path::AbstractString, sidecar::AbstractString)
    t = TOML.parsefile(sidecar)
    _isourrecord(t) || begin
        @warn "`$(basename(sidecar))` was not written by Provenancer (writer: "*
              "$(repr(get(t, "writer", nothing)))), so it is left alone and reports nothing; "*
              "delete or rename it to have the file recorded on its next use." maxlog=1 _id=Symbol(sidecar)
        return nothing
    end
    return _recordfromtoml(t, location = path, sidecar = sidecar)
end

# A parsed sidecar of ours as a `ProvenanceRecord`, read from `sidecar` for the file at `location`:
# the keys the writers use, any it does not know left behind. The source is read from `source`, or
# from `dataset`, which earlier records called it.
function _recordfromtoml(t::AbstractDict; location, sidecar)
    field(k) = get(t, k, nothing)
    identifies = field("identifies")
    return ProvenanceRecord(role = Symbol(something(field("role"), "habitat")),
                            source = something(field("source"),
                                               field("dataset"), "file"),
                            code = field("code"), path = field("file"),
                            url = something(field("final_url"), field("url"),
                                            Some(nothing)),
                            request = field("request"), job = field("job"),
                            fetched = field("fetched"), bytes = field("bytes"),
                            sha256 = field("sha256"),
                            doi = field("doi"), licence = field("licence"),
                            version = field("version"),
                            citation = field("citation"),
                            region = field("region"), title = field("title"),
                            primary = field("primary"),
                            releasedate = field("release_date"),
                            references = something(field("references"), []),
                            notices = something(field("notices"), []),
                            requestpage = field("request_page"),
                            modified = field("modified"),
                            history = field("history"),
                            identifies = isnothing(identifies) ? nothing :
                                         Symbol(identifies),
                            key = field("key"), namespace = field("namespace"),
                            catalogue = field("catalogue"),
                            deposit = _depositfrom(field("deposit"),
                                                   field("deposit_source")),
                            route = _routefrom(field("route"),
                                               field("route_identifier")),
                            location = location, sidecar = sidecar)
end

# An optional text field, `nothing` kept, anything else as a `String`.
_optstring(x) = isnothing(x) ? nothing : string(x)

# An optional text field whose blank is unknown too: `nothing` kept, blank text made `nothing`,
# anything else as a `String`.
function _opttext(x)
    isnothing(x) && return nothing
    s = string(x)
    return isempty(s) ? nothing : s
end

# An optional time field: a `DateTime` as it is, ISO 8601 text with or without its `Z` parsed.
_optdate(::Nothing) = nothing

_optdate(t::Dates.DateTime) = t

_optdate(s::AbstractString) = Dates.DateTime(chopsuffix(String(s), "Z"))

# References in the one shape a record holds them, from named tuples, a sidecar's tables or the
# facts of the works themselves.
function _references(refs::AbstractVector)
    return @NamedTuple{doi::String, text::String}[_reference(r) for r in refs]
end

# One reference in the shape a record holds it, from a named tuple or a sidecar's table.
_reference(r) = (doi = String(_entry(r, :doi)), text = String(_entry(r, :text)))

# Notices likewise.
function _notices(notices::AbstractVector)
    return @NamedTuple{strength::Symbol, text::String}[(strength = _strength(_entry(n,
                                                                                    :strength)),
                                                        text = String(_entry(n,
                                                                             :text)))
                                                       for n in notices]
end

# A notice's strength, one a provider can ask with.
function _strength(s)
    strength = Symbol(s)
    strength in NOTICE_STRENGTHS ||
        throw(ArgumentError("a notice's `strength` must be one of " *
                            "$(join(NOTICE_STRENGTHS, ", ")); got `$strength`."))
    return strength
end

# One entry of a reference or a notice, from a named tuple or from a sidecar's table.
_entry(x::NamedTuple, key::Symbol) = x[key]

_entry(x::AbstractDict, key::Symbol) = x[string(key)]

# A provider's text with its placeholders filled from when the file was fetched: `DD-MMM-YYYY` with
# the date, then `YYYY` with the year. Left as it is where the fetch time is not known.
_fillplaceholders(text::AbstractString, ::Nothing) = String(text)

function _fillplaceholders(text::AbstractString, fetched::Dates.DateTime)
    dated = replace(text, "DD-MMM-YYYY" => Dates.format(fetched, "dd-u-yyyy"))
    return replace(dated, "YYYY" => string(Dates.year(fetched)))
end

# Whether a record's own citation is cited: not where it has none; not where it is composed from its
# references - their texts joined with a space, as a catalogue row's is, in sidecars old and new -
# nor where a reference carries the record's own DOI, that reference being the source's citation;
# otherwise yes, before the references, which are cited alongside it.
function _citesown(r::ProvenanceRecord)
    isnothing(r.citation) && return false
    isempty(r.references) && return true
    r.citation == join((ref.text for ref in r.references), " ") && return false
    isnothing(r.doi) && return true
    return !any(ref -> lowercase(ref.doi) == lowercase(r.doi), r.references)
end

# What a record says is to be cited and carried, a line each: its own citation where `_citesown`
# says so, its deposit's where that is cited apart, the references, the notices, placeholders filled
# from its fetch time, and where the documentation quotes the request.
function _requestlines(r::ProvenanceRecord)
    cites = _citesown(r) ? [r.citation] : String[]
    append!(cites, _depositcitations(r.deposit))
    append!(cites,
            [_fillplaceholders(ref.text, r.fetched) for ref in r.references])
    lines = ["cite    " * c for c in cites]
    append!(lines,
            ["notice  ($(n.strength)) " * _fillplaceholders(n.text, r.fetched)
             for n in r.notices])
    isnothing(r.requestpage) || push!(lines, "asks    " * r.requestpage)
    return lines
end

# One input's line under its source in a provenance listing: its role, its layer where it has one,
# and its file, or where it was fetched from where it has no file.
function _recordline(r::ProvenanceRecord)
    parts = (string(r.role), something(r.code, ""),
             something(r.path, r.url, ""))
    return join(filter(!isempty, parts), "  ")
end
