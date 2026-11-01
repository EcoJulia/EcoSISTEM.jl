# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The record of where an input came from, and the reader of the sidecars the fetches write.

"""
    InputRecord(; role, source, code = nothing, path = nothing, url = nothing, request = nothing,
                region = nothing, job = nothing, fetched = nothing, bytes = nothing,
                sha256 = nothing, doi = "", licence = "", version = "", citation = "",
                title = "", primary = nothing, releasedate = "", references = [],
                notices = [], requestpage = "", modified = nothing, history = nothing)

The provenance of one input to a result: a file the package fetched, or a dataset the user handed
it that it never fetched - an occurrence download, a trait table, a phylogeny. One shape for both,
so that everything a result was built from can be listed together and a paper's references written
from the list. [`provenance`](@ref) reads the records the fetches wrote; the user writes one for
anything the package could not have known about.

# Fields

  - `role`: which part of the system the input entered by, a word of the package using the record
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
  - `doi`, `licence`, `version`, `citation`: the source's own, blank where unknown; `citation` is
    the text a paper prints.
  - `title`: the file's own title, the release's with the layer and whatever else tells the file
    apart from its siblings, blank where unknown.
  - `primary`: whether the file is the cited item or one of its files, so that the same bytes can
    be fetched again (`true`) - a file fetched as published from a source whose identifier names
    the data - or not (`false`): a file of a source cited by its paper, an extract a service cut to
    a request, a file the user gave; `nothing` where unknown.
  - `releasedate`: when the release was published, as ISO 8601 text, blank where unknown.
  - `modified`: when the provider last changed the file, as its server said, in UTC; `nothing`
    where unknown.
  - `history`: the file's own global `history` attribute where the fetch replaced it - a netCDF
    file the Climate Data Store converts says there when and in which temporary directory it did
    so, which would make two fetches of the same data differ; `nothing` where nothing was
    replaced.
  - `references`: the works the source's provider asks to be cited, each its `doi` and `text`,
    and `notices`: the statements it asks to accompany the data, each its `strength` (`:required`,
    `:requested` or `:optional`) and `text`. `DD-MMM-YYYY` in a text stands for the date the file
    was fetched and `YYYY` for its year, filled in wherever the record is shown or written.
  - `requestpage`: the address of the documentation's section quoting what the provider asks.
"""
struct InputRecord
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
    doi::String
    licence::String
    version::String
    citation::String
    region::Union{Nothing, String}
    title::String
    primary::Union{Nothing, Bool}
    releasedate::String
    references::Vector{@NamedTuple{doi::String, text::String}}
    notices::Vector{@NamedTuple{strength::Symbol, text::String}}
    requestpage::String
    modified::Union{Nothing, Dates.DateTime}
    history::Union{Nothing, String}

    function InputRecord(; role::Symbol, source::AbstractString,
                         code = nothing, path = nothing, url = nothing,
                         request = nothing, region = nothing, job = nothing,
                         fetched = nothing, bytes = nothing, sha256 = nothing,
                         doi = "", licence = "", version = "", citation = "",
                         title = "", primary::Union{Nothing, Bool} = nothing,
                         releasedate = "", references = [], notices = [],
                         requestpage = "", modified = nothing,
                         history = nothing)
        isnothing(path) || !isabspath(path) ||
            throw(ArgumentError("an `InputRecord` holds a file's name or a relative path, never " *
                                "an absolute one: it is meant to be committed and shared, and " *
                                "`$path` names a directory on one machine."))
        return new(role, String(source), _optstring(code), _optstring(path),
                   _optstring(url),
                   isnothing(request) ? nothing :
                   Dict{String, Any}(string(k) => v for (k, v) in request),
                   _optstring(job), _optdate(fetched),
                   isnothing(bytes) ? nothing : Int(bytes), _optstring(sha256),
                   String(doi), String(licence), String(version),
                   String(citation), _optstring(region), String(title), primary,
                   String(releasedate), _references(references),
                   _notices(notices), String(requestpage), _optdate(modified),
                   _optstring(history))
    end
end

# Value equality and hashing, not the `===` fallback: a record's references, notices and request
# body are containers, which the fallback compares by identity, so two records of one file read on
# two ranks or twice would differ.
function Base.:(==)(a::InputRecord, b::InputRecord)
    return all(getfield(a, f) == getfield(b, f)
               for f in fieldnames(InputRecord))
end

function Base.isequal(a::InputRecord, b::InputRecord)
    return all(isequal(getfield(a, f), getfield(b, f))
               for f in fieldnames(InputRecord))
end

function Base.hash(r::InputRecord, h::UInt)
    return hash(Tuple(getfield(r, f) for f in fieldnames(InputRecord)),
                hash(InputRecord, h))
end

function Base.show(io::IO, r::InputRecord)
    print(io, "InputRecord(", r.role, ", ", repr(r.source))
    isnothing(r.code) || print(io, ", ", repr(r.code))
    isnothing(r.path) || print(io, ", ", repr(r.path))
    isempty(r.doi) || print(io, ", doi = ", repr(r.doi))
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", r::InputRecord)
    println(io, "InputRecord (", r.role, ")")
    for (label, value) in (("source", r.source), ("title", r.title),
        ("code", r.code), ("file", r.path), ("url", r.url),
        ("region", r.region), ("job", r.job), ("fetched", r.fetched),
        ("bytes", r.bytes), ("sha256", r.sha256), ("primary", r.primary),
        ("doi", r.doi), ("licence", r.licence), ("version", r.version),
        ("released", r.releasedate), ("modified", r.modified),
        ("citation", isempty(r.references) ? r.citation : ""))
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
    provenance(path::AbstractString)

Return what is known about where a file came from: the [`InputRecord`](@ref) its provenance
sidecar holds - written beside every file the package fetched or first used, as
`<file>.provenance.toml` - or `nothing` for a file with none. A sidecar of that name that this
machinery did not write is left alone and reports `nothing` too, with a warning saying so, since
its fields mean whatever its writer meant.

# Arguments

  - `path`: the file, not its sidecar.
"""
function provenance(path::AbstractString)
    isfile(sidecarpath(path)) || return nothing
    return _readsidecar(sidecarpath(path))
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
    recordat(path::AbstractString)

Return the record of ours beside the file at `path`, or `nothing` - silently - where there is no
sidecar, or one another program wrote, or one that does not parse. For a read recording its
inputs, where a warning on every build would be noise and a malformed file of someone else's must
not stop it; [`provenance`](@ref) is the asking form, which warns.

# Arguments

  - `path`: the file, not its sidecar.
"""
function recordat(path::AbstractString)
    sidecar = sidecarpath(path)
    isfile(sidecar) || return nothing
    t = TOML.tryparsefile(sidecar)
    (t isa AbstractDict && _isourrecord(t)) || return nothing
    return _recordfromtoml(t)
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

# A sidecar's TOML as an `InputRecord`. A sidecar another program wrote under the same name is not
# read: its `role` may be prose and its `source` a sentence, so a record built from it would be
# wrong in ways no reader could see.
function _readsidecar(sidecar::AbstractString)
    t = TOML.parsefile(sidecar)
    _isourrecord(t) || begin
        @warn "`$(basename(sidecar))` was not written by Provenancer (writer: "*
              "$(repr(get(t, "writer", nothing)))), so it is left alone and reports nothing; "*
              "delete or rename it to have the file recorded on its next use." maxlog=1 _id=Symbol(sidecar)
        return nothing
    end
    return _recordfromtoml(t)
end

# A parsed sidecar of ours as an `InputRecord`: the keys the writers use, any it does not know left
# behind. The source is read from `source`, or from `dataset`, which earlier records called it.
function _recordfromtoml(t::AbstractDict)
    field(k) = get(t, k, nothing)
    return InputRecord(role = Symbol(something(field("role"), "habitat")),
                       source = something(field("source"), field("dataset"),
                                          "file"),
                       code = field("code"), path = field("file"),
                       url = something(field("final_url"), field("url"),
                                       Some(nothing)),
                       request = field("request"), job = field("job"),
                       fetched = field("fetched"), bytes = field("bytes"),
                       sha256 = field("sha256"),
                       doi = something(field("doi"), ""),
                       licence = something(field("licence"), ""),
                       version = something(field("version"), ""),
                       citation = something(field("citation"), ""),
                       region = field("region"),
                       title = something(field("title"), ""),
                       primary = field("primary"),
                       releasedate = something(field("release_date"), ""),
                       references = something(field("references"), []),
                       notices = something(field("notices"), []),
                       requestpage = something(field("request_page"), ""),
                       modified = field("modified"),
                       history = field("history"))
end

# An optional text field, `nothing` kept, anything else as a `String`.
_optstring(x) = isnothing(x) ? nothing : string(x)

# An optional time field: a `DateTime` as it is, ISO 8601 text with or without its `Z` parsed.
_optdate(::Nothing) = nothing

_optdate(t::Dates.DateTime) = t

_optdate(s::AbstractString) = Dates.DateTime(chopsuffix(String(s), "Z"))

# References in the one shape a record holds them, from named tuples or a sidecar's tables.
function _references(refs::AbstractVector)
    return @NamedTuple{doi::String, text::String}[(doi = String(_entry(r, :doi)),
                                                   text = String(_entry(r,
                                                                        :text)))
                                                  for r in refs]
end

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

# What a record says its provider asks, a line each: the references to cite and the notices to
# carry, placeholders filled from its fetch time, and where the documentation quotes the request. A
# record with no references but a citation of its own cites that.
function _requestlines(r::InputRecord)
    cites = isempty(r.references) ?
            (isempty(r.citation) ? String[] : [r.citation]) :
            [_fillplaceholders(ref.text, r.fetched) for ref in r.references]
    lines = ["cite    " * c for c in cites]
    append!(lines,
            ["notice  ($(n.strength)) " * _fillplaceholders(n.text, r.fetched)
             for n in r.notices])
    isempty(r.requestpage) || push!(lines, "asks    " * r.requestpage)
    return lines
end

# One input's line under its source in a provenance listing: its role, its layer where it has one,
# and its file, or where it was fetched from where it has no file.
function _recordline(r::InputRecord)
    parts = (string(r.role), something(r.code, ""),
             something(r.path, r.url, ""))
    return join(filter(!isempty, parts), "  ")
end
