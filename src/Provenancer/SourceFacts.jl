# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What is known of a source - a dataset, a paper, a download a service built - that a file's record
# and its registration are written from, whether the facts come from a catalogue row, from the
# metadata a DOI resolves to, or from the user's own hand.

"""
    SourceFacts(; source, doi = nothing, url = nothing, title = nothing, publisher = nothing,
                authors = [], identifiers = Dict(), releasedate = nothing, version = nothing,
                description = nothing, licence = nothing, identifies = :other,
                citation = nothing, references = [], notices = [], requestpage = nothing,
                key = nothing, namespace = nothing)
    SourceFacts(facts::SourceFacts; kw...)
    SourceFacts(record::ProvenanceRecord)

What is known of the source a file comes from, from which the file's provenance record and its
registration with a data pipeline are written: the facts a catalogue row holds of a dataset, the
facts the metadata of a DOI gives ([`lookup_doi`](@ref)), or facts the user writes for a source
that has neither. `nothing` and an empty list mean unknown throughout, and blank text given for a
field is taken as `nothing`. The second form is `facts`
with the fields named by keyword replaced and every other as it was: a looked-up source's facts
with the citation its provider asks for, or the works it asks to be cited, written in. The third
is what a file's [`ProvenanceRecord`](@ref) keeps of its source, so that the file can be
registered from its record alone.

# Fields

  - `source`: the source's name, as a record's `source` carries it - a catalogue key such as
    `WorldClim{BioClim}`, a title, or the user's own name for it.
  - `doi`: its DOI, bare (`10.5285/...`), where it has one.
  - `url`: its address: a landing page, or where its files are described.
  - `title`: its title.
  - `publisher`: who published it.
  - `authors`: its authors in order, a person as given names then family name, an organisation as
    it is written; and `identifiers`: the ORCID or ROR identifier of each author that has one, as a
    URL, by the author's name.
  - `releasedate`: when the release was published, as ISO 8601 text as far as it is known
    (`"2022"`, `"2022-06"`, `"2022-06-22"`).
  - `version`: the release's version as the source states it.
  - `description`: what the source says of itself, an abstract.
  - `licence`: the licence the data are under, as an identifier in the SPDX licence list where
    the source gives one, else the address of its terms.
  - `identifies`: what the source's identifier names, one of [`IDENTIFIES`](@ref
    Provenancer.IDENTIFIES): `:data`, the data themselves, so that a file fetched as published is
    the cited item or one of its files; `:paper`, a work describing them; `:software`; or
    `:other`.
  - `citation`: the text a paper prints for it, where the source gives one; composed from the
    authors, date, title, publisher and DOI where it does not.
  - `references`: the works its provider asks to be cited, each its `doi` and `text`; `notices`:
    the statements it asks to accompany the data, each its `strength` (one of
    [`NOTICE_STRENGTHS`](@ref Provenancer.NOTICE_STRENGTHS)) and `text`; and `requestpage`: the
    address of the documentation quoting what it asks.
  - `key`: the key its citation is filed under where it has no DOI, for a registry that takes a
    citation by key.
  - `namespace`: the data pipeline namespace its files are registered in - ASCII letters, digits,
    `_` and `-` - or `nothing` for the namespace a run gives as its default.
"""
struct SourceFacts
    source::String
    doi::Union{Nothing, String}
    url::Union{Nothing, String}
    title::Union{Nothing, String}
    publisher::Union{Nothing, String}
    authors::Vector{String}
    identifiers::Dict{String, String}
    releasedate::Union{Nothing, String}
    version::Union{Nothing, String}
    description::Union{Nothing, String}
    licence::Union{Nothing, String}
    identifies::Symbol
    citation::Union{Nothing, String}
    references::Vector{@NamedTuple{doi::String, text::String}}
    notices::Vector{@NamedTuple{strength::Symbol, text::String}}
    requestpage::Union{Nothing, String}
    key::Union{Nothing, String}
    namespace::Union{Nothing, String}

    function SourceFacts(; source::AbstractString, doi = nothing, url = nothing,
                         title = nothing, publisher = nothing,
                         authors = String[],
                         identifiers = Dict{String, String}(),
                         releasedate = nothing, version = nothing,
                         description = nothing, licence = nothing,
                         identifies::Symbol = :other, citation = nothing,
                         references = [], notices = [], requestpage = nothing,
                         key = nothing, namespace = nothing)
        identifies in IDENTIFIES ||
            throw(ArgumentError("`identifies` must be one of $(join(IDENTIFIES, ", ")); got " *
                                "`$identifies`."))
        ns = _opttext(namespace)
        isnothing(ns) || occursin(r"^[A-Za-z0-9_-]+$", ns) ||
            throw(ArgumentError("`namespace` names a data pipeline namespace, in ASCII letters, " *
                                "digits, `_` and `-`; got `$ns`."))
        return new(String(source), _opttext(doi), _opttext(url),
                   _opttext(title), _opttext(publisher),
                   String[string(a) for a in authors],
                   Dict{String, String}(string(k) => string(v)
                                        for (k, v) in identifiers),
                   _opttext(releasedate), _opttext(version),
                   _opttext(description), _opttext(licence), identifies,
                   _opttext(citation), _references(references),
                   _notices(notices), _opttext(requestpage), _opttext(key), ns)
    end
end

# The facts with the named fields replaced, through the keyword constructor, which checks them.
function SourceFacts(facts::SourceFacts; kw...)
    fields = NamedTuple{fieldnames(SourceFacts)}(Tuple(getfield(facts, f)
                                                       for f in fieldnames(SourceFacts)))
    return SourceFacts(; merge(fields, NamedTuple(kw))...)
end

# The facts a record keeps of its source: the entries `sourceentries` wrote and what
# `_factsentries` kept.
function SourceFacts(record::ProvenanceRecord)
    return SourceFacts(source = record.source, doi = record.doi,
                       url = record.url, title = record.title,
                       releasedate = record.releasedate,
                       version = record.version, licence = record.licence,
                       identifies = record.identifies,
                       citation = record.citation,
                       references = record.references, notices = record.notices,
                       requestpage = record.requestpage, key = record.key,
                       namespace = record.namespace)
end

# Value equality and hashing: the authors, identifiers, references and notices are containers,
# which the fallback compares by identity.
function Base.:(==)(a::SourceFacts, b::SourceFacts)
    return all(getfield(a, f) == getfield(b, f)
               for f in fieldnames(SourceFacts))
end

function Base.hash(f::SourceFacts, h::UInt)
    for name in fieldnames(SourceFacts)
        h = hash(getfield(f, name), h)
    end
    return h
end

function Base.show(io::IO, f::SourceFacts)
    print(io, "SourceFacts(", repr(f.source))
    isnothing(f.doi) || print(io, ", doi ", f.doi)
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", f::SourceFacts)
    println(io, "SourceFacts ", repr(f.source), " (", f.identifies, ")")
    for name in (:title, :doi, :url, :publisher, :releasedate, :version,
        :licence, :citation, :key, :namespace)
        value = getfield(f, name)
        isnothing(value) || println(io, "  ", rpad(string(name), 12), value)
    end
    isempty(f.authors) ||
        println(io, "  ", rpad("authors", 12), join(f.authors, "; "))
    isnothing(f.description) ||
        println(io, "  ", rpad("description", 12),
                length(f.description) > 70 ?
                first(f.description, 67) * "..." : f.description)
    for (name, list) in ((:references, f.references), (:notices, f.notices))
        isempty(list) || println(io, "  ", rpad(string(name), 12), length(list))
    end
    return nothing
end

"""
    IDENTIFIES

What a source's identifier may name: `:data`, the data themselves; `:paper`, a work describing
them; `:software`; or `:other`.
"""
const IDENTIFIES = (:data, :paper, :software, :other)

# == Functions ======================================================================================

"""
    isprimary(facts::SourceFacts, route::AbstractFileRoute)

Return whether a file of the source is the cited item or one of its files under the source's
identifier - what a data pipeline registers as primary. The source's version of record
([`VersionOfRecord`](@ref)) is, where the identifier names the data; a file a service generated to
a request ([`GeneratedOnRequest`](@ref)) is only where the service gave it an identifier of its own
and the facts are that file's; a file of any other origin ([`OtherOrigin`](@ref)) is not, since
nothing is claimed of how its bytes relate to the source's.

# Arguments

  - `facts`: what is known of the source.
  - `route`: the route the file came by.
"""
function isprimary(facts::SourceFacts, ::VersionOfRecord)
    return facts.identifies === :data
end

function isprimary(facts::SourceFacts, route::GeneratedOnRequest)
    return !isnothing(route.identifier) && facts.identifies === :data
end

isprimary(::SourceFacts, ::AbstractFileRoute) = false

"""
    sourceentries(facts::SourceFacts, route::AbstractFileRoute; title = facts.title,
                  code = nothing, catalogue = nothing)

Return what a file's provenance record carries of its source, as a dictionary of the sidecar's
keys, to be merged into the record of the file itself: the `source`, the file's `title`, the
`route` it came by, with the identifier a service gave it where it gave one, the `catalogue`
fetching it where one is, whether it is `primary` for that route ([`isprimary`](@ref Provenancer.isprimary)), and the DOI, licence,
version, release date, citation, references, notices and request page where known, each under
the key the sidecar uses. The citation is composed from the authors, date, title,
publisher and DOI where the facts give none of their own and have a DOI and authors to compose
one from.

# Arguments

  - `facts`: what is known of the source.
  - `route`: the route the file came by.
  - `title`: the file's own title, the source's by default; a source of several files gives
    each its own.
  - `code`: which layer of the source the file holds, where it is one.
  - `catalogue`: the name of the catalogue the file is fetched from, given by a package's fetch
    from a catalogue of its own - the package's name, in ASCII letters, digits, `_` and `-` - and
    `nothing` for any other file: with a fetch time, the mark of a file the package fetched from
    its catalogue ([`iscatalogued`](@ref Provenancer.iscatalogued)).
"""
function sourceentries(facts::SourceFacts, route::AbstractFileRoute;
                       title::Union{Nothing, AbstractString} = facts.title,
                       code = nothing,
                       catalogue::Union{Nothing, AbstractString} = nothing)
    catalogue = _opttext(catalogue)
    isnothing(catalogue) || occursin(r"^[A-Za-z0-9_-]+$", catalogue) ||
        throw(ArgumentError("`catalogue` names the catalogue the file is fetched from, in ASCII " *
                            "letters, digits, `_` and `-`; got `$catalogue`."))
    record = Dict{String, Any}("source" => facts.source,
                               "primary" => isprimary(facts, route))
    merge!(record, _routeentries(route))
    for (key, value) in (("title", title), ("catalogue", catalogue),
        ("doi", facts.doi), ("licence", facts.licence),
        ("version", facts.version), ("release_date", facts.releasedate),
        ("citation", _citationtext(facts)),
        ("request_page", facts.requestpage))
        isnothing(value) || (record[key] = String(value))
    end
    isempty(facts.references) ||
        (record["references"] = [Dict("doi" => r.doi, "text" => r.text)
                                 for r in facts.references])
    isempty(facts.notices) ||
        (record["notices"] = [Dict("strength" => string(n.strength),
                                   "text" => n.text)
                              for n in facts.notices])
    isnothing(code) || (record["code"] = string(code))
    return record
end

# The citation a record carries: the source's own, else one composed as a reference list writes a
# dataset - authors, year, title with its version, publisher, DOI - where there is a DOI and there
# are authors to compose it from, else none.
function _citationtext(facts::SourceFacts)
    isnothing(facts.citation) || return facts.citation
    (isnothing(facts.doi) || isempty(facts.authors)) && return nothing
    year = !isnothing(facts.releasedate) && length(facts.releasedate) >= 4 ?
           " (" * first(facts.releasedate, 4) * ")" : ""
    title = isnothing(facts.title) ? "" :
            " " * facts.title *
            (isnothing(facts.version) ? "" :
             " (version " * facts.version * ")") *
            "."
    publisher = isnothing(facts.publisher) ? "" : " " * facts.publisher * "."
    return join(facts.authors, ", ") * year * "." * title * publisher *
           " https://doi.org/" * facts.doi
end

# What a record of a file recorded from facts keeps of them beyond its entries, so that the file
# can be registered later from its record alone: what the identifier identifies, the citation key
# where there is no DOI, and the namespace, each only where it says something.
function _factsentries(facts::SourceFacts)
    record = Dict{String, Any}("identifies" => string(facts.identifies))
    isnothing(facts.key) || (record["key"] = facts.key)
    isnothing(facts.namespace) || (record["namespace"] = facts.namespace)
    return record
end

# The facts of a file's source as its record of ours at `sidecar` gives them back, for registering
# the file from its record, or `nothing` where there is no record of ours there.
function _factsat(sidecar::AbstractString)
    isfile(sidecar) || return nothing
    t = TOML.tryparsefile(sidecar)
    (t isa AbstractDict && _isourrecord(t)) || return nothing
    return SourceFacts(_recordfromtoml(t, location = nothing,
                                       sidecar = sidecar))
end
