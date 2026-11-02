# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What is known of a source - a dataset, a paper, a download a service built - that a file's record
# and its registration are written from, whether the facts come from a catalogue row, from the
# metadata a DOI resolves to, or from the user's own hand.

"""
    SourceFacts(; source, doi = "", url = "", title = "", publisher = "", authors = [],
                identifiers = Dict(), releasedate = "", version = "", description = "",
                licence = "", identifies = :other, citation = "", references = [],
                notices = [], requestpage = "", key = "")
    SourceFacts(facts::SourceFacts; kw...)

What is known of the source a file comes from, from which the file's provenance record and its
registration with a data pipeline are written: the facts a catalogue row holds of a dataset, the
facts the metadata of a DOI gives ([`lookup_doi`](@ref)), or facts the user writes for a source
that has neither. Blank text and empty lists mean unknown throughout. The second form is `facts`
with the fields named by keyword replaced and every other as it was: a looked-up source's facts
with the citation its provider asks for, or the works it asks to be cited, written in.

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
"""
struct SourceFacts
    source::String
    doi::String
    url::String
    title::String
    publisher::String
    authors::Vector{String}
    identifiers::Dict{String, String}
    releasedate::String
    version::String
    description::String
    licence::String
    identifies::Symbol
    citation::String
    references::Vector{@NamedTuple{doi::String, text::String}}
    notices::Vector{@NamedTuple{strength::Symbol, text::String}}
    requestpage::String
    key::String

    function SourceFacts(; source::AbstractString, doi = "", url = "",
                         title = "", publisher = "", authors = String[],
                         identifiers = Dict{String, String}(), releasedate = "",
                         version = "", description = "", licence = "",
                         identifies::Symbol = :other, citation = "",
                         references = [], notices = [], requestpage = "",
                         key = "")
        identifies in IDENTIFIES ||
            throw(ArgumentError("`identifies` must be one of $(join(IDENTIFIES, ", ")); got " *
                                "`$identifies`."))
        return new(String(source), String(doi), String(url), String(title),
                   String(publisher), String[string(a) for a in authors],
                   Dict{String, String}(string(k) => string(v)
                                        for (k, v) in identifiers),
                   String(releasedate), String(version), String(description),
                   String(licence), identifies, String(citation),
                   _references(references), _notices(notices),
                   String(requestpage), String(key))
    end
end

# The facts with the named fields replaced, through the keyword constructor, which checks them.
function SourceFacts(facts::SourceFacts; kw...)
    fields = NamedTuple{fieldnames(SourceFacts)}(Tuple(getfield(facts, f)
                                                       for f in fieldnames(SourceFacts)))
    return SourceFacts(; merge(fields, NamedTuple(kw))...)
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
    isempty(f.doi) || print(io, ", doi ", f.doi)
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", f::SourceFacts)
    println(io, "SourceFacts ", repr(f.source), " (", f.identifies, ")")
    for name in (:title, :doi, :url, :publisher, :releasedate, :version,
        :licence, :citation, :key)
        value = getfield(f, name)
        isempty(value) || println(io, "  ", rpad(string(name), 12), value)
    end
    isempty(f.authors) ||
        println(io, "  ", rpad("authors", 12), join(f.authors, "; "))
    isempty(f.description) ||
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

Return whether a file of the source is the cited item or one of its files, so that the same bytes
can be fetched again under the source's identifier - what a data pipeline registers as primary.
A file fetched as its provider publishes it ([`PublishedByProvider`](@ref)) is, where the
identifier names the data; an extract a service made ([`ExtractedByService`](@ref)) is only where
the service gave it an identifier of its own and the facts are that extract's; a file the user
gave ([`SuppliedByUser`](@ref)) never is, since what was done to it is not known.

# Arguments

  - `facts`: what is known of the source.
  - `route`: the route the file came by.
"""
function isprimary(facts::SourceFacts, ::PublishedByProvider)
    return facts.identifies === :data
end

function isprimary(facts::SourceFacts, route::ExtractedByService)
    return !isnothing(route.identifier) && facts.identifies === :data
end

isprimary(::SourceFacts, ::AbstractFileRoute) = false

"""
    sourceentries(facts::SourceFacts, route::AbstractFileRoute; title = facts.title,
                  code = nothing)

Return what a file's provenance record carries of its source, as a dictionary of the sidecar's
keys, to be merged into the record of the file itself: the `source`, the file's `title`, whether
it is `primary` for the route it came by ([`isprimary`](@ref Provenancer.isprimary)), and the
DOI, licence, version, release date, citation, references, notices and request page where known,
each under the key the sidecar uses. The citation is composed from the authors, date, title,
publisher and DOI where the facts give none of their own and have a DOI and authors to compose
one from.

# Arguments

  - `facts`: what is known of the source.
  - `route`: the route the file came by.
  - `title`: the file's own title, the source's by default; a source of several files gives
    each its own.
  - `code`: which layer of the source the file holds, where it is one.
"""
function sourceentries(facts::SourceFacts, route::AbstractFileRoute;
                       title::AbstractString = facts.title, code = nothing)
    record = Dict{String, Any}("source" => facts.source,
                               "title" => String(title),
                               "primary" => isprimary(facts, route))
    for (key, value) in (("doi", facts.doi), ("licence", facts.licence),
        ("version", facts.version), ("release_date", facts.releasedate),
        ("citation", _citationtext(facts)),
        ("request_page", facts.requestpage))
        isempty(value) || (record[key] = value)
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
    isempty(facts.citation) || return facts.citation
    (isempty(facts.doi) || isempty(facts.authors)) && return ""
    year = length(facts.releasedate) >= 4 ?
           " (" * first(facts.releasedate, 4) * ")" : ""
    title = isempty(facts.title) ? "" :
            " " * facts.title *
            (isempty(facts.version) ? "" : " (version " * facts.version * ")") *
            "."
    publisher = isempty(facts.publisher) ? "" : " " * facts.publisher * "."
    return join(facts.authors, ", ") * year * "." * title * publisher *
           " https://doi.org/" * facts.doi
end
