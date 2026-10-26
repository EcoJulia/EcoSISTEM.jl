# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The record of where a published input came from, and the reader of the records the fetches write.

using Unitful

"""
    InputRecord(; role, dataset, code = nothing, path = nothing, url = nothing, request = nothing,
                region = nothing, job = nothing, fetched = nothing, bytes = nothing,
                sha256 = nothing, doi = "", licence = "", version = "", citation = "",
                title = "", primary = nothing, releasedate = "", references = [],
                notices = [], requestpage = "")

The provenance of one published input to an experiment: a file this package fetched, or a dataset
you handed it that it never fetched - an occurrence download, a trait table, a phylogeny. One shape
for both, so that everything an ecosystem was built from can be listed together and a paper's
references written from the list. [`provenance`](@ref) reads the records the package wrote; you
write one for anything it could not have known about.

# Fields

  - `role`: which part of the system the input entered by - `:habitat` for a layer, `:region` for a
    shape or a named region's outline, `:species`, `:phylogeny`, `:dispersal`, `:demography` or
    `:intervention` for what you attached, `:abundance` for records a starting population was
    placed from, `:state` for a saved run a forward run starts from, `:software` for the package
    itself.
  - `dataset`: the source's name - a catalogue key such as `TwentyCR` or `WorldClim{BioClim}`, or
    your own name for a dataset of yours.
  - `code`: which layer of the source, where the input is one.
  - `path`: the file's own name, or a path relative to wherever its record sits - never an absolute
    one, since a record is meant to be committed and shared.
  - `url`: where it was fetched from, where it was.
  - `request`: the request body sent to a service, as it was sent.
  - `region`: the named region a request's area was taken from, as the
    [`NaturalEarthSpec`](@ref) naming it, where it was given as one.
  - `job`: the identifier a service gave the request.
  - `fetched`: when it was fetched, in UTC.
  - `bytes`, `sha256`: the file's size and checksum at fetch, which is what tells a truncated or
    replaced copy from the original.
  - `doi`, `licence`, `version`, `citation`: the dataset's own, blank where unknown; `citation` is
    the text a paper prints.
  - `title`: the file's own title, the release's with the layer and whatever else tells the file
    apart from its siblings, blank where unknown.
  - `primary`: whether the file is the provider's own as published (`true`), or derived from it, as
    an extract a service cut to a request is (`false`); `nothing` where unknown.
  - `releasedate`: when the release was published, as ISO 8601 text, blank where unknown.
  - `references`: the works the dataset's provider asks to be cited, each its `doi` and `text`,
    and `notices`: the statements it asks to accompany the data, each its `strength` (`:required`,
    `:requested` or `:optional`) and `text`. `DD-MMM-YYYY` in a text stands for the date the file
    was fetched and `YYYY` for its year, filled in wherever the record is shown or written.
  - `requestpage`: the address of the documentation's section quoting what the provider asks.
"""
struct InputRecord
    role::Symbol
    dataset::String
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

    function InputRecord(; role::Symbol, dataset::AbstractString,
                         code = nothing, path = nothing, url = nothing,
                         request = nothing, region = nothing, job = nothing,
                         fetched = nothing, bytes = nothing, sha256 = nothing,
                         doi = "", licence = "", version = "", citation = "",
                         title = "", primary::Union{Nothing, Bool} = nothing,
                         releasedate = "", references = [], notices = [],
                         requestpage = "")
        role in _INPUT_ROLES ||
            throw(ArgumentError("`role` must be one of $(join(_INPUT_ROLES, ", ")); got `$role`."))
        isnothing(path) || !isabspath(path) ||
            throw(ArgumentError("an `InputRecord` holds a file's name or a relative path, never " *
                                "an absolute one: it is meant to be committed and shared, and " *
                                "`$path` names a directory on one machine."))
        return new(role, String(dataset), _optstring(code), _optstring(path),
                   _optstring(url),
                   isnothing(request) ? nothing :
                   Dict{String, Any}(string(k) => v for (k, v) in request),
                   _optstring(job), _optdate(fetched),
                   isnothing(bytes) ? nothing : Int(bytes), _optstring(sha256),
                   String(doi), String(licence), String(version),
                   String(citation), _optstring(region), String(title), primary,
                   String(releasedate), _references(references),
                   _notices(notices), String(requestpage))
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
    print(io, "InputRecord(", r.role, ", ", repr(r.dataset))
    isnothing(r.code) || print(io, ", ", repr(r.code))
    isnothing(r.path) || print(io, ", ", repr(r.path))
    isempty(r.doi) || print(io, ", doi = ", repr(r.doi))
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", r::InputRecord)
    println(io, "InputRecord (", r.role, ")")
    for (label, value) in (("dataset", r.dataset), ("title", r.title),
        ("code", r.code), ("file", r.path), ("url", r.url),
        ("region", r.region), ("job", r.job), ("fetched", r.fetched),
        ("bytes", r.bytes), ("sha256", r.sha256), ("primary", r.primary),
        ("doi", r.doi), ("licence", r.licence), ("version", r.version),
        ("released", r.releasedate),
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

"""
    Provenance

What a study area, a habitat or an ecosystem can say about how it came to be, returned by
[`provenance`](@ref): the software, the published inputs, the grid and, for an ecosystem, the run.

# Fields

  - `software`: a named tuple of the `package`, its `version`, the `doi` every version of it is
    cited by, and the `julia` version.
  - `inputs`: an [`InputRecord`](@ref) per input, one per file, in a fixed order.
  - `grid`: a named tuple of the grid's `crs` (`nothing` for a synthetic grid), its `cellsize`, the
    `extent` its cells cover, and how many `cells` it has and how many of them are `active`.
  - `run`: for an ecosystem, a named tuple of its `seed`, its `epoch` (`nothing` without one), the
    run `calendar` its dates were placed by and the time `elapsed`; `nothing` for a study area or a
    habitat.
"""
struct Provenance
    software::NamedTuple
    inputs::Vector{InputRecord}
    grid::NamedTuple
    run::Union{Nothing, NamedTuple}
end

function Base.show(io::IO, p::Provenance)
    n = length(p.inputs)
    return print(io, "Provenance(", p.software.package, " ", p.software.version,
                 ", $(n) input$(n == 1 ? "" : "s"), $(p.grid.active) of $(p.grid.cells) cells",
                 isnothing(p.run) ? "" : ", seed $(p.run.seed)", ")")
end

# The software, grid and run on a line each, then the inputs grouped by dataset: each dataset's DOI
# and citation once, above the files it supplied.
function Base.show(io::IO, ::MIME"text/plain", p::Provenance)
    s = p.software
    println(io, "Provenance")
    println(io, "  software  ", s.package, " ", s.version,
            isempty(s.doi) ? "" : " (doi $(s.doi))", ", Julia ", s.julia)
    println(io,
            "  grid      $(p.grid.active) of $(p.grid.cells) cells active, ",
            "cells of $(p.grid.cellsize)",
            isnothing(p.grid.crs) ? ", synthetic" :
            ", crs $(_crsname(p.grid.crs))")
    isnothing(p.run) ||
        println(io, "  run       seed $(p.run.seed), ",
                _elapsedname(p.run.elapsed), " elapsed",
                isnothing(p.run.epoch) ? "" : " from $(p.run.epoch)",
                _calendarnote(p.run.calendar))
    if isempty(p.inputs)
        print(io, "  inputs    none")
        return nothing
    end
    println(io, "  inputs")
    for dataset in unique(r.dataset for r in p.inputs)
        group = filter(r -> r.dataset == dataset, p.inputs)
        dois = unique(r.doi for r in group if !isempty(r.doi))
        println(io, "    ", dataset,
                isempty(dois) ? "" : " (doi $(join(dois, ", ")))")
        for line in unique(Iterators.flatten(_requestlines(r) for r in group))
            println(io, "      ", line)
        end
        foreach(r -> println(io, "      ", _recordline(r)), group)
    end
    print(io, "  modified  ", _modification(p))
    return nothing
end

# == Functions ======================================================================================

"""
    provenance(path::AbstractString)
    provenance(spec::RasterSpec)
    provenance(spec::ConstructedRasterSpec)
    provenance(spec::AbstractShapeSpec)

Return what is known about where an input came from.

For a file, the [`InputRecord`](@ref) its provenance sidecar holds - written beside every file this
package fetched or first used, as `<file>.provenance.toml` - or `nothing` for a file with none. A
sidecar of that name that this package did not write is left alone and reports `nothing` too, with
a warning saying so, since its fields mean whatever its writer meant.

For a raster spec, one entry per file the spec reads, in order: the file's record, or `nothing`
where the file is absent or has no record. Nothing is fetched to answer; use
[`fetchfiles`](@ref EcoSISTEM.fetchfiles) first for a spec whose files are not there yet. A
combination of layers lists its data members' entries in order, a synthetic member having none.

For a shape spec, the same for the files it outlines: a vector file's record, or a named region's
Natural Earth zip, with the source's licence and citation and the version the zip states filled in
where the record lacks them. A combination of shapes lists its members' in order.

# Arguments

  - `path`: the file, not its sidecar.
  - `spec`: the spec whose files to report on.
"""
function provenance(path::AbstractString)
    isfile(_sidecarpath(path)) || return nothing
    return _readsidecar(_sidecarpath(path))
end

# The parts of the system a published input can enter by.
const _INPUT_ROLES = (:habitat, :region, :species, :phylogeny, :dispersal,
                      :demography, :intervention, :abundance, :state,
                      :software)

# A `provenance` keyword's value as a vector of records: one record, or a vector of them. Anything
# else is refused by name - the untyped fallback exists to give a better message than a
# `MethodError` naming this function.
_recordvector(record::InputRecord) = [record]

_recordvector(records::AbstractVector{InputRecord}) = collect(records)

function _recordvector(x)
    return throw(ArgumentError("`provenance` takes an `InputRecord` or a vector of them; got a " *
                               "$(typeof(x))."))
end

# An optional text field, `nothing` kept, anything else as a `String`.
_optstring(x) = isnothing(x) ? nothing : string(x)

# An optional time field: a `DateTime` as it is, ISO 8601 text with or without its `Z` parsed.
_optdate(::Nothing) = nothing

_optdate(t::Dates.DateTime) = t

_optdate(s::AbstractString) = Dates.DateTime(chopsuffix(String(s), "Z"))

# Whether a parsed sidecar is one this package wrote: its `writer` says so.
function _ourrecord(t::AbstractDict)
    return startswith(string(get(t, "writer", "")), "EcoSISTEM")
end

# A sidecar's TOML as an `InputRecord`. A sidecar another program wrote under the same name is not
# read: its `role` may be prose and its `source` a sentence, so a record built from it would be
# wrong in ways no reader could see.
function _readsidecar(sidecar::AbstractString)
    t = TOML.parsefile(sidecar)
    _ourrecord(t) || begin
        @warn "`$(basename(sidecar))` was not written by EcoSISTEM (writer: "*
              "$(repr(get(t, "writer", nothing)))), so it is left alone and reports nothing; "*
              "delete or rename it to have the file recorded on its next use." maxlog=1 _id=Symbol(sidecar)
        return nothing
    end
    return _recordfromtoml(t)
end

# The record of ours beside `path`, or `nothing` - silently - where there is no sidecar, or one
# another program wrote, or one that does not parse. For a read recording its inputs, where a
# warning on every build would be noise and a malformed file of someone else's must not stop it.
function _ourrecordat(path::AbstractString)
    sidecar = _sidecarpath(path)
    isfile(sidecar) || return nothing
    t = TOML.tryparsefile(sidecar)
    (t isa AbstractDict && _ourrecord(t)) || return nothing
    return _recordfromtoml(t)
end

# A parsed sidecar of ours as an `InputRecord`: the keys the writers use, any it does not know left
# behind.
function _recordfromtoml(t::AbstractDict)
    field(k) = get(t, k, nothing)
    return InputRecord(role = Symbol(something(field("role"), "habitat")),
                       dataset = something(field("source"), field("dataset"),
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
                       requestpage = something(field("request_page"), ""))
end

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
    strength in _NOTICE_STRENGTHS ||
        throw(ArgumentError("a notice's `strength` must be one of " *
                            "$(join(_NOTICE_STRENGTHS, ", ")); got `$strength`."))
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

# What was done to the inputs to build on them, for the statement a modified dataset's licence or
# provider asks for: each was read for this grid and resampled or rasterised onto its cells.
function _modification(p::Provenance)
    crs = isnothing(p.grid.crs) ? "" : " in $(_crsname(p.grid.crs))"
    return "contains modified data: each input was read for this grid's extent and resampled " *
           "or rasterised onto its cells of $(p.grid.cellsize)$(crs)."
end

# One input's line under its dataset in a provenance listing: its role, its layer where it has one,
# and its file, or where it was fetched from where it has no file.
function _recordline(r::InputRecord)
    parts = (string(r.role), something(r.code, ""),
             something(r.path, r.url, ""))
    return join(filter(!isempty, parts), "  ")
end

# A run's elapsed time in the unit it is most easily read in: seconds under a day, days under a
# year, and years beyond.
function _elapsedname(t::Unitful.Time)
    t < 1.0u"d" && return string(round(typeof(1.0u"s"), t, digits = 1))
    t < 1.0u"yr" && return string(round(typeof(1.0u"d"), t, digits = 1))
    return string(round(typeof(1.0u"yr"), t, digits = 2))
end
