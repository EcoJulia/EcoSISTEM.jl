# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The record of where a published input came from, and the reader of the records the fetches write.

using Unitful

"""
    InputRecord(; role, dataset, code = nothing, path = nothing, url = nothing, request = nothing,
                job = nothing, fetched = nothing, bytes = nothing, sha256 = nothing, doi = "",
                licence = "", version = "", citation = "")

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
  - `job`: the identifier a service gave the request.
  - `fetched`: when it was fetched, in UTC.
  - `bytes`, `sha256`: the file's size and checksum at fetch, which is what tells a truncated or
    replaced copy from the original.
  - `doi`, `licence`, `version`, `citation`: the dataset's own, blank where unknown; `citation` is
    the text a paper prints.
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

    function InputRecord(; role::Symbol, dataset::AbstractString,
                         code = nothing, path = nothing, url = nothing,
                         request = nothing, job = nothing, fetched = nothing,
                         bytes = nothing, sha256 = nothing, doi = "",
                         licence = "", version = "", citation = "")
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
                   String(citation))
    end
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
    for (label, value) in (("dataset", r.dataset), ("code", r.code),
        ("file", r.path), ("url", r.url), ("job", r.job),
        ("fetched", r.fetched), ("bytes", r.bytes),
        ("sha256", r.sha256), ("doi", r.doi),
        ("licence", r.licence), ("version", r.version),
        ("citation", r.citation))
        (isnothing(value) || (value isa AbstractString && isempty(value))) &&
            continue
        println(io, "  ", rpad(label, 9), value)
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
  - `run`: for an ecosystem, a named tuple of its `seed`, its `epoch` (`nothing` without one) and
    the time `elapsed`; `nothing` for a study area or a habitat.
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
                isnothing(p.run.epoch) ? "" : " from $(p.run.epoch)")
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
        for citation in unique(r.citation
                               for r in group
                               if !isempty(r.citation))
            println(io, "      cite  ", citation)
        end
        foreach(r -> println(io, "      ", _recordline(r)), group)
    end
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
                       citation = something(field("citation"), ""))
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
