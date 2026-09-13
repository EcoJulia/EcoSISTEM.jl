# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The record of where a published input came from, and the reader of the records the fetches write.

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
    `:intervention` for what you attached, `:software` for the package itself.
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

# The parts of the system a published input can enter by.
const _INPUT_ROLES = (:habitat, :region, :species, :phylogeny, :dispersal,
                      :demography, :intervention, :software)

# An optional text field, `nothing` kept, anything else as a `String`.
_optstring(x) = isnothing(x) ? nothing : string(x)

# An optional time field: a `DateTime` as it is, ISO 8601 text with or without its `Z` parsed.
_optdate(::Nothing) = nothing

_optdate(t::Dates.DateTime) = t

_optdate(s::AbstractString) = Dates.DateTime(chopsuffix(String(s), "Z"))

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

# == Functions ======================================================================================

"""
    provenance(path::AbstractString)
    provenance(spec::RasterSpec)

Return what is known about where an input came from.

For a file, the [`InputRecord`](@ref) its provenance sidecar holds - written beside every file this
package fetched, as `<file>.provenance.toml` - or `nothing` for a file with none, which is one it
never fetched.

For a raster spec, one entry per file the spec reads, in order: the file's record, or `nothing`
where the file is absent or has no record. Nothing is fetched to answer; use
[`fetchfiles`](@ref EcoSISTEM.fetchfiles) first for a spec whose files are not there yet.

# Arguments

  - `path`: the file, not its sidecar.
  - `spec`: the spec whose files to report on.
"""
function provenance(path::AbstractString)
    isfile(_sidecarpath(path)) || return nothing
    return _readsidecar(_sidecarpath(path))
end

# A sidecar's TOML as an `InputRecord`: the keys the writers use, any it does not know left behind.
function _readsidecar(sidecar::AbstractString)
    t = TOML.parsefile(sidecar)
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
