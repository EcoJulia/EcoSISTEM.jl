# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The account of what a result was built from - the software, its inputs and the context the package
# using this machinery adds - and the TOML file it is written as.

"""
    AbstractProvenanceContext

What a package adds to the account of a result beyond the software and its inputs: for EcoSISTEM
the grid a model sits on and the run, for another package whatever its results depend on. A
subtype is a struct of those things; the generic display and TOML list its fields by name, and a
subtype gives [`contextsummary`](@ref Provenancer.contextsummary),
[`contextlines`](@ref Provenancer.contextlines), [`contexttables`](@ref Provenancer.contexttables)
and [`modification`](@ref Provenancer.modification) methods where that does not suit it.
"""
abstract type AbstractProvenanceContext end

"""
    NoContext()

The context of a result whose account is its software and inputs alone.
"""
struct NoContext <: AbstractProvenanceContext end

"""
    SoftwareRecord(package, version, doi, julia)

The software a result came from, as its account records it. [`softwarerecord`](@ref
Provenancer.softwarerecord) builds one for a package.

# Fields

  - `package`: the package's name.
  - `version`: its version, as text.
  - `doi`: the DOI every version of it is cited by, or `""` where it has none.
  - `julia`: the Julia version it ran on, as text.
"""
struct SoftwareRecord
    package::String
    version::String
    doi::String
    julia::String
end

"""
    Provenance(software::SoftwareRecord, inputs, context = NoContext())

What a result can say about how it came to be: the software, the inputs it was built from and the
context the package using this machinery adds - for EcoSISTEM's study areas, habitats and
ecosystems the grid and, for an ecosystem, the run. Shown as a listing and written as TOML by
[`write_provenance`](@ref).

# Fields

  - `software`: a [`SoftwareRecord`](@ref).
  - `inputs`: a [`ProvenanceRecord`](@ref) per input, one per file, in a fixed order.
  - `context`: an [`AbstractProvenanceContext`](@ref Provenancer.AbstractProvenanceContext), the
    package's own; `NoContext()` where there is nothing to add.
"""
struct Provenance
    software::SoftwareRecord
    inputs::Vector{ProvenanceRecord}
    # Abstract by design: the context is whatever the package using this machinery defines, and
    # nothing reads it on a hot path - a `Provenance` is assembled once per recording - so the
    # field is left abstract rather than making the struct parametric, which would turn every
    # `Provenance` in a signature or a field into a `UnionAll`.
    context::AbstractProvenanceContext

    function Provenance(software::SoftwareRecord,
                        inputs::AbstractVector{ProvenanceRecord},
                        context::AbstractProvenanceContext = NoContext())
        return new(software, collect(inputs), context)
    end
end

function Base.show(io::IO, p::Provenance)
    n = length(p.inputs)
    return print(io, "Provenance(", p.software.package, " ", p.software.version,
                 ", $(n) input$(n == 1 ? "" : "s")", contextsummary(p.context),
                 ")")
end

# The software and the context on a line each, then the inputs grouped by source: each source's
# DOI and citation once, above the files it supplied; then, where the context says the inputs were
# modified, how.
function Base.show(io::IO, ::MIME"text/plain", p::Provenance)
    s = p.software
    println(io, "Provenance")
    println(io, "  software  ", s.package, " ", s.version,
            isempty(s.doi) ? "" : " (doi $(s.doi))", ", Julia ", s.julia)
    for line in contextlines(p.context)
        println(io, "  ", line)
    end
    if isempty(p.inputs)
        print(io, "  inputs    none")
        return nothing
    end
    println(io, "  inputs")
    for source in unique(r.source for r in p.inputs)
        group = filter(r -> r.source == source, p.inputs)
        dois = unique(r.doi for r in group if !isnothing(r.doi))
        println(io, "    ", source,
                isempty(dois) ? "" : " (doi $(join(dois, ", ")))")
        for line in unique(Iterators.flatten(_requestlines(r) for r in group))
            println(io, "      ", line)
        end
        foreach(r -> println(io, "      ", _recordline(r)), group)
    end
    modified = modification(p.context)
    isempty(modified) || print(io, "  modified  ", modified)
    return nothing
end

# == Functions ======================================================================================

"""
    write_provenance(path::AbstractString, p::Provenance)

Write a [`Provenance`](@ref) to `path` as TOML, and return the path: the software, the context's
tables, and every input with its DOI, licence, version, citation, file name and checksum - the
file a paper's repository commits beside its figures, so a reader can see what made them.

An input is written as its record holds it, named by its file and never by a directory, so nothing
in the file is particular to the machine that wrote it. Fetch times are written as TOML dates, and
anything TOML has no form for as its printed text.

# Arguments

  - `path`: the file to write; its directory is created if it does not exist.
  - `p`: the provenance to write.
"""
function write_provenance(path::AbstractString, p::Provenance)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        return TOML.print(io, _provenancetoml(p), sorted = true)
    end
    return path
end

"""
    contextsummary(context::AbstractProvenanceContext)

Return what a [`Provenance`](@ref)'s one-line display adds for its context, after the count of
inputs: text beginning with a comma, or `""` for nothing, the default.

# Arguments

  - `context`: the result's context.
"""
contextsummary(::AbstractProvenanceContext) = ""

"""
    contextlines(context::AbstractProvenanceContext)

Return the lines a [`Provenance`](@ref)'s listing gives its context, after the software's line and
before the inputs: each a label padded to ten characters and a value. By default one per field of
the context with something to say, the field's name as the label.

# Arguments

  - `context`: the result's context.
"""
function contextlines(context::AbstractProvenanceContext)
    return [rpad(string(f), 10) * string(getfield(context, f))
            for f in fieldnames(typeof(context))
            if !_tomlblank(getfield(context, f))]
end

"""
    contexttables(context::AbstractProvenanceContext)

Return the tables a [`Provenance`](@ref)'s TOML file holds for its context, keyed by name, beside
`software` and `inputs`. By default one per field of the context with something to say, the field's
name as the key and its value as [`tomlvalue`](@ref Provenancer.tomlvalue) writes it.

# Arguments

  - `context`: the result's context.
"""
function contexttables(context::AbstractProvenanceContext)
    return Dict{String, Any}(string(f) => tomlvalue(getfield(context, f))
                             for f in fieldnames(typeof(context))
                             if !_tomlblank(getfield(context, f)))
end

"""
    modification(context::AbstractProvenanceContext)

Return the statement of what was done to the inputs to build on them, for a licence or a provider
that asks for one, or `""` where nothing is to be said, the default. Written as a
[`Provenance`](@ref)'s `modification` and shown as its `modified` line, where it has inputs.

# Arguments

  - `context`: the result's context.
"""
modification(::AbstractProvenanceContext) = ""

"""
    softwarerecord(mod::Module)

Return the [`SoftwareRecord`](@ref) of the package `mod` belongs to: its name, its version, the
concept DOI its `codemeta.json` declares as its `identifier` - the Zenodo DOI every version is
cited by, `""` where the file or the DOI is missing - and the Julia running.

# Arguments

  - `mod`: the package's module, or any module within it.
"""
function softwarerecord(mod::Module)
    root = Base.moduleroot(mod)
    return SoftwareRecord(string(nameof(root)),
                          string(something(pkgversion(root), "")),
                          _conceptdoi(root), string(VERSION))
end

"""
    tomltable(x)

Return the fields of a named tuple, a [`SoftwareRecord`](@ref) or a [`ProvenanceRecord`](@ref) as
a TOML table, leaving out those with nothing to say - `nothing`, empty text or an empty list - each
value as [`tomlvalue`](@ref Provenancer.tomlvalue) writes it. A record's references and notices
are written with their placeholders filled from its fetch time, its route as its sidecar keeps it,
and where the file and its record were when it was read not at all: the table names the file
alone, as the sidecar does.

# Arguments

  - `x`: what to write.
"""
function tomltable(nt::NamedTuple)
    return Dict{String, Any}(string(k) => tomlvalue(v)
                             for (k, v) in pairs(nt) if !_tomlblank(v))
end

function tomltable(s::SoftwareRecord)
    return Dict{String, Any}(string(k) => tomlvalue(getfield(s, k))
                             for k in fieldnames(SoftwareRecord)
                             if !_tomlblank(getfield(s, k)))
end

function tomltable(r::ProvenanceRecord)
    function written(k)
        v = getfield(r, k)
        k in (:references, :notices) || return v
        return [merge(x, (text = _fillplaceholders(x.text, r.fetched),))
                for x in v]
    end
    table = Dict{String, Any}(string(k) => tomlvalue(written(k))
                              for k in fieldnames(ProvenanceRecord)
                              if !(k in (:route, :location, :sidecar)) &&
        !_tomlblank(getfield(r, k)))
    return merge!(table, _routeentries(r.route))
end

"""
    tomlvalue(v)

Return `v` in a form TOML holds: text, booleans, floats, dates and integers that fit TOML's own as
they are, a seed too large for a signed integer as text, dictionaries and named tuples as tables,
vectors as arrays of these, and anything else - a quantity, a coordinate reference system, a
symbol, another calendar's date - as its printed text. A package writing a context adds methods
for what it wants written otherwise.

# Arguments

  - `v`: the value.
"""
tomlvalue(v::AbstractString) = String(v)

tomlvalue(v::Bool) = v

tomlvalue(v::Union{Float64, Dates.Date, Dates.DateTime}) = v

function tomlvalue(v::Integer)
    return typemin(Int64) <= v <= typemax(Int64) ? Int64(v) : string(v)
end

function tomlvalue(v::AbstractDict)
    return Dict{String, Any}(string(k) => tomlvalue(x) for (k, x) in v)
end

tomlvalue(v::AbstractVector) = Any[tomlvalue(x) for x in v]

function tomlvalue(v::NamedTuple)
    return Dict{String, Any}(string(k) => tomlvalue(x) for (k, x) in pairs(v))
end

# Deliberately untyped, as the fallback for everything else.
tomlvalue(v) = string(v)

# The concept DOI the `codemeta.json` of the package `mod` belongs to declares as its `identifier`,
# bare, or `""` where the file is missing or the identifier is no DOI.
function _conceptdoi(mod::Module)
    path = pkgdir(Base.moduleroot(mod), "codemeta.json")
    (isnothing(path) || !isfile(path)) && return ""
    meta = JSON.parse(read(path, String))
    identifier = meta isa AbstractDict ? get(meta, "identifier", nothing) :
                 nothing
    identifier isa AbstractString || return ""
    return something(_trydoi(identifier), "")
end

# A `Provenance` as the tables its TOML file holds: `software`, `inputs` (a table per input), the
# context's tables, and `modification` where the context states one and there are inputs.
function _provenancetoml(p::Provenance)
    toml = Dict{String, Any}("software" => tomltable(p.software),
                             "inputs" => Dict{String, Any}[tomltable(r)
                                               for r in p.inputs])
    merge!(toml, contexttables(p.context))
    modified = modification(p.context)
    (isempty(p.inputs) || isempty(modified)) ||
        (toml["modification"] = modified)
    return toml
end

# Whether a field has nothing to say: `nothing`, empty text, or an empty list.
function _tomlblank(v)
    return isnothing(v) ||
           (v isa Union{AbstractString, AbstractVector} && isempty(v))
end
