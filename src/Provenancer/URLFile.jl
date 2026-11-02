# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A file fetched from a URL on first use: what a package with no catalogue of its own describes a
# download by, and records with whatever it knows of the source.

"""
    URLFile(url::AbstractString, path::AbstractString; facts = nothing)

A file downloaded from `url` to `path` on first use by [`fetchfile`](@ref Provenancer.fetchfile),
which writes its provenance record beside it: what the download itself says, and what `facts` say
of the source. Nothing is fetched at construction. It implements
[`IsFetchable`](@ref Provenancer.IsFetchable), and its files come by [`PublishedByProvider`](@ref),
as the provider serves them.

# Fields

  - `url`: where to fetch the file from.
  - `path`: where it lives once fetched.
  - `facts`: what is known of its source, as [`SourceFacts`](@ref), or `nothing`, which records
    the download alone.
"""
struct URLFile
    url::String
    path::String
    facts::Union{Nothing, SourceFacts}

    function URLFile(url::AbstractString, path::AbstractString;
                     facts::Union{Nothing, SourceFacts} = nothing)
        return new(String(url), String(path), facts)
    end
end

@traitimpl IsFetchable{URLFile}

route(::URLFile) = PublishedByProvider()

function Base.show(io::IO, f::URLFile)
    print(io, "URLFile(", repr(f.url), ", ", repr(f.path))
    isnothing(f.facts) || print(io, ", facts = ", f.facts)
    return print(io, ")")
end
