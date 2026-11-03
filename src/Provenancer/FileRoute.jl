# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The routes a file comes by, and the trait a thing that fetches or finds a file implements.

"""
    AbstractFileRoute

The route a file came by - how its bytes came to be - which decides what its record says it is
and how a registry is told of it: [`VersionOfRecord`](@ref), [`GeneratedOnRequest`](@ref) or [`OtherOrigin`](@ref). Every
type implementing [`IsFetchable`](@ref Provenancer.IsFetchable) says which by a
[`route`](@ref Provenancer.route) method.
"""
abstract type AbstractFileRoute end

"""
    VersionOfRecord()

The route of a file that is its source's version of record: its bytes as the source issued them,
however they reached you - fetched from a URL of the provider's or by the source's own resolution,
or sent by its author. A file of a source whose identifier names the data is then the source's
own, the cited item or one of its files.
"""
struct VersionOfRecord <: AbstractFileRoute end

"""
    GeneratedOnRequest(identifier = nothing)

The route of a file a service generated to a request, whose bytes did not exist until they were
asked for - a box of months cut from a reanalysis, an occurrence download built to a predicate.
The request identifies what was asked for, and the
service's own identifier for the result, where it gives one - the DOI an occurrence download is
minted - identifies the file, and makes it the cited item itself; `nothing` says the service gave
none, as the Climate Data Store does not.

# Fields

  - `identifier`: the service's identifier for the file, or `nothing`.
"""
struct GeneratedOnRequest <: AbstractFileRoute
    identifier::Union{Nothing, String}

    function GeneratedOnRequest(identifier::Union{Nothing, AbstractString} = nothing)
        return new(isnothing(identifier) ? nothing : String(identifier))
    end
end

function Base.show(io::IO, route::GeneratedOnRequest)
    return print(io, "GeneratedOnRequest(",
                 isnothing(route.identifier) ? "" : repr(route.identifier), ")")
end

"""
    OtherOrigin()

The route of a file of any other origin: one named only by its path, one made or altered by you,
or one recorded before routes were. Nothing is claimed of how its bytes relate to its source's, so
it is not taken as the cited item ([`isprimary`](@ref Provenancer.isprimary)).
"""
struct OtherOrigin <: AbstractFileRoute end

@traitdef IsFetchable{X}

@doc """
    IsFetchable{X}

The trait of a type `X` that fetches or finds a file and can say which route it comes by: an entry
describing one file - a download from a URL, a request to a service - or a dataset type whose source
resolves its own files, into another package's cache. A type implements it with
`@traitimpl IsFetchable{X}` and defines [`route`](@ref Provenancer.route) for it - on the value for
an entry, on the type itself for a dataset - and whatever the package using this machinery asks of
its fetchable things.
""" IsFetchable

"""
    route(x)
    route(::Type{S})

Return the [`AbstractFileRoute`](@ref Provenancer.AbstractFileRoute) the files described by `x`
come by. Defined by every type implementing [`IsFetchable`](@ref Provenancer.IsFetchable): on the
value for an entry describing one file, on the type itself for a dataset type whose source resolves
its own files.

# Arguments

  - `x`: the entry describing the file, or the dataset type.
"""
function route end

# A route as a record's sidecar keeps it: its name under `route`, and under `route_identifier` the
# identifier a service gave an extract, where it gave one.
function _routeentries(route::AbstractFileRoute)
    entries = Dict{String, Any}("route" => string(nameof(typeof(route))))
    if route isa GeneratedOnRequest && !isnothing(route.identifier)
        entries["route_identifier"] = route.identifier
    end
    return entries
end

_routeentries(::Nothing) = Dict{String, Any}()

# The route a sidecar's `route` and `route_identifier` keys name, or `nothing` for a record written
# before routes were recorded, or one naming a route this version does not know.
function _routefrom(name, identifier)
    isnothing(name) && return nothing
    name == "VersionOfRecord" && return VersionOfRecord()
    name == "GeneratedOnRequest" && return GeneratedOnRequest(identifier)
    name == "OtherOrigin" && return OtherOrigin()
    return nothing
end
