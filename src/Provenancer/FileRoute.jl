# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The routes a file comes by, and the trait a thing that fetches or finds a file implements.

"""
    AbstractFileRoute

The route a file came by, which decides what its record says it is and how a registry is told of
it: [`PublishedByProvider`](@ref), [`ExtractedByService`](@ref) or [`SuppliedByUser`](@ref). Every
type implementing [`IsFetchable`](@ref Provenancer.IsFetchable) says which by a
[`route`](@ref Provenancer.route) method.
"""
abstract type AbstractFileRoute end

"""
    PublishedByProvider()

The route of a file fetched as its provider publishes it, from a URL of the provider's or by its
source's own resolution. A file of a source whose identifier names the data is then the source's
own, the cited item or one of its files, so that the same bytes can be fetched again.
"""
struct PublishedByProvider <: AbstractFileRoute end

"""
    ExtractedByService(identifier = nothing)

The route of a file a service made to a request - a box of months cut from a reanalysis, an
occurrence download built to a predicate. The request identifies what was asked for, and the
service's own identifier for the result, where it gives one - the DOI an occurrence download is
minted - identifies the file, and makes it the cited item itself; `nothing` says the service gave
none, as the Climate Data Store does not.

# Fields

  - `identifier`: the service's identifier for the file, or `nothing`.
"""
struct ExtractedByService <: AbstractFileRoute
    identifier::Union{Nothing, String}

    function ExtractedByService(identifier::Union{Nothing, AbstractString} = nothing)
        return new(isnothing(identifier) ? nothing : String(identifier))
    end
end

function Base.show(io::IO, route::ExtractedByService)
    return print(io, "ExtractedByService(",
                 isnothing(route.identifier) ? "" : repr(route.identifier), ")")
end

"""
    SuppliedByUser()

The route of a file the user named, which the package never fetched and knows only as the user
describes it.
"""
struct SuppliedByUser <: AbstractFileRoute end

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
    if route isa ExtractedByService && !isnothing(route.identifier)
        entries["route_identifier"] = route.identifier
    end
    return entries
end

_routeentries(::Nothing) = Dict{String, Any}()

# The route a sidecar's `route` and `route_identifier` keys name, or `nothing` for a record written
# before routes were recorded, or one naming a route this version does not know.
function _routefrom(name, identifier)
    isnothing(name) && return nothing
    name == "PublishedByProvider" && return PublishedByProvider()
    name == "ExtractedByService" && return ExtractedByService(identifier)
    name == "SuppliedByUser" && return SuppliedByUser()
    return nothing
end
