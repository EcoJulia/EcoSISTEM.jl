# SPDX-License-Identifier: LGPL-3.0-or-later

"""
    Provenancer

The provenance machinery EcoSISTEM runs on, written to depend on nothing of EcoSISTEM's so that
any package can use it: a record of where a file came from ([`InputRecord`](@ref)), written beside
the file as `<file>.provenance.toml` and read back by [`provenance`](@ref); a download that
resumes, locks and records itself ([`download!`](@ref Provenancer.download!),
[`writesidecar`](@ref Provenancer.writesidecar)); the names a file and its record take in a FAIR data
pipeline run ([`Track`](@ref)); and the account of what a result was built from
([`Provenance`](@ref)): the software, its inputs and whatever context the package using it adds.

A package using it implements the [`IsFetchable`](@ref Provenancer.IsFetchable) trait for each
kind of thing it fetches, says which [`AbstractFileRoute`](@ref Provenancer.AbstractFileRoute) the
file comes by, and subtypes [`AbstractProvenanceContext`](@ref Provenancer.AbstractProvenanceContext)
for what it adds to a result's account. Nothing here reads a configuration of the package using it:
the cache directory and the software record are asked for by module.
"""
module Provenancer

using Dates: Dates
using Downloads: Downloads
using FileWatching: mkpidlock
using Preferences: Preferences
using SHA: SHA
using Scratch: get_scratch!
using SimpleTraits
using TOML: TOML

# The types first, each with its constructors and `show`: the record of one input and the reader of
# the sidecars holding them; the account of a result and the context a package adds to it; a
# pipeline run as the route files take; the routes a file comes by and the trait a fetchable thing
# implements. Then the functions: the asset cache, the download and the sidecar, and what a tracked
# file is registered as.
include("InputRecord.jl")
include("Provenance.jl")
include("Track.jl")
include("FileRoute.jl")
include("assets.jl")
include("naming.jl")

export InputRecord, Provenance, Track, NoContext, SoftwareRecord,
       PublishedByProvider, ExtractedByService, SuppliedByUser,
       provenance, write_provenance

public AbstractProvenanceContext, AbstractFileRoute, IsFetchable, route,
       NOTICE_STRENGTHS, contextsummary, contextlines, contexttables,
       modification, softwarerecord, tomltable, tomlvalue

public assetdir, set_assetdir!, fetchlocked, download!, sidecarpath,
       recordpresent, writesidecar, verifyfile, recordat, baserecord,
       downloadrecord

end
