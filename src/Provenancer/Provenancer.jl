# SPDX-License-Identifier: LGPL-3.0-or-later

"""
    Provenancer

The provenance machinery EcoSISTEM runs on, written to depend on nothing of EcoSISTEM's so that
any package can use it: a record of where a file came from ([`ProvenanceRecord`](@ref)), written
beside the file as `<file>.provenance.toml` by every fetch, or by [`record_provenance`](@ref) for a
file the user has, and read back by [`provenance`](@ref); what is known of a source
([`SourceFacts`](@ref)), from a catalogue, from the metadata its DOI resolves to
([`lookup_doi`](@ref)) or by hand; a download that resumes, locks and records itself
([`URLFile`](@ref), [`fetchfile`](@ref Provenancer.fetchfile)); the names a file and its record
take in a FAIR data pipeline run ([`Track`](@ref)); and the account of what a result was built
from ([`Provenance`](@ref)): the software, its inputs and whatever context the package using it
adds.

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
using JSON: JSON
using SHA: SHA
using Scratch: get_scratch!
using SimpleTraits
using TOML: TOML

# The types first, each with its constructors and `show`: the routes a file comes by and the trait
# a fetchable thing implements; the record of one file, which keeps its route, and the reader of
# the sidecars holding them; the account of a result and the context a package adds to it; a
# pipeline run as the route files take. What is known of a source, and a file fetched from a URL
# with it. Then the functions: the asset cache, the download and the sidecar, what a tracked file
# is registered as, and what a DOI resolves to.
include("FileRoute.jl")
include("ProvenanceRecord.jl")
include("Provenance.jl")
include("Track.jl")
include("SourceFacts.jl")
include("URLFile.jl")
include("assets.jl")
include("naming.jl")
include("doi.jl")

export ProvenanceRecord, Provenance, Track, NoContext, SoftwareRecord,
       PublishedByProvider, ExtractedByService, SuppliedByUser, SourceFacts,
       URLFile, provenance, write_provenance, lookup_doi, record_provenance

public AbstractProvenanceContext, AbstractFileRoute, IsFetchable, route,
       NOTICE_STRENGTHS, IDENTIFIES, isprimary, sourceentries, contextsummary,
       contextlines, contexttables, modification, softwarerecord, tomltable,
       tomlvalue

public assetdir, set_assetdir!, fetchlocked, download!, fetchfile, sidecarpath,
       recordpresent, writesidecar, verifyfile, recordat, baserecord,
       downloadrecord, iscatalogued

end
