# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a tracked spec's files are registered as with a FAIR data pipeline run, the part that knows a
# spec and the catalogue: each file a spec reads and the route it comes by, the source its files
# are registered with - the catalogue row's facts - the facts the catalogue gives a file now, the
# registration keywords built from those facts, and how a spec keeps its track. The data product names, the version
# stamping and the config lines are `Provenancer`'s (`src/Provenancer/naming.jl`); the calls that
# ask a registry and register are the DataPipeline extension's.

# Each file a spec reads, where it is or would be, with the route it comes by - `VersionOfRecord()`
# for a file fetched as its source issues it (by the source's own resolution, or from a URL),
# `GeneratedOnRequest()` for one a service generates to a request, `OtherOrigin()` for a file the
# user names by its path - the layer `code` it holds where one is known, and the `url` it is fetched from where
# it has one.
_specroutes(spec::RasterSpec) = _specroutes(spec, spec.files)

function _specroutes(spec::RasterSpec, ::Nothing)
    route = _fileroute(spec.source)
    files = _sourcefiles(spec.source, spec.code; _readkw(spec)...)
    isnothing(files) &&
        return [(path = p, route = route, code = nothing, url = nothing)
                for p in _specpaths(spec)]
    return [(path = f.path, route = route, code = f.code, url = f.url)
            for f in files]
end

function _specroutes(spec::RasterSpec, files::AbstractVector)
    code = spec.code isa AbstractVector ? nothing : spec.code
    return [(path = _localpath(f), route = _fileroute(f), code = code,
             url = _entryurl(f))
            for f in files]
end

function _specroutes(spec::ShapeSpec)
    return [(path = _localpath(spec.path), route = _fileroute(spec.path),
             code = nothing, url = _entryurl(spec.path))]
end

function _specroutes(spec::NaturalEarthSpec)
    url = _nesource(_checklevel(spec.level)).url
    return [(path = p, route = Provenancer.route(spec), code = nothing,
             url = url)
            for p in _specpaths(spec)]
end

function _specroutes(authority::WFOPlantList)
    asset = _wfoasset(authority)
    return [(path = _localpath(asset), route = Provenancer.route(authority),
             code = nothing, url = asset.url)]
end

# The route a file comes by, asked of whatever describes it: a path the user gave; an entry of a
# spec's `files`, a `CachedAsset` or a `CDSRequest`, which says of itself; or the dataset type of a
# source resolving its own files. Anything else is refused by name: implementing `IsFetchable` is
# the way in.
_fileroute(::AbstractString) = OtherOrigin()

# A file given by its record comes by the route the record keeps, or as one the user named where
# the record predates routes.
_fileroute(record::ProvenanceRecord) = something(record.route, OtherOrigin())

@traitfn function _fileroute(entry::E) where {E; IsFetchable{E}}
    return Provenancer.route(entry)
end

@traitfn _fileroute(entry::E) where {E; !IsFetchable{E}} = _notfetchable(E)

@traitfn function _fileroute(::Type{S}) where {S; IsFetchable{S}}
    return Provenancer.route(S)
end

@traitfn _fileroute(::Type{S}) where {S; !IsFetchable{S}} = _notfetchable(S)

function _notfetchable(T)
    return throw(ArgumentError("`$T` does not say which route its files come by: a `files` " *
                               "entry is a path or a value of a type marked " *
                               "`@traitimpl Provenancer.IsFetchable{$T}` with a " *
                               "`Provenancer.route` method, and a source resolving its own " *
                               "files is such a type."))
end

_entryurl(asset::CachedAsset) = asset.url

_entryurl(record::ProvenanceRecord) = record.url

_entryurl(::Any) = nothing

# The source a tracked spec's files are registered with - its facts, and whether its files state
# their own versions rather than the source giving one (Natural Earth's zips, each naming its
# own) - from the catalogue row of the spec's source. A spec of a file no row describes is
# refused: such a file enters a spec through its record, which `Provenancer.record_provenance`
# writes and, given the run's track, registers.
function _trackedsource(spec)
    rec = _catalogued(spec)
    isnothing(rec) &&
        error("`$(basename(first(_specpaths(spec))))`, read by $spec, is of no catalogued " *
              "source, so what it is cannot be registered: a file of your own enters a spec " *
              "through its record, `Provenancer.record_provenance(path, doi or facts, " *
              "role = ..., track = ...)`, which registers it with the source the facts name.")
    return (facts = SourceFacts(rec), fileversioned = _fileversioned(rec))
end

# The catalogue row of the source a spec reads, or `nothing` for a spec of a file no row
# describes.
_catalogued(spec::RasterSpec) = _datasetrecord(spec.source)

_catalogued(::NaturalEarthSpec) = _datasetrecord(NaturalEarthLevel)

_catalogued(::ShapeSpec) = nothing

_catalogued(::WFOPlantList) = _datasetrecord(WFOPlantList)

# What the catalogue says now of one file of a tracked spec, as its record holds it: the catalogue's
# facts for the route it comes by, where it is fetched from, the release a track gives a file the
# user gave, and the version a file states of itself. `entry` is one of `_specroutes(spec)`.
function _currentfacts(spec::RasterSpec, entry)
    layer = isnothing(entry.code) ? nothing :
            _layerornothing(spec.source, entry.code)
    facts = _cataloguerecord(spec.source, layer, entry.path,
                             route = entry.route)
    isnothing(entry.url) || (facts["url"] = entry.url)
    return _trackrelease!(facts, spec.track, entry.route)
end

# The release a track names for the files the user gives, written over a record's version and
# release date, which only such a file's track may say.
function _trackrelease!(facts::AbstractDict, track, route::AbstractFileRoute)
    (route isa OtherOrigin && !isnothing(track)) || return facts
    isnothing(track.release_version) ||
        (facts["version"] = track.release_version)
    isnothing(track.release_date) ||
        (facts["release_date"] = track.release_date)
    return facts
end

function _currentfacts(spec::NaturalEarthSpec, entry)
    facts = _cataloguerecord(NaturalEarthLevel, nothing, entry.path,
                             route = entry.route)
    facts["url"] = entry.url
    version = _fileversion(datasetinfo(NaturalEarthLevel), entry.path)
    isnothing(version) || (facts["version"] = version)
    return facts
end

function _currentfacts(::WFOPlantList, entry)
    facts = _cataloguerecord(WFOPlantList, nothing, entry.path,
                             route = entry.route)
    facts["url"] = entry.url
    return facts
end

# A spec's files fetched by its own route, as an untracked spec's are, for a tracked fetch run to
# register after.
_fetchspec(spec::RasterSpec) = _fetchedpaths(spec, spec.files)

function _fetchspec(spec::NaturalEarthSpec)
    return [assetpath(_nesource(_checklevel(spec.level)))]
end

_fetchspec(spec::ShapeSpec) = [String(_resolvepath(spec.path))]

_fetchspec(authority::WFOPlantList) = [assetpath(_wfoasset(authority))]

# The provenance record of a file the user gave a spec, written beside it where it has none: what
# the catalogue says of it, the release its track names, and no fetch time, which is not known.
# `entry` is the file's, one of `_specroutes(spec)`.
function _recordgiven(spec::RasterSpec, entry)
    isfile(entry.path) ||
        error("`$(entry.path)`, a file given to $spec, is not there, so it cannot be registered.")
    Provenancer.recordpresent(entry.path) do
        return merge!(Provenancer.baserecord(entry.path, _givenrole(spec)),
                      _currentfacts(spec, entry))
    end
    return entry.path
end

# The role a file the user gave a spec is recorded under: a layer's.
_givenrole(::RasterSpec) = :habitat

# A source's layer by its code, or `nothing` where the catalogue has none of that code.
function _layerornothing(source::Type, code)
    return try
        layerinfo(source, code)
    catch
        nothing
    end
end

# Whether a source's files state their own versions rather than its row giving one: Natural
# Earth's zips, each naming its own.
_fileversioned(rec::DatasetRecord) = rec.format === :Shapefile

# The release date a file is registered with: a stamped version's day; else its record's, the
# row's or the one its track gave; else, for a file stating its own version (`fileversioned`),
# when its provider last changed it, or failing that when it was fetched; else `nothing`, the
# citation's date standing.
function _registereddate(r::ProvenanceRecord, fileversioned::Bool,
                         stamp::Union{Nothing, Dates.Date})
    isnothing(stamp) || return string(stamp)
    isnothing(r.releasedate) || return r.releasedate
    if fileversioned
        day = something(r.modified, r.fetched, Some(nothing))
        isnothing(day) || return string(Dates.Date(day))
    end
    return nothing
end

function _registereddate(r::ProvenanceRecord, rec::DatasetRecord,
                         stamp::Union{Nothing, Dates.Date})
    return _registereddate(r, _fileversioned(rec), stamp)
end

# Whether a file's registration needs a date its record does not hold: one stating its own version,
# whose provider's time was not recorded when it was fetched.
function _needsdate(r::ProvenanceRecord, fileversioned::Bool)
    return isnothing(r.releasedate) && fileversioned && isnothing(r.modified)
end

function _needsdate(r::ProvenanceRecord, rec::DatasetRecord)
    return _needsdate(r, _fileversioned(rec))
end

# How a source is named to the registry: by its DOI, or, for one identified by its web page, by
# the key of its citation - a catalogue row's request name, under which `citations.bib` has it.
function _sourcename(facts::SourceFacts)
    return isnothing(facts.doi) ? (bibtex = facts.key,) : (doi = facts.doi,)
end

# Text or `nothing` where it is blank, as the registration keywords take an absent fact and a
# catalogue row writes one.
_ornothing(s::AbstractString) = isempty(s) ? nothing : String(s)

_ornothing(::Nothing) = nothing

# The keywords a file is registered with, as `link_register!` and `link_write!` take them, from its
# record and what is known of its source - a catalogue row's facts, or a DOI's: the source, the
# file's title, whether it is the source's own, the licence, the release, and the documentation's
# section quoting what the provider asks - for an extract also the request it was made to.
# `version` and `date` are the release's.
function _sourcekeywords(r::ProvenanceRecord, rec::DatasetRecord,
                         route::AbstractFileRoute; kw...)
    return _sourcekeywords(r, SourceFacts(rec), route; kw...)
end

function _sourcekeywords(r::ProvenanceRecord, facts::SourceFacts,
                         route::AbstractFileRoute;
                         version::AbstractString,
                         date::Union{Nothing, AbstractString})
    description = route isa GeneratedOnRequest ? _extractdescription(r) :
                  r.requestpage
    return (; _sourcename(facts)..., title = r.title,
            primary = something(r.primary, false),
            licence = something(r.licence, facts.licence, Some(nothing)),
            release_version = String(version), release_date = date,
            description = description)
end

# What an extract is, for its registration: the service's dataset and the request it was made to,
# where the record holds one, then where the documentation quotes what the provider asks.
function _extractdescription(r::ProvenanceRecord)
    text = if isnothing(r.request)
        "An extract of $(r.source)."
    else
        request = join(("$k = $(repr(r.request[k]))"
                        for k in sort!(collect(keys(r.request)))), ", ")
        "An extract of $(r.source) made to the request $request."
    end
    return isnothing(r.requestpage) ? text :
           text * " What its provider asks: " * r.requestpage
end

# The keywords the deposit a source's files are derived from is registered with: the dataset itself,
# file-less, its own where its identifier names the data. Its release is the row's, else its
# citation's in `citations.bib` (Natural Earth's), never a file's: the registry keeps a deposit's
# first release, so it must not depend on which file is registered first.
_depositkeywords(rec::DatasetRecord) = _depositkeywords(SourceFacts(rec))

function _depositkeywords(facts::SourceFacts)
    return (; _sourcename(facts)..., title = facts.title,
            primary = facts.identifies === :data,
            licence = facts.licence,
            release_version = isnothing(facts.version) ? nothing :
                              Provenancer._fullversion(facts.version),
            release_date = facts.releasedate,
            description = facts.requestpage)
end

# Whether a track says which release a spec's files come from, which only files the user gives
# may be told.
_hasrelease(::Nothing) = false

function _hasrelease(track::Track)
    return !isnothing(track.release_version) ||
           !isnothing(track.release_date)
end

# Refuse a track saying which release the files come from, given to a spec whose files are not all
# the user's own: for a fetched file the catalogue's row and the provider say it.
function _refusefetchedrelease(track, files)
    _hasrelease(track) || return nothing
    files isa AbstractVector && all(f -> f isa AbstractString, files) &&
        return nothing
    return error("`release_version` and `release_date` of a `Track` say which release files you " *
                 "give come from, by `file`, `files` or `directory`; this spec's files are fetched " *
                 "for it, and its catalogue row says which release they are. Leave them out.")
end
