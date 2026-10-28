# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a tracked spec's files are registered as with a FAIR data pipeline run, decided from plain
# data so that it is tested with no registry: the data product names, the route each file comes by,
# the facts of its source, the version a file takes whose bytes differ from those registered, the
# config lines a run lacks, and why a file is still to do. The calls that ask a registry and
# register are the DataPipeline extension's.

# The stem a file's data products are named by: its name without its extension, which the data
# store adds.
_filestem(path::AbstractString) = first(splitext(basename(path)))

# The data products a tracked file gives rise to under its track's root: the file itself, its
# provenance record, and the file-less deposit an extract's files are derived from.
function _dataname(track::Track, path::AbstractString)
    return track.root * "/data/" * _filestem(path)
end

function _sidecarname(track::Track, path::AbstractString)
    return track.root * "/provenance/" * _filestem(path)
end

_depositname(track::Track) = track.root * "/deposit"

# Refuse, before anything is fetched, a file whose stem cannot end a data product name.
function _checkstems(paths)
    for path in paths
        stem = _filestem(path)
        occursin(_PRODUCTNAME, stem) && !occursin('/', stem) ||
            error("`$(basename(path))` cannot be tracked: a data product is named by the file's " *
                  "name without its extension, `$stem`, which may hold only ASCII letters, " *
                  "digits, `_`, `-` and `.`, not beginning or ending with `.`. Rename the file.")
    end
    return nothing
end

# Each file a spec reads, where it is or would be, with the route it comes by - `:published` for a
# file fetched as its provider publishes it (by its source's own resolution, or from a URL),
# `:extract` for one a service makes to a request, `:given` for a file the user names - the layer
# `code` it holds where one is known, and the `url` it is fetched from where it has one.
_specroutes(spec::RasterSpec) = _specroutes(spec, spec.files)

function _specroutes(spec::RasterSpec, ::Nothing)
    files = _sourcefiles(spec.source, spec.code; _readkw(spec)...)
    isnothing(files) &&
        return [(path = p, route = :published, code = nothing, url = nothing)
                for p in _specpaths(spec)]
    return [(path = f.path, route = :published, code = f.code, url = f.url)
            for f in files]
end

function _specroutes(spec::RasterSpec, files::AbstractVector)
    code = spec.code isa AbstractVector ? nothing : spec.code
    return [(path = _localpath(f), route = _entryroute(f), code = code,
             url = _entryurl(f))
            for f in files]
end

function _specroutes(spec::NaturalEarthSpec)
    url = _nesource(_checklevel(spec.level)).url
    return [(path = p, route = :published, code = nothing, url = url)
            for p in _specpaths(spec)]
end

_entryroute(::AbstractString) = :given

_entryroute(::CachedAsset) = :published

_entryroute(::CDSRequest) = :extract

_entryurl(asset::CachedAsset) = asset.url

_entryurl(::Any) = nothing

# What the catalogue says now of one file of a tracked spec, as its record holds it: the catalogue's
# facts for the route it comes by, where it is fetched from, the release a track gives a file the
# user gave, and the version a file states of itself. `entry` is one of `_specroutes(spec)`.
function _currentfacts(spec::RasterSpec, entry)
    layer = isnothing(entry.code) ? nothing :
            _layerornothing(spec.source, entry.code)
    facts = _cataloguerecord(spec.source, layer, entry.path,
                             route = entry.route)
    isnothing(entry.url) || (facts["url"] = entry.url)
    track = spec.track
    if entry.route === :given && !isnothing(track)
        isnothing(track.release_version) ||
            (facts["version"] = track.release_version)
        isnothing(track.release_date) ||
            (facts["release_date"] = track.release_date)
    end
    return facts
end

function _currentfacts(spec::NaturalEarthSpec, entry)
    facts = _cataloguerecord(NaturalEarthLevel, nothing, entry.path,
                             route = :published)
    facts["url"] = entry.url
    version = _fileversion(datasetinfo(NaturalEarthLevel), entry.path)
    isnothing(version) || (facts["version"] = version)
    return facts
end

# A file's record of ours as it would be with `facts` in place of what it says of them, read without
# writing it; `nothing` where the file has no record of ours.
function _recordwith(path::AbstractString, facts::AbstractDict)
    sidecar = _sidecarpath(path)
    isfile(sidecar) || return nothing
    t = TOML.tryparsefile(sidecar)
    (t isa AbstractDict && _ourrecord(t)) || return nothing
    return _recordfromtoml(merge(t, facts))
end

# Bring a file's record of ours up to date with `facts`, what the catalogue says of it now, before it
# is registered, rewriting it only where they differ: a record written by an earlier version of this
# package may lack facts the registration needs, or hold ones the catalogue has since corrected,
# and the record and the registration must agree. What the fetch itself recorded - when, the size
# and checksum, the request, its job - is kept. Returns the record, or `nothing` for a file with no
# record of ours.
function _refreshrecord!(path::AbstractString, facts::AbstractDict)
    sidecar = _sidecarpath(path)
    isfile(sidecar) || return nothing
    t = TOML.tryparsefile(sidecar)
    (t isa AbstractDict && _ourrecord(t)) || return nothing
    if any(get(t, k, nothing) != v for (k, v) in facts)
        merge!(t, facts)
        t["writer"] = "EcoSISTEM $(pkgversion(EcoSISTEM))"
        _writesidecar(path, t)
    end
    return _ourrecordat(path)
end

# A spec's files fetched by its own route, as an untracked spec's are, for a tracked fetch run to
# register after.
_fetchspec(spec::RasterSpec) = _fetchedpaths(spec, spec.files)

function _fetchspec(spec::NaturalEarthSpec)
    return [assetpath(_nesource(_checklevel(spec.level)))]
end

# The provenance record of a file the user gave a spec, written beside it where it has none: what
# the catalogue says of it, the release its track names, and no fetch time, which is not known.
# `entry` is the file's, one of `_specroutes(spec)`.
function _recordgiven(spec::RasterSpec, entry)
    isfile(entry.path) ||
        error("`$(entry.path)`, a file given to $spec, is not there, so it cannot be registered.")
    _recordpresent(entry.path) do
        return merge!(_baserecord(entry.path, :habitat),
                      _currentfacts(spec, entry))
    end
    return entry.path
end

# A source's layer by its code, or `nothing` where the catalogue has none of that code.
function _layerornothing(source::Type, code)
    return try
        layerinfo(source, code)
    catch
        nothing
    end
end

# A release version as the registry writes it, in full: `2.1` as `2.1.0`, and `1.0.0` for none.
function _fullversion(version::AbstractString)
    return isempty(version) ? "1.0.0" :
           string(VersionNumber(version))
end

# The release version a file is registered under: the version its bytes are registered under
# already (`held`), so a third fetch of a changed file adds nothing; else its source's own; or,
# where that is taken by other bytes - its provider changed the file without a new release -
# the same with build metadata `+YYYYMMDD` from `day`.
function _registeredversion(version::AbstractString,
                            held::Union{Nothing, AbstractString}, taken::Bool,
                            day::Dates.Date)
    isnothing(held) || return String(held)
    full = _fullversion(version)
    taken || return full
    return full * "+" * Dates.format(day, "yyyymmdd")
end

# Whether a registered release version is a stamped one, a provider's change dated by its day.
_isstamped(version::AbstractString) = occursin('+', version)

# The day a file's bytes are dated by when they take a stamp: when its provider last changed it,
# else when it was fetched, else `today`, the registration's day, for a file the user gave.
function _stampday(r::InputRecord, today::Dates.Date)
    return Dates.Date(something(r.modified, r.fetched, today))
end

# Whether a source's files state their own versions rather than its row giving one: Natural
# Earth's zips, each naming its own.
_fileversioned(rec::DatasetRecord) = rec.format === :Shapefile

# The release date a file is registered with: a stamped version's day; else its record's, the
# row's or the one its track gave; else, for a file stating its own version, when its provider last
# changed it, or failing that when it was fetched; else `nothing`, the citation's date standing.
function _registereddate(r::InputRecord, rec::DatasetRecord,
                         stamp::Union{Nothing, Dates.Date})
    isnothing(stamp) || return string(stamp)
    isempty(r.releasedate) || return r.releasedate
    if _fileversioned(rec)
        day = something(r.modified, r.fetched, Some(nothing))
        isnothing(day) || return string(Dates.Date(day))
    end
    return nothing
end

# Whether a file's registration needs a date its record does not hold: one stating its own version,
# whose provider's time was not recorded when it was fetched.
function _needsdate(r::InputRecord, rec::DatasetRecord)
    return isempty(r.releasedate) && _fileversioned(rec) &&
           isnothing(r.modified)
end

# How a source is named to the registry: by its DOI, or, for a row identified by its web page, by
# the key of its citation in `citations.bib`, which is the row's request name.
function _sourcename(rec::DatasetRecord)
    return isempty(rec.doi) ? (bibtex = rec.request,) : (doi = rec.doi,)
end

# Text or `nothing` where it is blank, as the registration keywords take an absent fact.
_ornothing(s::AbstractString) = isempty(s) ? nothing : String(s)

_ornothing(::Nothing) = nothing

# The keywords a file is registered with, as `link_register!` and `link_write!` take them, from its
# record and its source's row: the source, the file's title, whether it is the source's own, the
# licence, the release, and the documentation's section quoting what the provider asks - for an
# extract also the request it was made to. `version` and `date` are the release's.
function _sourcekeywords(r::InputRecord, rec::DatasetRecord, route::Symbol;
                         version::AbstractString,
                         date::Union{Nothing, AbstractString})
    description = route === :extract ? _extractdescription(r) :
                  _ornothing(r.requestpage)
    return (; _sourcename(rec)..., title = r.title,
            primary = something(r.primary, false),
            licence = _ornothing(something(_ornothing(r.licence), rec.licence)),
            release_version = String(version), release_date = date,
            description = description)
end

# What an extract is, for its registration: the service's dataset and the request it was made to,
# then where the documentation quotes what the provider asks.
function _extractdescription(r::InputRecord)
    request = isnothing(r.request) ? "" :
              join(("$k = $(repr(r.request[k]))"
                    for k in sort!(collect(keys(r.request)))),
                   ", ")
    text = "An extract of $(r.dataset) made to the request $request."
    return isempty(r.requestpage) ? text :
           text * " What its provider asks: " * r.requestpage
end

# The keywords the deposit a source's files are derived from is registered with: the dataset itself,
# file-less, its own where its identifier names the data. Its release is the row's, else its
# citation's in `citations.bib` (Natural Earth's), never a file's: the registry keeps a deposit's
# first release, so it must not depend on which file is registered first.
function _depositkeywords(rec::DatasetRecord)
    return (; _sourcename(rec)..., title = rec.title,
            primary = rec.identifies === :data,
            licence = _ornothing(rec.licence),
            release_version = isempty(rec.version) ? nothing :
                              _fullversion(rec.version),
            release_date = _ornothing(rec.releasedate),
            description = _ornothing(_requestpage(rec)))
end

# The file type of a file's data product, its extension, as a `write:` entry names it.
_filetype(path::AbstractString) = lstrip(last(splitext(path)), '.')

# What a fetch config needs for a track's files: its source's namespace as the run's default for
# inputs and outputs, since the run registers that source's deposit and nothing else, and `write:`
# entries for the files, of `filetype`, and their provenance records.
function _writelines(track::Track, namespace::AbstractString,
                     filetype::AbstractString)
    return join(["run_metadata:", "  default_input_namespace: $namespace",
                    "  default_output_namespace: $namespace", "write:",
                    "- data_product: $(track.root)/data/*",
                    "  description: A file of the source tracked under $(track.root)",
                    "  file_type: $filetype",
                    "- data_product: $(track.root)/provenance/*",
                    "  description: EcoSISTEM's provenance record of a file under $(track.root)",
                    "  file_type: toml"], "\n")
end

# The `read:` entry a model config needs for a track's files and their records, both in their
# source's namespace: `*` matches one segment of a name, so `<root>/*/*` takes the files and the
# records and leaves the deposit, which has no file.
function _readlines(track::Track, namespace::AbstractString)
    return join(["read:", "- data_product: $(track.root)/*/*", "  use:",
                    "    namespace: $namespace"], "\n")
end

# Why a fetch run still has a file to do, from what its registry answered, or `nothing` where it has
# none: whether the file is `present`; the version its bytes hold under its name (`held`, `nothing`
# for none); whether that version records the facts the catalogue gives now (`same`); whether its
# record's bytes are registered (`recorded`); and whether its name is registered at all (`named`).
function _fetchreason(; present::Bool, held::Union{Nothing, AbstractString},
                      same::Bool, recorded::Bool, named::Bool)
    present || return "not held: the fetch run fetches and registers it"
    if isnothing(held)
        named &&
            return "registered with other bytes, its provider having changed it: the fetch " *
                   "run registers these under a stamped version"
        return "not registered: the fetch run registers it"
    end
    same ||
        return "registered with other facts, a catalogue change: the fetch run registers a " *
               "new version"
    recorded ||
        return "its provenance record is not registered: the fetch run registers it"
    return nothing
end

# Why a model run cannot read a file, from whether its registry holds it (`registered`, `nothing`
# where no `read:` entry covers its name, so the registry could not be asked).
function _modelreason(registered::Union{Nothing, Bool}, track::Track,
                      namespace::AbstractString)
    isnothing(registered) &&
        return "no `read:` entry covers it; add\n" *
               _readlines(track, namespace)
    registered &&
        return "registered since this job was prepared: `fair pull`, then `fair run` again"
    return "not registered in namespace `$namespace`: run the fetch run first"
end
