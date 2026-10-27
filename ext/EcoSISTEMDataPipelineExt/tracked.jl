# SPDX-License-Identifier: LGPL-3.0-or-later

# --- A tracked spec's files through a FAIR data pipeline run ---------------------------------------
#
# The route methods `src/extensions.jl` declares, for a `Track` holding a run's handle. A fetch run -
# one whose working config writes a file's provenance record - fetches a spec's files by their
# untracked route and registers them as the registry can derive them truly: its source's file-less
# deposit, by `link_register!` in the source's namespace, is the run's one input, and every file -
# downloaded, extracted or given - and every file's record are its outputs, by `link_write!`. The
# registry derives each output of a run from all of the run's inputs, so a fetch run registers one
# source's files and refuses another's. Any other read takes each file and its record from the run's
# store by `link_read!`, never fetching. What each file is registered as is decided by the plain
# functions of `src/registration.jl`; these ask the registry and make the calls.

# What has been done for each run in this process, so it is done once: its authors and citations
# read, its deposit registered, each file put in place. Keyed by the handle itself.
const _DONE = IdDict{Any, Set{String}}()

# The deposit each fetch run registers files of, by its handle: the one source a run's outputs may
# be derived from.
const _DEPOSIT = IdDict{Any, String}()

_isdone(handle, key::AbstractString) = key in get(_DONE, handle, Set{String}())

function _markdone!(handle, key::AbstractString)
    return push!(get!(Set{String}, _DONE, handle), key)
end

# Read the catalogue's authors and citations into the run, once. Without them every author a
# citation names by an ORCID or ROR identifier reads as unregistered, and every question about a
# source answers false.
function _readcitations!(handle)
    _isdone(handle, "citations") && return nothing
    dir = joinpath(pkgdir(EcoSISTEM), "data", "catalogue")
    register_authors!(handle, joinpath(dir, "authors.toml"))
    register_citations!(handle, joinpath(dir, "citations.bib"))
    _markdone!(handle, "citations")
    return nothing
end

# The catalogue row of the source a tracked spec reads, which names its namespace and its facts. A
# file no row describes has no source to register it with.
function _sourcerow(spec)
    rec = EcoSISTEM._datasetrecord(_sourcetype(spec))
    isnothing(rec) &&
        error("$spec reads no catalogued source, so its files have no source to be registered " *
              "with; read it untracked.")
    return rec
end

_sourcetype(spec::RasterSpec) = spec.source

_sourcetype(::NaturalEarthSpec) = EcoSISTEM.NaturalEarthLevel

# The deposit a tracked spec's files are derived from, named in its source's namespace.
_depositof(track, rec) = rec.namespace * ":" * EcoSISTEM._depositname(track)

# Claim a fetch run for one deposit, refusing a second: every output of a run is derived from every
# input, so a second source's deposit would be shown as what the first's files came from.
function _claimdeposit!(track, rec)
    deposit = _depositof(track, rec)
    claimed = get!(_DEPOSIT, track.handle, deposit)
    claimed == deposit ||
        error("This fetch run registers the files of `$claimed`, so it cannot register " *
              "`$deposit`'s too: the registry derives every output of a run from all of its " *
              "inputs, so each file would be shown as derived from both. Fetch them in a run of " *
              "their own.")
    return deposit
end

# The name a file is registered, and asked about, under: the label of an output, which the run's
# `write:` entry resolves in the namespace it gives.
function _registrationname(track, path)
    return ConfigDataProduct(EcoSISTEM._dataname(track, path))
end

# The data products the bytes at `path` are registered as under the name `name`, in whatever
# namespace the config put them, and with `newest` only where they are that data product's newest
# version.
function _registeredas(handle, path, name; newest::Bool)
    id = DataPipeline.identify(handle, path)
    return [dp
            for obj in id.objects
            for dp in obj.data_products
            if dp.name == name && !(newest && dp.version != dp.newest)]
end

# The release version of the source the bytes of a data file are registered under as `name`, or
# `nothing` where they are not: a data product holding them with no source is not a registration of
# a file of ours, which always has one, but what a registration that failed part way left behind.
function _heldversion(handle, path, name; newest::Bool)
    for dp in _registeredas(handle, path, name, newest = newest)
        isnothing(dp.source) || return something(dp.source.version, "1.0.0")
    end
    return nothing
end

# Whether the bytes of a provenance record are registered as `name`, a data product with no source.
function _isrecorded(handle, path, name)
    return !isempty(_registeredas(handle, path, name, newest = true))
end

# Whether the newest version of a file's data product records the facts its catalogue row gives
# now, the release version being the one its bytes hold: a catalogue change makes a new version.
# The release date is given though not compared, since a source cannot be named without one.
function _samefacts(track, rec, r, path, route, held)
    kw = EcoSISTEM._sourcekeywords(r, rec, route, version = held,
                                   date = EcoSISTEM._registereddate(r, rec,
                                                                    nothing))
    return DataPipeline.isregistered(track.handle,
                                     _registrationname(track, path); kw...)
end

# Whether a file is registered already as it would be now, with its record: here, with a record of
# ours, its bytes the newest version of its data product with the facts the catalogue gives now, and
# its record's bytes registered too - so a fetch run with nothing else to do fetches and registers
# nothing. `entry` is the file's, one of `_specroutes(spec)`.
function _isregistered(track, rec, spec, entry)
    path, route = entry.path, entry.route
    isfile(path) || return false
    r = EcoSISTEM._recordwith(path, EcoSISTEM._currentfacts(spec, entry))
    isnothing(r) && return false
    handle = track.handle
    held = _heldversion(handle, path, EcoSISTEM._dataname(track, path),
                        newest = true)
    isnothing(held) && return false
    _samefacts(track, rec, r, path, route, held) || return false
    return _isrecorded(handle, EcoSISTEM._sidecarpath(path),
                       EcoSISTEM._sidecarname(track, path))
end

# Whether the run's working config writes both of a file's outputs, the file and its record.
function _writable(track, path)
    return DataPipeline.canwrite(track.handle,
                                 EcoSISTEM._dataname(track, path)) &&
           DataPipeline.canwrite(track.handle,
                                 EcoSISTEM._sidecarname(track, path))
end

# Refuse, before anything is fetched, a fetch run whose working config does not write what it would
# register for these files, writing out the entries to add.
function _checkwritable(track, rec, entries)
    lacking = filter(e -> !_writable(track, e.path), entries)
    isempty(lacking) && return nothing
    return error("This run's working config does not write what fetching " *
                 "$(join((basename(e.path) for e in lacking), ", ")) registers. Give it:\n" *
                 EcoSISTEM._writelines(track, rec.namespace,
                                       EcoSISTEM._filetype(first(lacking).path)))
end

# The version and release date a file is registered under, and its record rewritten where they or
# the provider's time it was dated by are new to it: the version its bytes hold already; else its
# source's own; else, that being held by other bytes, its source's own stamped with the day its
# provider changed it - asked of the server where the record does not say, if the server can
# vouch for these bytes - or failing that the day it was fetched or registered.
function _versionfor!(track, rec, r, path, route)
    handle = track.handle
    name = EcoSISTEM._dataname(track, path)
    held = _heldversion(handle, path, name, newest = false)
    isnothing(held) ||
        return (version = held,
                date = EcoSISTEM._registereddate(r, rec, nothing), record = r)
    full = EcoSISTEM._fullversion(r.version)
    kw = EcoSISTEM._sourcekeywords(r, rec, route, version = full,
                                   date = EcoSISTEM._registereddate(r, rec,
                                                                    nothing))
    taken = DataPipeline.isregistered(handle, _registrationname(track, path);
                                      kw...)
    taken || return (version = full,
            date = EcoSISTEM._registereddate(r, rec, nothing),
            record = r)
    if isnothing(r.modified) && !isnothing(r.url)
        modified = EcoSISTEM._headprovidertime(r.url, path)
        isnothing(modified) || (r = _rerecord!(path, "modified" => modified))
    end
    day = EcoSISTEM._stampday(r, Dates.today())
    version = EcoSISTEM._registeredversion(r.version, nothing, true, day)
    date = EcoSISTEM._registereddate(r, rec, day)
    r = _rerecord!(path, "version" => version, "release_date" => date)
    return (version = version, date = date, record = r)
end

# Rewrite a file's record of ours with `facts` in place of what it said, and return it read back.
function _rerecord!(path, facts::Pair...)
    sidecar = EcoSISTEM._sidecarpath(path)
    record = TOML.parsefile(sidecar)
    for (key, value) in facts
        record[key] = value
    end
    EcoSISTEM._writesidecar(path, record)
    return EcoSISTEM._ourrecordat(path)
end

# Register one file of a tracked spec with its run as an output, derived from its source's deposit -
# the deposit registered first, once a run - its record first brought up to date with what the
# catalogue says of it now; then that record, an output too. `entry` is the file's, one of
# `_specroutes(spec)`.
function _register!(track, spec, rec, entry)
    handle = track.handle
    path, route = entry.path, entry.route
    isfile(path) ||
        error("`$path` was fetched for $spec, yet is not where it was expected, so it cannot be " *
              "registered.")
    route === :given && EcoSISTEM._recordgiven(spec, entry)
    r = EcoSISTEM._refreshrecord!(path, EcoSISTEM._currentfacts(spec, entry))
    isnothing(r) &&
        error("`$(basename(path))` has no provenance record of EcoSISTEM's, so what it is cannot " *
              "be registered: delete its `.provenance.toml` to have one written.")
    if EcoSISTEM._needsdate(r, rec) && !isnothing(r.url)
        modified = EcoSISTEM._headprovidertime(r.url, path)
        isnothing(modified) || (r = _rerecord!(path, "modified" => modified))
    end
    release = _versionfor!(track, rec, r, path, route)
    kw = EcoSISTEM._sourcekeywords(release.record, rec, route,
                                   version = release.version,
                                   date = release.date)
    deposit = _claimdeposit!(track, rec)
    if !_isdone(handle, "deposit")
        link_register!(handle, deposit;
                       EcoSISTEM._depositkeywords(rec, date = release.date)...)
        _markdone!(handle, "deposit")
    end
    store = route === :published ? release.record.url : nothing
    link_write!(handle, EcoSISTEM._dataname(track, path); path = path,
                original_store = store, kw...)
    link_write!(handle, EcoSISTEM._sidecarname(track, path),
                path = EcoSISTEM._sidecarpath(path))
    return path
end

# A fetch run's call for a tracked spec: refuse a second source's files and a config that cannot
# write what it would register, return at once where every file is registered already as it would
# be now, and otherwise fetch the files by their untracked route and register each.
function EcoSISTEM._trackedfetch(track::Track{<:DataRegistryHandle}, spec)
    rec = _sourcerow(spec)
    entries = EcoSISTEM._specroutes(spec)
    paths = [e.path for e in entries]
    EcoSISTEM._checkstems(paths)
    _claimdeposit!(track, rec)
    _checkwritable(track, rec, entries)
    _readcitations!(track.handle)
    all(e -> _isregistered(track, rec, spec, e), entries) && return paths
    EcoSISTEM._fetchspec(spec)
    for e in entries
        _register!(track, spec, rec, e)
    end
    return paths
end

# Put one file and its record in place from the run's store, refusing a different file in the way
# rather than replacing it, and replacing a record of ours, whose file has matched.
function _providefile!(track, rec, path)
    handle = track.handle
    name = EcoSISTEM._dataname(track, path)
    sidecarname = EcoSISTEM._sidecarname(track, path)
    (DataPipeline.canread(handle, name) &&
     DataPipeline.canread(handle, sidecarname)) ||
        error("This run's working config does not read `$name` and its provenance record, " *
              "which a tracked spec reads only from the run's store. Add to its config:\n" *
              EcoSISTEM._readlines(track, rec.namespace))
    mkpath(dirname(abspath(path)))
    try
        link_read!(handle, name, path = path)
    catch e
        # DataPipeline refuses a file in the way, of other bytes than the run's, by this message
        # alone; anything else it refuses is passed on as it said it.
        (e isa ReadWriteException && occursin(" exists and is not ", e.msg)) ||
            rethrow()
        error("`$path` is there with other bytes than this run's `$name`, so it is not read in " *
              "its place: delete it and its `.provenance.toml` to have the run's copy put there.")
    end
    sidecar = EcoSISTEM._sidecarpath(path)
    if isfile(sidecar)
        t = TOML.tryparsefile(sidecar)
        (t isa AbstractDict && EcoSISTEM._ourrecord(t)) ||
            error("`$(basename(sidecar))` was not written by EcoSISTEM, so it is not replaced by " *
                  "this run's record `$sidecarname`: move it aside.")
        link_read!(handle, sidecarname, path = sidecar, force = true)
    else
        link_read!(handle, sidecarname, path = sidecar)
    end
    return nothing
end

# Any other read of a tracked spec: each file and its record put in place from the run's store,
# once a run.
function EcoSISTEM._trackedprovide!(track::Track{<:DataRegistryHandle}, spec,
                                    paths)
    rec = _sourcerow(spec)
    for path in paths
        key = "provided:" * path
        _isdone(track.handle, key) && continue
        _providefile!(track, rec, path)
        _markdone!(track.handle, key)
    end
    return nothing
end

# Why one file of a tracked spec is still to do in its run, or `nothing`: nothing where the config
# reads it; in a fetch run, what the registry holds of it; in a model run, why the config cannot.
# `entry` is the file's, one of `_specroutes(spec)`.
function _missingreason(track, rec, spec, entry)
    handle = track.handle
    path, route = entry.path, entry.route
    name = EcoSISTEM._dataname(track, path)
    sidecarname = EcoSISTEM._sidecarname(track, path)
    (DataPipeline.canread(handle, name) &&
     DataPipeline.canread(handle, sidecarname)) &&
        return nothing
    if DataPipeline.canwrite(handle, sidecarname)
        present = isfile(path)
        r = present ?
            EcoSISTEM._recordwith(path, EcoSISTEM._currentfacts(spec, entry)) :
            nothing
        held = isnothing(r) ? nothing :
               _heldversion(handle, path, name, newest = true)
        same = !isnothing(held) && _samefacts(track, rec, r, path, route, held)
        recorded = present &&
                   _isrecorded(handle, EcoSISTEM._sidecarpath(path),
                               sidecarname)
        named = DataPipeline.isregistered(handle,
                                          _registrationname(track, path))
        return EcoSISTEM._fetchreason(present = present, held = held,
                                      same = same,
                                      recorded = recorded, named = named)
    end
    registered = try
        DataPipeline.isregistered(handle, ConfigDataProduct(name))
    catch e
        e isa ConfigFileException || rethrow()
        nothing
    end
    return EcoSISTEM._modelreason(registered, track, rec.namespace)
end

# What `missingfiles` reports of a tracked spec: its files still to do in its run, and why.
function EcoSISTEM._trackedmissing(track::Track{<:DataRegistryHandle}, spec)
    rec = _sourcerow(spec)
    _readcitations!(track.handle)
    todo = @NamedTuple{file::String, reason::String}[]
    for e in EcoSISTEM._specroutes(spec)
        reason = _missingreason(track, rec, spec, e)
        isnothing(reason) || push!(todo, (file = e.path, reason = reason))
    end
    return todo
end

# Whether the run's working config writes what a fetch run registers for each of a tracked spec's
# files, for a dry run.
function EcoSISTEM._trackedcovers(track::Track{<:DataRegistryHandle}, spec)
    return [_writable(track, e.path) for e in EcoSISTEM._specroutes(spec)]
end
