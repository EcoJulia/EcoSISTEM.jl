# SPDX-License-Identifier: LGPL-3.0-or-later

# --- Tracked files through a FAIR data pipeline run --------------------------------------------------
#
# The route methods `src/extensions.jl` declares for a spec whose `Track` holds a run's handle, and
# the two `Provenancer` declares for a file of the user's own that `record_provenance` records. A fetch
# run - one whose working config writes a file's provenance record - fetches a spec's files by their
# untracked route and registers them as the registry can derive them truly: its source's file-less
# deposit, by `link_register!` in the source's namespace, is the run's one input, and every file -
# downloaded, extracted or given - and every file's record are its outputs, by `link_write!`. The
# registry derives each output of a run from all of the run's inputs, so a fetch run registers one
# source's files and refuses another's. Any other read takes each file and its record from the run's
# store by `link_read!`, never fetching. What each file is registered as is decided by the plain
# functions of `src/registration.jl` and `src/Provenancer/naming.jl`; these ask the registry and
# make the calls. The source a spec's files are registered with is its catalogue row's; a file of
# the user's own is registered from the facts its own record keeps.

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

# The deposit a tracked file is derived from, named in its source's namespace, or in the run's
# default input namespace for a source naming none.
function _depositof(track, namespace::Union{Nothing, AbstractString})
    name = Provenancer._depositname(track)
    return isnothing(namespace) ? name : namespace * ":" * name
end

# The name a file registered as its source itself is registered under, in its source's namespace,
# or the run's default input namespace for a source naming none: no `write:` entry resolves it.
function _fileof(track, namespace::Union{Nothing, AbstractString}, path)
    name = Provenancer._dataname(track, path)
    return isnothing(namespace) ? name : namespace * ":" * name
end

# Claim a fetch run for one source, refusing a second: every output of a run is derived from every
# input, so a second source - another deposit, or another file registered as its source - would be
# shown as what the first's files came from. A source is claimed by its deposit's name, or by the
# file's own where the file is the deposit.
function _claimsource!(track, facts::SourceFacts, path)
    source = _sourceclaim(facts.deposit, track, facts.namespace, path)
    claimed = get!(_DEPOSIT, track.handle, source)
    claimed == source ||
        error("This fetch run registers `$claimed` as its source, so it cannot register " *
              "`$source` too: the registry derives every output of a run from all of its " *
              "inputs, so each file would be shown as derived from both. Fetch them in a run of " *
              "their own.")
    return source
end

# What a run claims for a source: its deposit's name, or the file's own where the file is the deposit.
function _sourceclaim(::AbstractDeposit, track, namespace, path)
    return _depositof(track, namespace)
end

function _sourceclaim(::DepositIsFile, track, namespace, path)
    return _fileof(track, namespace, path)
end

# Why a source's files cannot be registered with it: its facts name it neither by a DOI nor by a
# citation key, and the registry takes nothing else as a source.
function _nosource(facts)
    return "`$(facts.source)`'s facts name neither a DOI nor the key of a citation, so the " *
           "registry cannot be told where its files come from. Give them a `doi`, or the `key` " *
           "of a citation registered with the run. A file of your own making is an output of " *
           "the run that made it, registered with DataPipeline's `link_write!`."
end

# One file to register, read or ask about, with everything the calls need of it: its `path`, the
# `route` it comes by, the layer `code` it holds and the `url` it came from where known, where its
# `sidecar` is, and `current`, what is known of it now that its record is brought up to date with -
# the catalogue's facts for a spec's file, the release its track names for one the user gives.
function _item(entry, current::AbstractDict;
               sidecar::AbstractString = Provenancer.sidecarpath(entry.path))
    return (; entry..., sidecar = sidecar, current = current)
end

function _items(spec)
    return [_item(e, EcoSISTEM._currentfacts(spec, e))
            for e in EcoSISTEM._specroutes(spec)]
end

# The name a file is registered, and asked about, under: the label of an output, which the run's
# `write:` entry resolves in the namespace it gives.
function _registrationname(track, path)
    return ConfigDataProduct(Provenancer._dataname(track, path))
end

# The same for a file whose source's `facts` are known: the label of an output, or the file's own
# name where it is registered as its source itself, which no `write:` entry resolves.
function _registrationname(track, path, facts::SourceFacts)
    return _registrationname(facts.deposit, track, path, facts.namespace)
end

function _registrationname(::AbstractDeposit, track, path, namespace)
    return _registrationname(track, path)
end

function _registrationname(::DepositIsFile, track, path, namespace)
    return _fileof(track, namespace, path)
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

# Whether the newest version of a file's data product records the facts known of its source now
# (`src`, as `EcoSISTEM._trackedsource` gives it), the release version being the one its bytes
# hold: a change in the facts makes a new version. The release date is given though not compared,
# since a source cannot be named without one.
function _samefacts(track, src, r, path, route, held)
    kw = EcoSISTEM._sourcekeywords(r, src.facts, route, version = held,
                                   date = EcoSISTEM._registereddate(r,
                                                                    src.fileversioned,
                                                                    nothing))
    return DataPipeline.isregistered(track.handle,
                                     _registrationname(track, path, src.facts);
                                     kw...)
end

# Whether a file is registered already as it would be now, with its record: here, with a record of
# ours, its bytes the newest version of its data product with the facts known now, and its record's
# bytes registered too - so a fetch run with nothing else to do fetches and registers nothing.
function _isregistered(track, src, item)
    path, route = item.path, item.route
    isfile(path) || return false
    r = Provenancer._recordwith(path, item.current, sidecar = item.sidecar)
    isnothing(r) && return false
    handle = track.handle
    held = _heldversion(handle, path, Provenancer._dataname(track, path),
                        newest = true)
    isnothing(held) && return false
    _samefacts(track, src, r, path, route, held) || return false
    return _isrecorded(handle, item.sidecar,
                       Provenancer._sidecarname(track, path))
end

# Whether the run's working config writes a file's record: what marks a fetch run, however the
# file's source is registered.
function _writable(track, path)
    return DataPipeline.canwrite(track.handle,
                                 Provenancer._sidecarname(track, path))
end

# Whether it writes everything registering the file writes: the file and its record, or the record
# alone for a file registered as its source itself, which is the run's input.
function _writesall(::AbstractDeposit, track, path)
    return _writable(track, path) &&
           DataPipeline.canwrite(track.handle,
                                 Provenancer._dataname(track, path))
end

_writesall(::DepositIsFile, track, path) = _writable(track, path)

# Whether the run's working config reads both, the file and its record.
function _readable(track, path)
    return DataPipeline.canread(track.handle,
                                Provenancer._dataname(track, path)) &&
           DataPipeline.canread(track.handle,
                                Provenancer._sidecarname(track, path))
end

# Refuse, before anything is fetched, a fetch run whose working config does not write what it would
# register for these files, writing out the entries to add.
function _checkwritable(track, facts::SourceFacts, items)
    lacking = filter(i -> !_writesall(facts.deposit, track, i.path), items)
    isempty(lacking) && return nothing
    return error("This run's working config does not write what fetching " *
                 "$(join((basename(i.path) for i in lacking), ", ")) registers. Give it:\n" *
                 Provenancer._writelines(track, facts.namespace,
                                         Provenancer._filetype(first(lacking).path),
                                         facts.deposit))
end

# The version and release date a file is registered under, and its record rewritten where they or
# the provider's time it was dated by are new to it: the version its bytes hold already; else its
# source's own; else, that being held by other bytes, its source's own stamped with the day its
# provider changed it - asked of the server where the record does not say, if the server can
# vouch for these bytes - or failing that the day it was fetched or registered.
function _versionfor!(track, src, r, path, route, sidecar)
    handle = track.handle
    name = Provenancer._dataname(track, path)
    held = _heldversion(handle, path, name, newest = false)
    isnothing(held) ||
        return (version = held,
                date = EcoSISTEM._registereddate(r, src.fileversioned, nothing),
                record = r)
    full = Provenancer._fullversion(r.version)
    kw = EcoSISTEM._sourcekeywords(r, src.facts, route, version = full,
                                   date = EcoSISTEM._registereddate(r,
                                                                    src.fileversioned,
                                                                    nothing))
    taken = DataPipeline.isregistered(handle,
                                      _registrationname(track, path, src.facts);
                                      kw...)
    taken || return (version = full,
            date = EcoSISTEM._registereddate(r, src.fileversioned, nothing),
            record = r)
    if isnothing(r.modified) && !isnothing(r.url)
        modified = Provenancer._headprovidertime(r.url, path)
        isnothing(modified) ||
            (r = _rerecord!(path, "modified" => modified, sidecar = sidecar))
    end
    day = Provenancer._stampday(r, Dates.today())
    version = Provenancer._registeredversion(r.version, nothing, true, day)
    date = EcoSISTEM._registereddate(r, src.fileversioned, day)
    r = _rerecord!(path, "version" => version, "release_date" => date,
                   sidecar = sidecar)
    return (version = version, date = date, record = r)
end

# Rewrite a file's record of ours at `sidecar` with `facts` in place of what it said, and return it
# read back.
function _rerecord!(path, facts::Pair...; sidecar::AbstractString)
    record = TOML.parsefile(sidecar)
    for (key, value) in facts
        record[key] = value
    end
    Provenancer.writesidecar(path, record, sidecar = sidecar)
    return Provenancer.recordat(path, sidecar = sidecar)
end

# Register one tracked file with its run as its source's deposit says, its record first brought up
# to date with what is known of it now; then that record, an output of the run.
function _register!(track, src, item)
    handle = track.handle
    path, route, sidecar = item.path, item.route, item.sidecar
    isfile(path) ||
        error("`$path` was to be registered, yet is not where it was expected, so it cannot be.")
    r = Provenancer._refreshrecord!(path, item.current, sidecar = sidecar)
    isnothing(r) &&
        error("`$(basename(path))` has no provenance record of EcoSISTEM's, so what it is cannot " *
              "be registered: delete its `.provenance.toml` to have one written, or record it " *
              "with `Provenancer.record_provenance`.")
    if EcoSISTEM._needsdate(r, src.fileversioned) && !isnothing(r.url)
        modified = Provenancer._headprovidertime(r.url, path)
        isnothing(modified) ||
            (r = _rerecord!(path, "modified" => modified, sidecar = sidecar))
    end
    release = _versionfor!(track, src, r, path, route, sidecar)
    kw = EcoSISTEM._sourcekeywords(release.record, src.facts, route,
                                   version = release.version,
                                   date = release.date)
    _registerfile!(src.facts.deposit, track, src, path, release.record.url, kw)
    link_write!(handle, Provenancer._sidecarname(track, path), path = sidecar)
    return path
end

# Register a file, its source's keywords `kw` given: as an output of the run derived from its
# source's deposit, the deposit registered first, once a run, from the source's facts or those of a
# deposit cited apart; or, where the file is the deposit, the file itself as the run's one input.
# `url` is where the file was fetched from, as its record says: a provider's address for a
# download, the one the user gave for a file they fetched themselves, none for an extract or a file
# given by hand.
function _registerfile!(deposit::AbstractDeposit, track, src, path, url, kw)
    handle = track.handle
    name = _claimsource!(track, src.facts, path)
    if !_isdone(handle, "deposit")
        facts = Provenancer._depositfacts(deposit, src.facts)
        link_register!(handle, name; EcoSISTEM._depositkeywords(facts)...)
        _markdone!(handle, "deposit")
    end
    link_write!(handle, Provenancer._dataname(track, path); path = path,
                original_store = url, kw...)
    return nothing
end

function _registerfile!(::DepositIsFile, track, src, path, url, kw)
    name = _claimsource!(track, src.facts, path)
    link_register!(track.handle, name; path = path, original_store = url, kw...)
    return nothing
end

# A fetch run's call for a tracked spec: refuse a source the registry cannot be told of, a second
# source's files and a config that cannot write what it would register; return at once where every
# file is registered already as it would be now; otherwise fetch the files by their untracked route,
# record each the user gave, and register each.
function EcoSISTEM._trackedfetch(track::Track{<:DataRegistryHandle}, spec)
    src = EcoSISTEM._trackedsource(spec)
    Provenancer._namessource(src.facts) || error(_nosource(src.facts))
    items = _items(spec)
    paths = [i.path for i in items]
    Provenancer._checkstems(paths)
    foreach(i -> _claimsource!(track, src.facts, i.path), items)
    _checkwritable(track, src.facts, items)
    _readcitations!(track.handle)
    all(i -> _isregistered(track, src, i), items) && return paths
    EcoSISTEM._fetchspec(spec)
    for i in items
        i.route isa OtherOrigin && EcoSISTEM._recordgiven(spec, i)
        _register!(track, src, i)
    end
    return paths
end

# Put one file and its record in place from the run's store, refusing a different file in the way
# rather than replacing it, and replacing a record of ours, whose file has matched. `namespace` is
# the source's, for the config entry a run that cannot read them is told to add.
function _providefile!(track, namespace::Union{Nothing, AbstractString}, path,
                       sidecar)
    handle = track.handle
    name = Provenancer._dataname(track, path)
    sidecarname = Provenancer._sidecarname(track, path)
    _readable(track, path) ||
        error("This run's working config does not read `$name` and its provenance record, " *
              "which a tracked file is read only from the run's store. Add to its config:\n" *
              Provenancer._readentry(track, namespace))
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
    if isfile(sidecar)
        t = TOML.tryparsefile(sidecar)
        (t isa AbstractDict && Provenancer._isourrecord(t)) ||
            error("`$(basename(sidecar))` was not written by EcoSISTEM, so it is not replaced by " *
                  "this run's record `$sidecarname`: move it aside.")
        link_read!(handle, sidecarname, path = sidecar, force = true)
    else
        mkpath(dirname(abspath(sidecar)))
        link_read!(handle, sidecarname, path = sidecar)
    end
    return nothing
end

# Any other read of a tracked spec: each file and its record put in place from the run's store,
# once a run.
function EcoSISTEM._trackedprovide!(track::Track{<:DataRegistryHandle}, spec,
                                    paths)
    rec = EcoSISTEM._catalogued(spec)
    namespace = isnothing(rec) ? nothing : EcoSISTEM._ornothing(rec.namespace)
    for path in paths
        key = "provided:" * path
        _isdone(track.handle, key) && continue
        _providefile!(track, namespace, path, Provenancer.sidecarpath(path))
        _markdone!(track.handle, key)
    end
    return nothing
end

# Why one tracked file is still to do in its run, or `nothing`: nothing where the config reads it;
# in a fetch run, what the registry holds of it; in a model run, why the config cannot. `src` is the
# source as `EcoSISTEM._trackedsource` gives it for a spec's file, or a file of the user's own
# record's facts, `nothing` for one with no record of ours yet.
function _missingreason(track, src, item; fetched::Bool = true)
    handle = track.handle
    path, route, sidecar = item.path, item.route, item.sidecar
    name = Provenancer._dataname(track, path)
    sidecarname = Provenancer._sidecarname(track, path)
    _readable(track, path) && return nothing
    if DataPipeline.canwrite(handle, sidecarname)
        present = isfile(path)
        isnothing(src) &&
            return present ?
                   "not recorded: `Provenancer.record_provenance(path, doi or facts, role = ...)` " *
                   "records it, and the fetch run registers it" :
                   "not there: a file of your own is not fetched, so put it in place and " *
                   "record it"
        r = present ?
            Provenancer._recordwith(path, item.current, sidecar = sidecar) :
            nothing
        held = isnothing(r) ? nothing :
               _heldversion(handle, path, name, newest = true)
        same = !isnothing(held) && _samefacts(track, src, r, path, route, held)
        recorded = present && _isrecorded(handle, sidecar, sidecarname)
        named = DataPipeline.isregistered(handle,
                                          _registrationname(track, path,
                                                            src.facts))
        return Provenancer._fetchreason(present = present, held = held,
                                        same = same,
                                        recorded = recorded, named = named,
                                        deposit = src.facts.deposit,
                                        fetched = fetched)
    end
    registered = try
        DataPipeline.isregistered(handle, ConfigDataProduct(name))
    catch e
        e isa ConfigFileException || rethrow()
        nothing
    end
    namespace = isnothing(src) ? nothing : src.facts.namespace
    return Provenancer._modelreason(registered, track, namespace)
end

# What `missingfiles` reports of a tracked spec: its files still to do in its run, and why.
function EcoSISTEM._trackedmissing(track::Track{<:DataRegistryHandle}, spec)
    src = EcoSISTEM._trackedsource(spec)
    _readcitations!(track.handle)
    todo = @NamedTuple{file::String, reason::String}[]
    for i in _items(spec)
        reason = _missingreason(track, src, i)
        isnothing(reason) || push!(todo, (file = i.path, reason = reason))
    end
    return todo
end

# Register a file recorded already - by `record_provenance`, or by a download - from its record `r`,
# as an output of the fetch run derived from its source's deposit, unless it is registered already as
# it would be now; return its record as it then stands. `route` is the one it was recorded with,
# `code` the layer it holds and `url` where it was fetched from, where known.
function _registerrecorded!(track, r, route; code = nothing, url = nothing)
    path, sidecar = r.location, r.sidecar
    facts = SourceFacts(r)
    Provenancer._namessource(facts) || error(_nosource(facts))
    src = (facts = facts, fileversioned = false)
    entry = (path = path, route = route, code = code, url = url)
    item = _item(entry,
                 EcoSISTEM._trackrelease!(Dict{String, Any}(), track, route),
                 sidecar = sidecar)
    _claimsource!(track, facts, path)
    _checkwritable(track, facts, [item])
    _readcitations!(track.handle)
    _isregistered(track, src, item) && return r
    _register!(track, src, item)
    return Provenancer.recordat(path, sidecar = sidecar)
end

# Put the run's copy of a recorded file and its record in place from the store, once a run, and
# return the record.
function _providerecorded!(track, path, sidecar)
    key = "provided:" * path
    if !_isdone(track.handle, key)
        _providefile!(track, nothing, path, sidecar)
        _markdone!(track.handle, key)
    end
    return Provenancer.recordat(path, sidecar = sidecar)
end

# Why a run cannot take a recorded file: its working config neither writes nor reads the file and
# its record, with the lines a fetch run and a model run each need, the first as the source's
# `deposit` has it registered.
function _neitherrun(track, path, namespace, deposit)
    return "This run's working config neither writes nor reads " *
           "`$(Provenancer._dataname(track, path))` and its provenance record. A fetch " *
           "run, which records and registers the file, needs:\n" *
           Provenancer._writelines(track, namespace,
                                   Provenancer._filetype(path), deposit) *
           "\nand a model run, which reads them from the run's store, needs:\n" *
           Provenancer._readentry(track, namespace)
end

# What `record_provenance` does for a file of the user's own given a run's track: in a fetch run - one
# whose working config writes the file's record - the file is recorded as an untracked call records
# it, once, and registered from its record's facts; in a model run - one whose config reads the
# record - the run's copy of the file and its record are put in place from the store; a run whose
# config does neither is refused, with the lines each kind needs. `source` is the facts, or the DOI
# whose metadata gives them; the rest are `record_provenance`'s keywords.
function Provenancer._trackedrecord(track::Track{<:DataRegistryHandle}, path,
                                    source; role::Symbol,
                                    sidecar::AbstractString,
                                    route::EcoSISTEM.AbstractFileRoute = VersionOfRecord(),
                                    kw...)
    path = String(path)
    Provenancer._checkstems([path])
    if _writable(track, path)
        isfile(path) ||
            error("`$path` is not there, so it cannot be recorded or registered: a file of your " *
                  "own is not fetched.")
        r = Provenancer.record_provenance(path, source; role = role,
                                          route = route, sidecar = sidecar,
                                          kw...)
        return _registerrecorded!(track, r, route,
                                  code = get(kw, :code, nothing),
                                  url = get(kw, :url, nothing))
    end
    _readable(track, path) && return _providerecorded!(track, path, sidecar)
    known = source isa SourceFacts ? source : Provenancer._factsat(sidecar)
    isnothing(known) &&
        return error(_neitherrun(track, path, nothing,
                                 DepositSharesIdentifier()))
    return error(_neitherrun(track, path, known.namespace, known.deposit))
end

# What `fetchfile` does for a download given a run's track, as `_trackedrecord` does for a file of
# the user's own: in a fetch run the file is fetched as an untracked call fetches it and registered
# from its record's facts, which must name its source; in a model run the run's copy of the file and
# its record are put in place from the store; a run whose config does neither is refused.
function Provenancer._trackedfetchfile(track::Track{<:DataRegistryHandle},
                                       file::Provenancer.URLFile; role::Symbol,
                                       verify::Bool = false)
    path = file.path
    Provenancer._checkstems([path])
    facts = file.facts
    isnothing(facts) &&
        error("`$(basename(path))` has no facts of its source, so the registry cannot be told " *
              "where it comes from: give the `URLFile` its `facts`.")
    Provenancer._namessource(facts) || error(_nosource(facts))
    if _writable(track, path)
        _claimsource!(track, facts, path)
        _checkwritable(track, facts, [(path = path,)])
        r = Provenancer.fetchfile(file, role = role, verify = verify)
        return _registerrecorded!(track, r, Provenancer.route(file),
                                  url = file.url)
    end
    _readable(track, path) &&
        return _providerecorded!(track, path, Provenancer.sidecarpath(path))
    return error(_neitherrun(track, path, facts.namespace, facts.deposit))
end

# What `missingfiles` reports of a tracked file of the user's own, or of a download: what its run's
# working config makes of it, from the record beside it where there is one, else from the `facts`
# given, a download's, which is fetched where it has a `url` and given by hand where it has none. The citations are read first, as for a spec, since a registered file's facts
# are compared by its citation key where it has no DOI.
function Provenancer._trackedinputmissing(track::Track{<:DataRegistryHandle},
                                          path::AbstractString;
                                          route::EcoSISTEM.AbstractFileRoute,
                                          sidecar::AbstractString,
                                          facts::Union{Nothing, SourceFacts} = nothing,
                                          url::Union{Nothing, AbstractString} = nothing)
    known = something(Provenancer._factsat(sidecar), facts, Some(nothing))
    src = isnothing(known) ? nothing : (facts = known, fileversioned = false)
    entry = (path = String(path), route = route, code = nothing, url = nothing)
    item = _item(entry,
                 EcoSISTEM._trackrelease!(Dict{String, Any}(), track, route),
                 sidecar = sidecar)
    _readcitations!(track.handle)
    todo = @NamedTuple{file::String, reason::String}[]
    reason = _missingreason(track, src, item, fetched = !isnothing(url))
    isnothing(reason) || push!(todo, (file = String(path), reason = reason))
    return todo
end
