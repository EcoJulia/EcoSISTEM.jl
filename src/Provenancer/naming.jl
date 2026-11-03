# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a tracked file and its record are registered as with a FAIR data pipeline run, decided from
# plain data so that it is tested with no registry: the data product names under a track's root,
# the version a file takes whose bytes differ from those registered, the config lines a run lacks,
# why a file is still to do, and a record brought up to date with what is known of it now. The
# calls that ask a registry and register are the pipeline extension's of the package using this,
# which is why nothing here is public: the extension and its tests reach these by name.

# The stem a file's data products are named by: its name without its extension, which the data
# store adds.
_filestem(path::AbstractString) = first(splitext(basename(path)))

# The data product names a tracked file gives rise to under its track's root: the file itself
# (`<root>/data/<stem>`), its provenance record (`<root>/provenance/<stem>`), and the file-less
# deposit its source's files are derived from (`<root>/deposit`).
function _dataname(track::Track, path::AbstractString)
    return track.root * "/data/" * _filestem(path)
end

function _sidecarname(track::Track, path::AbstractString)
    return track.root * "/provenance/" * _filestem(path)
end

_depositname(track::Track) = track.root * "/deposit"

# Whether facts name their source as a registry takes one, by a DOI or by the key of a citation;
# files of a source naming neither cannot be registered with it.
_namessource(facts::SourceFacts) = !isempty(facts.doi) || !isempty(facts.key)

# Refuse, before anything is fetched, a file whose stem cannot end a data product name: one holding
# anything but ASCII letters, digits, `_`, `-` and `.`, or beginning or ending with `.`.
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

# A file's record of ours as it would be with `facts` (sidecar keys and values) in place of what it
# says of them, read without writing it; `nothing` where the file has no record of ours.
function _recordwith(path::AbstractString, facts::AbstractDict;
                     sidecar::AbstractString = sidecarpath(path))
    isfile(sidecar) || return nothing
    t = TOML.tryparsefile(sidecar)
    (t isa AbstractDict && _isourrecord(t)) || return nothing
    return _recordfromtoml(merge(t, facts), location = String(path),
                           sidecar = sidecar)
end

# Bring a file's record of ours up to date with `facts`, what is known of it now, before it is
# registered, rewriting it only where they differ: a record written by an earlier version may lack
# facts a registration needs, or hold ones since corrected, and the record and the registration
# must agree. What the fetch itself recorded - when, the size and checksum, the request, its job -
# is kept. Returns the record, or `nothing` for a file with no record of ours.
function _refreshrecord!(path::AbstractString, facts::AbstractDict;
                         sidecar::AbstractString = sidecarpath(path))
    isfile(sidecar) || return nothing
    t = TOML.tryparsefile(sidecar)
    (t isa AbstractDict && _isourrecord(t)) || return nothing
    if any(get(t, k, nothing) != v for (k, v) in facts)
        merge!(t, facts)
        t["writer"] = _writer()
        writesidecar(path, t, sidecar = sidecar)
    end
    return recordat(path, sidecar = sidecar)
end

# A release version as the registry writes it, in full: `2.1` as `2.1.0`, and `1.0.0` for none.
function _fullversion(version::AbstractString)
    return isempty(version) ? "1.0.0" :
           string(VersionNumber(version))
end

# The release version a file is registered under: the version its bytes are registered under
# already (`held`), so a third fetch of a changed file adds nothing; else its source's own
# (`version`, in full); or, where that is `taken` by other bytes - its provider changed the file
# without a new release - the same with build metadata `+YYYYMMDD` from `day`.
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
function _stampday(r::ProvenanceRecord, today::Dates.Date)
    return Dates.Date(something(r.modified, r.fetched, today))
end

# The file type of a file's data product, its extension without the dot, as a `write:` entry names
# it.
_filetype(path::AbstractString) = lstrip(last(splitext(path)), '.')

# What a fetch config needs for a track's files, as the lines to add: its source's namespace as the
# run's default for inputs and outputs, since the run registers that source's deposit and nothing
# else, and `write:` entries for the files, of `filetype`, and their provenance records. A source
# with no namespace of its own is registered in whatever the config's defaults are.
function _writelines(track::Track, namespace::AbstractString,
                     filetype::AbstractString)
    defaults = isempty(namespace) ? String[] :
               ["run_metadata:", "  default_input_namespace: $namespace",
        "  default_output_namespace: $namespace"]
    return join([defaults..., "write:",
                    "- data_product: $(track.root)/data/*",
                    "  description: A file of the source tracked under $(track.root)",
                    "  file_type: $filetype",
                    "- data_product: $(track.root)/provenance/*",
                    "  description: The provenance record of a file under $(track.root)",
                    "  file_type: toml"], "\n")
end

# The `read:` entry a model config needs for a track's files and their records, both in their
# source's namespace: `*` matches one segment of a name, so `<root>/*/*` takes the files and the
# records and leaves the deposit, which has no file. A source with no namespace of its own was
# registered in the namespace its fetch run's config gave as the default, which the entry's `use:`
# block must name.
function _readentry(track::Track, namespace::AbstractString)
    use = isempty(namespace) ?
          ["  use:", "    namespace: <the fetch run's default namespace>"] :
          ["  use:", "    namespace: $namespace"]
    return join(["read:", "- data_product: $(track.root)/*/*", use...], "\n")
end

# Why a fetch run still has a file to do, from what its registry answered - whether the file is
# `present`, the version its bytes are `held` under its name, whether that version records the
# `same` facts as known now, whether its record's bytes are `recorded`, and whether its name is
# registered (`named`) at all - or `nothing` where it has none.
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
        return "registered with other facts, a change in what is known of it: the fetch run " *
               "registers a new version"
    recorded ||
        return "its provenance record is not registered: the fetch run registers it"
    return nothing
end

# Why a model run cannot read a file, from whether its registry holds it: `nothing` where no
# `read:` entry covers its name, so the registry could not be asked.
function _modelreason(registered::Union{Nothing, Bool}, track::Track,
                      namespace::AbstractString)
    isnothing(registered) &&
        return "no `read:` entry covers it; add\n" *
               _readentry(track, namespace)
    registered &&
        return "registered since this job was prepared: `fair pull`, then `fair run` again"
    return "not registered in namespace `$namespace`: run the fetch run first"
end
