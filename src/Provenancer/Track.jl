# SPDX-License-Identifier: LGPL-3.0-or-later

# A data pipeline run as the route files take: the run, and the name they are filed under.

# A data product name, as DataPipeline.jl checks one: ASCII letters, digits, `_`, `-` and `.`, in
# parts joined by `/`, no part beginning or ending with `.`.
const _PRODUCTNAME = r"^[A-Za-z0-9_-](?:[A-Za-z0-9_.-]*[A-Za-z0-9_-])?(?:/[A-Za-z0-9_-](?:[A-Za-z0-9_.-]*[A-Za-z0-9_-])?)*$"

"""
    Track(handle, root::AbstractString; release_version = nothing, release_date = nothing)

Name a FAIR data pipeline run as the route a set of files takes, given as the `track` keyword of
whatever reads them: each file is the data product `<root>/data/<stem>`, the stem its name without
its extension, and its provenance record `<root>/provenance/<stem>`, which a fetch run registers,
and which any other read takes from the run's store, never fetching. A `handle` of `nothing`
tracks nothing, so one script runs with the pipeline or without it by the handle alone.

# Fields

  - `handle`: the run, as `DataPipeline.initialise()` returns it, or `nothing` for none.
  - `root`: the leading part of every data product name the files are filed under, as `"era5"`
    or `"inputs/era5"` - ASCII letters, digits, `_`, `-` and `.`, in parts joined by `/`, no part
    beginning or ending with `.`. Checked whatever the handle, so a script refuses a bad root with
    the pipeline or without it.
  - `release_version`, `release_date`: the release the files the user gives come from, where
    their source's own facts cannot say - a series whose files are of one release - as a version,
    numbers with no leading zero as the registry keeps them (`"4.8"` or `v"4.8"` for a release
    called 4.08), and a date (a `Date`, or `"2024"`, `"2024-05"` or `"2024-05-30"`); `nothing`
    for the source's own version and the date its citation gives. Only for files the user gives:
    a file fetched for them carries its own release.
"""
struct Track{H}
    handle::H
    root::String
    release_version::Union{Nothing, String}
    release_date::Union{Nothing, String}

    function Track(handle::H, root::AbstractString;
                   release_version::Union{Nothing, AbstractString,
                                          VersionNumber} = nothing,
                   release_date::Union{Nothing, AbstractString, Dates.Date} = nothing) where {H}
        occursin(_PRODUCTNAME, root) ||
            error("`$root` cannot begin a data product name: a root holds ASCII letters, digits, " *
                  "`_`, `-` and `.`, in parts joined by `/`, no part beginning or ending with `.`.")
        return new{H}(handle, String(root),
                      _releaseversiontext(release_version),
                      _releasedatetext(release_date))
    end
end

# The handle is the run itself, which has nothing short to say, so only whether there is one shows.
function Base.show(io::IO, track::Track)
    run = isnothing(track.handle) ? "nothing" : "..."
    print(io, "Track(", run, ", ", repr(track.root))
    for (name, value) in (("release_version", track.release_version),
                          ("release_date", track.release_date))
        isnothing(value) || print(io, ", ", name, " = ", repr(value))
    end
    return print(io, ")")
end

# Whether a track follows a run at all: one with no handle tracks nothing, so whatever holds it
# takes the untracked route.
_istracked(::Nothing) = false

_istracked(track::Track) = !isnothing(track.handle)

# Why a track cannot be followed: its handle is not a pipeline run the package using this can
# use. The extension of that package supplies each tracked route's method for its handle.
function _nopipeline(track::Track)
    return "A `Track` follows a FAIR data pipeline run, the handle `DataPipeline.initialise()` " *
           "returns with DataPipeline loaded; a `$(typeof(track.handle))` is not one."
end

# A release version as a `Track` keeps it, text a registry takes as it stands, or `nothing`.
_releaseversiontext(::Nothing) = nothing

_releaseversiontext(v::VersionNumber) = string(v)

function _releaseversiontext(v::AbstractString)
    _isregistryversion(v) ||
        error("`release_version = $(repr(v))` is not a version a registry keeps as it is " *
              "written: give its numbers with no leading zero, `\"4.8\"` or `v\"4.8\"` for " *
              "a release called 4.08.")
    return String(v)
end

# Whether text is a version as a data pipeline registry keeps one: it reads as a version, and
# written in full - a missing minor or patch added - it reads back as it was written, which a
# number with a leading zero, as `08`, does not.
function _isregistryversion(v::AbstractString)
    parsed = tryparse(VersionNumber, v)
    isnothing(parsed) && return false
    numbers, rest = match(r"^([^-+]*)(.*)$", v)
    return any(string(parsed) == numbers * added * rest
               for added in ("", ".0", ".0.0"))
end

# A release date as a `Track` keeps it, ISO 8601 text to the precision given, or `nothing`.
_releasedatetext(::Nothing) = nothing

_releasedatetext(d::Dates.Date) = string(d)

function _releasedatetext(d::AbstractString)
    occursin(r"^\d{4}(-\d{2}(-\d{2})?)?$", d) ||
        error("`release_date = $(repr(d))` is not a date: write it `YYYY`, `YYYY-MM` or " *
              "`YYYY-MM-DD`, or give a `Date`.")
    return String(d)
end
