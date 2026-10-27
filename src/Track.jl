# SPDX-License-Identifier: LGPL-3.0-or-later

# A data pipeline run as the route a spec's files take: the run, and the name they are filed under.

# A data product name, as DataPipeline.jl checks one: ASCII letters, digits, `_`, `-` and `.`, in
# parts joined by `/`, no part beginning or ending with `.`.
const _PRODUCTNAME = r"^[A-Za-z0-9_-](?:[A-Za-z0-9_.-]*[A-Za-z0-9_-])?(?:/[A-Za-z0-9_-](?:[A-Za-z0-9_.-]*[A-Za-z0-9_-])?)*$"

"""
    Track(handle, root::AbstractString; release_version = nothing, release_date = nothing)

Name a FAIR data pipeline run as the route a spec's files take, given as the `track` keyword of
[`SourceSpec`](@ref) and [`NaturalEarthSpec`](@ref): each file is the data product
`<root>/data/<stem>`, the stem its name without its extension, and its provenance record
`<root>/provenance/<stem>`, which [`fetchfiles`](@ref) registers in a fetch run, and which any
other read takes from the run's store, never fetching. A `handle` of `nothing` tracks nothing, so
one script runs with the pipeline or without it by the handle alone.

# Fields

  - `handle`: the run, as `DataPipeline.initialise()` returns it, or `nothing` for none.
  - `root`: the leading part of every data product name the spec's files are filed under, as
    `"era5"` or `"inputs/era5"` - ASCII letters, digits, `_`, `-` and `.`, in parts joined by `/`,
    no part beginning or ending with `.`. Checked whatever the handle, so a script refuses a bad
    root with the pipeline or without it.
  - `release_version`, `release_date`: the release the files you give a spec come from, where
    its catalogue row cannot say - CRU TS's writes `4`, the series, where your files are of one
    release - as a version, numbers with no leading zero as the registry keeps them (`"4.8"` or
    `v"4.8"` for CRU TS 4.08), and a date (a `Date`, or `"2024"`, `"2024-05"` or
    `"2024-05-30"`); `nothing` for the row's version and the date its citation gives. Only for
    files given by `file`, `files` or `directory`: a spec whose files are fetched for it refuses
    them, its row and its provider saying which release they are.
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

# What a spec stores of its `track` keyword: the track, or `nothing` for no tracking, whether that was
# given as `nothing` or as a track with no run.
_storedtrack(::Nothing) = nothing

_storedtrack(::Track{Nothing}) = nothing

_storedtrack(track::Track) = track
