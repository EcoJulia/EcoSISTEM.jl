# SPDX-License-Identifier: LGPL-3.0-or-later

# A data pipeline run as the route a spec's files take: the run, and the name they are filed under.

# A data product name, as DataPipeline.jl checks one: ASCII letters, digits, `_`, `-` and `.`, in
# parts joined by `/`, no part beginning or ending with `.`.
const _PRODUCTNAME = r"^[A-Za-z0-9_-](?:[A-Za-z0-9_.-]*[A-Za-z0-9_-])?(?:/[A-Za-z0-9_-](?:[A-Za-z0-9_.-]*[A-Za-z0-9_-])?)*$"

"""
    Track(handle, root::AbstractString)

Name a FAIR data pipeline run as the route a spec's files take, given as the `track` keyword of
[`SourceSpec`](@ref) and [`NaturalEarthSpec`](@ref): each file is the data product `<root>/<file
name>`, which [`fetchfiles`](@ref) registers in a fetch run, and which any other read takes from the
run's store, never fetching. A `handle` of `nothing` tracks nothing, so one script runs with the
pipeline or without it by the handle alone.

# Fields

  - `handle`: the run, as `DataPipeline.initialise()` returns it, or `nothing` for none.
  - `root`: the leading part of every data product name the spec's files are filed under, as
    `"era5"` or `"inputs/era5"` - ASCII letters, digits, `_`, `-` and `.`, in parts joined by `/`,
    no part beginning or ending with `.`. Checked whatever the handle, so a script refuses a bad
    root with the pipeline or without it.
"""
struct Track{H}
    handle::H
    root::String

    function Track(handle::H, root::AbstractString) where {H}
        occursin(_PRODUCTNAME, root) ||
            error("`$root` cannot begin a data product name: a root holds ASCII letters, digits, " *
                  "`_`, `-` and `.`, in parts joined by `/`, no part beginning or ending with `.`.")
        return new{H}(handle, String(root))
    end
end

# The handle is the run itself, which has nothing short to say, so only whether there is one shows.
function Base.show(io::IO, track::Track)
    run = isnothing(track.handle) ? "nothing" : "..."
    return print(io, "Track(", run, ", ", repr(track.root), ")")
end

# What a spec stores of its `track` keyword: the track, or `nothing` for no tracking, whether that was
# given as `nothing` or as a track with no run.
_storedtrack(::Nothing) = nothing

_storedtrack(::Track{Nothing}) = nothing

_storedtrack(track::Track) = track
