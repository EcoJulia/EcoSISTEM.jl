# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A recorded run's diversity measurements, alongside the ecosystem that produced them.

using DataFrames
using Feather
using Missings
using Unitful

"""
    DiversitySet

Diversity measurements recorded from a [`CachedEcosystem`](@ref) over a run, held alongside the
folder its cached abundances were written to.

# Fields

  - `data`: the measurements so far, or `missing` before any have been taken.
  - `folder`: where the cached abundances live, which is also where results are read back from.
  - `times`: every timepoint diversity is wanted at. [`gettimes`](@ref) narrows this to the ones not
    yet computed.
"""
mutable struct DiversitySet
    data::Union{DataFrame, Missing}
    folder::String
    times::Vector{Unitful.Time}
end
# == Functions ==================================================================================

"""
    DiversitySet(cache::CachedEcosystem, times::Vector{T}) where T <: Unitful.Time

Construct a `DiversitySet` from a [`CachedEcosystem`](@ref), initialising it
with the cache output folder and a vector of timepoints `times` for which
diversity is to be recorded.
"""
function DiversitySet(cache::CachedEcosystem,
                      times::Vector{T}) where {T <: Unitful.Time}
    return DiversitySet(missing, cache.abundances.outputfolder, times)
end

"""
    updatesimulation!(cache::CachedEcosystem, tm::Unitful.Time)

Trigger the computation and caching of ecosystem abundances at timepoint `tm` in
a [`CachedEcosystem`](@ref).
"""
function updatesimulation!(cache::CachedEcosystem, tm::Unitful.Time)
    return abundances(cache, tm)
end

"""
    gettimes(diversityset::DiversitySet)

Return the timepoints in a [`DiversitySet`](@ref) for which diversity has not yet been calculated.
Where a previously saved Feather file of results is found, only times beyond the latest it records
are returned, so a run can be resumed rather than repeated.

# Arguments

  - `diversityset`: the set to ask, whose `folder` is searched for saved results.
"""
function gettimes(diversityset::DiversitySet)
    file = _searchdir(diversityset.folder, ".feather")
    (ismissing(diversityset.data) && isempty(file)) && return diversityset.times
    if ismissing(diversityset.data)
        # A saved file's `time` column is stored bare, so it is given back its unit on the way in and
        # every comparison below is between `Unitful.Time`s.
        loaded = DataFrame(Feather.read(joinpath(diversityset.folder,
                                                 first(file))))
        loaded[!, :type_name] .= ""
        loaded[!, :time] = loaded[!, :time] .* 1s
        diversityset.data = loaded
    end
    latest = maximum(diversityset.data[!, :time])
    return diversityset.times[diversityset.times .> latest]
end
