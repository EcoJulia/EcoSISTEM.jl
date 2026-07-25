# SPDX-License-Identifier: LGPL-3.0-or-later

using Feather
using DataFrames
using Missings

mutable struct DiversitySet
    data::Union{DataFrame, Missing}
    folder::String
    times::Vector{Unitful.Time}
end

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

Return the timepoints in a [`DiversitySet`](@ref) for which diversity has not
yet been calculated. If a previously saved Feather file of results is found,
only times beyond the latest recorded time are returned.
"""
function gettimes(diversityset::DiversitySet)
    file = searchdir(diversityset.folder, ".feather")
    if ismissing(diversityset.data) & isempty(file)
        return diversityset.times
    else
        if ismissing(diversityset.data) & !isempty(file)
            diversityset.data = Feather.read(joinpath(diversityset.folder,
                                                      file[1]))
            diversityset.data[:type_name] = ""
            diversityset.data[:time] *= 1s
        end
        latesttime = maximum(diversityset.data[:time])
        newtimes = diversityset.times[diversityset.times .> latesttime]
        return newtimes
    end
end

import DataFrames.append!
"""
    append!(diversityset::DiversitySet, dat::DataFrame)

Append a `DataFrame` of diversity results `dat` to the data stored in a
[`DiversitySet`](@ref).
"""
function append!(diversityset::DiversitySet, dat::DataFrame)
    return append!(diversityset.data, dat)
end
