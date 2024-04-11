using Feather
using DataFrames
using Missings

mutable struct DiversitySet
    data::Union{DataFrame, Missing}
    folder::String
    times::Vector{Unitful.Time}
end

function DiversitySet(cache::CachedEcosystem,
                      times::Vector{T}) where {T <: Unitful.Time}
    return DiversitySet(missing, cache.abundances.outputfolder, times)
end

function updatesimulation!(cache::CachedEcosystem, tm::Unitful.Time)
    return abundances(cache, tm)
end

function gettimes(div::DiversitySet)
    file = searchdir(div.folder, ".feather")
    if ismissing(div.data) & isempty(file)
        return div.times
    else
        if ismissing(div.data) & !isempty(file)
            div.data = Feather.read(joinpath(div.folder, file[1]))
            div.data[:type_name] = ""
            div.data[:time] *= 1s
        end
        latesttime = maximum(div.data[:time])
        newtimes = div.times[div.times .> latesttime]
        return newtimes
    end
end

import DataFrames.append!
function append!(div::DiversitySet, dat::DataFrame)
    return append!(div.data, dat)
end
