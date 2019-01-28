using Feather
using DataFrames
using Missings

mutable struct DiversitySet
    data::Union{DataFrame, Missing}
    folder::String
    times::Vector{Unitful.Time}
end

function DiversitySet(cache::CachedEcosystem, times::Vector{T}) where T <: Unitful.Time
    return DiversitySet(missing, cache.abundances.outputfolder, times)
end
GLOBAL_typedict["DiversitySet"] = DiversitySet

function updatesimulation!(cache::CachedEcosystem, tm::Unitful.Time)
    abundances(cache, tm)
end
GLOBAL_funcdict["updatesimulation!"] = updatesimulation!
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
GLOBAL_funcdict["gettimes"] = gettimes
import DataFrames.append!
function append!(div::DiversitySet, dat::DataFrame)
    append!(div.data, dat)
end
