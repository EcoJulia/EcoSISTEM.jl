using Unitful
using Missings
using JLD2
using Random

function checkfile(::String, ::Missing)
    return false
end

function checkfile(file::String, tm::Int)
    return !isempty(searchdir(file, string(tm, ".jld2")))
end

function loadfile(file::String, tm::Int, dim::Tuple)
    a = @load joinpath(file, searchdir(file, string(tm, ".jld2"))[1]) abuns
    return GridLandscape(a, dim)
end

function clearcache(cache::CachedEcosystem)
    files = searchdir(cache.abundances.outputfolder, ".jld2")
    rm.(joinpath.(cache.abundances.outputfolder, files))
    len = length(files)
    return  "$len files cleared"
end


function _abundances(cache::CachedEcosystem, tm::Unitful.Time)
    yr = mod(tm, 1year) == 0year ? Int(ustrip(uconvert(year, tm))) : missing
    if ismissing(cache.abundances.matrix[tm])
        if checkfile(cache.abundances.outputfolder, yr)
            cache.abundances.matrix[tm] = loadfile(cache.abundances.outputfolder,
                                                    yr, (length(cache.spplist.names),
                                                    _getdimension(cache.abenv.habitat) ...))
            seed!(cache.abundances.matrix[tm].rngs)
            return tm, cache.abundances.matrix[tm]
        else
            newtm, abun =  _abundances(cache, tm - cache.abundances.saveinterval)
            if (newtm > 2 * cache.abundances.saveinterval)
                cache.abundances.matrix[(newtm - 2 * cache.abundances.saveinterval)] = missing
            end
        end
    else
        return tm, cache.abundances.matrix[tm]
    end
    simulate!(cache, newtm, cache.abundances.saveinterval)
    if !ismissing(yr)
        @save joinpath(cache.abundances.outputfolder, string(yr, ".jld2")) abuns = SavedLandscape(cache.abundances.matrix[tm])
    end
    _abundances(cache, newtm + cache.abundances.saveinterval)

end

"""
    abundances(cache::CachedEcosystem, tm::Unitful.Time)

Function to extract abundances for an ecosystem, `cache`, at a certain point
in time, `tm`. If the abundances for that time are missing from the ecosystem,
then the function checks on disk for the last saved version and simulates
forward.
"""
function abundances(cache::CachedEcosystem, tm::Unitful.Time)
    return _abundances(cache, tm)[2]
end
