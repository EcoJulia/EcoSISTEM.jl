using Unitful
using Missings
using JLD

searchdir(path,key) = filter(x->contains(x,key), readdir(path))

function checkfile(file::String, tm::Int)
    return !isempty(searchdir(file, string(tm, ".jld")))
end

function checkfile(::String, ::Missing)
    return false
end

function loadfile(file::String, tm::Int)
    return load(searchdir(file, string(tm, ".jld"))[1], string(tm))
end

function _abundances(cache::CachedEcosystem, tm::Unitful.Time)
    yr = mod(tm, 1year) == 0year ? Int(ustrip(uconvert(year, tm))) : missing
    if ismissing(cache.abundances.matrix[tm])
        if checkfile(cache.abundances.outputfolder, yr)
            cache.abundances.matrix[tm] = loadfile(cache.abundances.outputfolder,
                                                    yr)
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
        save(string(yr, ".jld"),
        string(yr), cache.abundances.matrix[tm])
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
