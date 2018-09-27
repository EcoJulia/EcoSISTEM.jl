using Unitful
using Missings
function abundances(cache::CachedEcosystem, tm::Unitful.Time)
    if ismissing(cache.abundances.matrix[tm])
        if checkfile(cache.abundances.outputfolder, tm)
            cache.abundances.matrix[tm] = loadfile(cache.abundances.outputfile,
                                                    tm)
            return tm
        else
            newtm =  abundances(cache, tm - cache.abundances.saveinterval)
        end
    else
        return tm
    end
    difftm = tm - newtm
    if difftm == 0s
        return cache.abundances.matrix[tm]
    else
        simulate!(cache, difftm, cache.abundances.saveinterval)
        return cache.abundances.matrix[tm]
    end
end


searchdir(path,key) = filter(x->contains(x,key), readdir(path))

function checkfile(file::String, tm::Unitful.Time)
    checktm = ustrip(tm)
    return !isempty(searchdir(file, string(checktm)))
end

function loadfile(file::String, tm::Unitful.Time)
    checktm = ustrip(tm)
    return load(searchdir(file, string(tm)), string(tm))
end
