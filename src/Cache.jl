# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Missings
using JLD2
using Random

"""
    checkfile(::String, ::Missing)

Check whether a cache file exists for a given timepoint. Always returns `false`
when the timepoint is `missing`.
"""
function checkfile(::String, ::Missing)
    return false
end

"""
    checkfile(file::String, tm::Int)

Check whether a JLD2 cache file exists in the folder `file` for timestep `tm`.
Returns `true` if a matching file is found.
"""
function checkfile(file::String, tm::Int)
    return !isempty(searchdir(file, string(tm, ".jld2")))
end

"""
    loadfile(file::String, tm::Int, dim::Tuple)

Load a cached [`GridLandscape`](@ref) from the folder `file` for timestep `tm`,
reshaping the abundance matrix to dimensions `dim`.
"""
function loadfile(file::String, tm::Int, dim::Tuple)
    @load joinpath(file, searchdir(file, string(tm, ".jld2"))[1]) abuns
    # Restore the task-local RNG so the resumed run continues a reproducible stream
    copy!(Random.default_rng(), abuns.rng)
    return GridLandscape(abuns, dim)
end

"""
    clearcache(cache::CachedEcosystem)

Delete all JLD2 cache files from the output folder of a
[`CachedEcosystem`](@ref). Returns a string reporting how many files were
removed.
"""
function clearcache(cache::CachedEcosystem)
    files = searchdir(cache.abundances.outputfolder, ".jld2")
    rm.(joinpath.(cache.abundances.outputfolder, files))
    len = length(files)
    return "$len files cleared"
end

function _abundances(cache::CachedEcosystem, tm::Unitful.Time)
    yr = mod(tm, 1year) == 0year ? Int(ustrip(uconvert(year, tm))) : missing
    if ismissing(cache.abundances.matrix[tm])
        if checkfile(cache.abundances.outputfolder, yr)
            cache.abundances.matrix[tm] = loadfile(cache.abundances.outputfolder,
                                                   yr,
                                                   (length(cache.spplist.names),
                                                    _getdimension(cache.abenv.habitat)...))
            return tm, cache.abundances.matrix[tm]
        else
            newtm, abun = _abundances(cache, tm - cache.abundances.saveinterval)
            if (newtm > 2 * cache.abundances.saveinterval)
                cache.abundances.matrix[(newtm - 2 * cache.abundances.saveinterval)] = missing
            end
        end
    else
        return tm, cache.abundances.matrix[tm]
    end
    simulate!(cache, newtm, cache.abundances.saveinterval)
    if !ismissing(yr)
        @save joinpath(cache.abundances.outputfolder, string(yr, ".jld2")) abuns=SavedLandscape(cache.abundances.matrix[tm])
    end
    return _abundances(cache, newtm + cache.abundances.saveinterval)
end

"""
    abundances(cache::CachedEcosystem, tm::Unitful.Time)

Extract abundances for an ecosystem, `cache`, at a certain point in time, `tm`.
If the abundances for that time are missing from the ecosystem, then the
function checks on disk for the last saved version and simulates forward.
"""
function abundances(cache::CachedEcosystem, tm::Unitful.Time)
    return _abundances(cache, tm)[2]
end
