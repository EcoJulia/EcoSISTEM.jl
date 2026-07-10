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
    checkfile(file::String, idx::Int)

Check whether a JLD2 checkpoint file exists in the folder `file` for checkpoint
index `idx`. Returns `true` if `<idx>.jld2` is present.
"""
function checkfile(file::String, idx::Int)
    return isfile(joinpath(file, string(idx, ".jld2")))
end

"""
    loadfile(cache::CachedEcosystem, file::String, idx::Int, dim::Tuple)

Load a cached [`GridLandscape`](@ref) from the folder `file` for checkpoint index
`idx`, reshaping the abundance matrix to dimensions `dim`. The saved per-species
RNG streams are restored into `cache.rngs` so the resumed run continues a
reproducible random stream.
"""
function loadfile(cache::CachedEcosystem, file::String, idx::Int, dim::Tuple)
    @load joinpath(file, string(idx, ".jld2")) abuns
    # Restore the per-species RNG streams for a reproducible resumed run
    cache.rngs .= copy.(abuns.rngs)
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
    timestep = cache.abundances.timestep
    saveinterval = cache.abundances.saveinterval
    # Checkpoints are written to disk only at multiples of `saveinterval`,
    # indexed by how many save intervals have elapsed; the simulation itself
    # always advances by `timestep`, so results are independent of `saveinterval`.
    issave = iszero(mod(tm, saveinterval))
    idx = issave ? Int(round(uconvert(NoUnits, tm / saveinterval))) : missing
    if ismissing(cache.abundances.matrix[tm])
        if checkfile(cache.abundances.outputfolder, idx)
            cache.abundances.matrix[tm] = loadfile(cache,
                                                   cache.abundances.outputfolder,
                                                   idx,
                                                   (length(cache.spplist.names),
                                                    _getdimension(cache.abenv.habitat)...))
            return tm, cache.abundances.matrix[tm]
        else
            newtm, abun = _abundances(cache, tm - timestep)
            if (newtm > 2 * timestep)
                cache.abundances.matrix[(newtm - 2 * timestep)] = missing
            end
        end
    else
        return tm, cache.abundances.matrix[tm]
    end
    simulate!(cache, newtm, timestep)
    if issave
        @save joinpath(cache.abundances.outputfolder, string(idx, ".jld2")) abuns=SavedLandscape(cache.abundances.matrix[tm],
                                                                                                 cache.rngs)
    end
    return _abundances(cache, newtm + timestep)
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
