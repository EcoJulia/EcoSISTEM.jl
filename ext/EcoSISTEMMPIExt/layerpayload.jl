# SPDX-License-Identifier: LGPL-3.0-or-later
#
# One rank's build handed to the others. The root builds and broadcasts the result - a small
# skeleton by `MPI.bcast`, then each array in place by `MPI.Bcast!` - and every other rank rebuilds
# it from those, so the data is read once however many ranks there are. Every rank must make each
# call, in the same order, or the ones that did wait forever.

using Base.ScopedValues: ScopedValue, with

# The rank that builds.
const _BUILDROOT = 0

# Whether this rank is already inside a shared build. A nested one would broadcast from the root
# alone and leave every other rank waiting, so it is refused.
const _INSHAREDBUILD = ScopedValue(false)

# One layer, built on `area`, whose cache records what the root read.
function EcoSISTEM._sharedlayer(build, area::EcoSISTEM.StudyArea)
    return _sharedbuild(build, EcoSISTEM._layerpayload,
                        area.report.cache) do skeleton, arrays
        return EcoSISTEM._rebuildlayer(skeleton, arrays)
    end
end

# One study area's report, rebuilt on each rank around that rank's own specs, constraints and cache.
function EcoSISTEM._sharedreport(build, layers::NamedTuple, cons::NamedTuple,
                                 cache::EcoSISTEM.LayerCache)
    return _sharedbuild(build, EcoSISTEM._reportpayload,
                        cache) do skeleton, arrays
        return EcoSISTEM._rebuildreport(skeleton, arrays, layers, cons, cache)
    end
end

# The root runs `build` and takes its result apart with `pack`; every other rank receives the parts
# and puts them together with `rebuild`. Then the root's input records replace everyone else's in
# `cache`.
function _sharedbuild(rebuild, build, pack, cache)
    _refusenesting()
    comm = MPI.COMM_WORLD
    isroot = MPI.Comm_rank(comm) == _BUILDROOT
    outcome = isroot ? _rootbuild(build, pack) : nothing
    header = MPI.bcast(isroot ? outcome.header : nothing, _BUILDROOT, comm)
    header.ok || _buildfailed(header, outcome)
    arrays = isroot ? outcome.payload.arrays :
             map(EcoSISTEM._allocatepayload, header.descriptors)
    for array in arrays
        MPI.Bcast!(array, _BUILDROOT, comm)
    end
    _shareinputs!(cache, comm)
    return isroot ? outcome.value : rebuild(header.skeleton, arrays)
end

# Raised before anything is broadcast, so that on the root it reaches the enclosing build, which
# hands it to every rank.
function _refusenesting()
    _INSHAREDBUILD[] || return nothing
    return error("a build shared between MPI ranks was started inside another one, which would " *
                 "leave every other rank waiting.")
end

# The root's build and its payload, or the exception that stopped either, caught so that it can be
# broadcast before it is raised.
function _rootbuild(build, pack)
    try
        value = with(build, _INSHAREDBUILD => true)
        payload = pack(value)
        return (value = value, payload = payload,
                header = (ok = true, skeleton = payload.skeleton,
                          descriptors = payload.descriptors))
    catch err
        return (exception = err,
                header = (ok = false, message = sprint(showerror, err)))
    end
end

# A failed build stops every rank: the root with its own exception, the rest naming it.
function _buildfailed(header::NamedTuple, ::Nothing)
    return error("the root MPI rank failed to build what every rank was waiting for:\n" *
                 header.message)
end

_buildfailed(header::NamedTuple, outcome::NamedTuple) = throw(outcome.exception)

# The root's input records replace every other rank's, since only the root read anything. An area
# whose cache has been dropped records none.
_shareinputs!(::Nothing, comm) = nothing

function _shareinputs!(cache::EcoSISTEM.LayerCache, comm)
    isroot = MPI.Comm_rank(comm) == _BUILDROOT
    inputs = MPI.bcast(isroot ? cache.inputs : nothing, _BUILDROOT, comm)
    isroot && return cache
    empty!(cache.inputs)
    merge!(cache.inputs, inputs)
    return cache
end
