# SPDX-License-Identifier: LGPL-3.0-or-later

using Missings
using AxisArrays
using Random

struct SavedLandscape
    matrix::Matrix{Int64}
    rngs::Vector{Random.Xoshiro}
end

"""
    GridLandscape

Ecosystem abundances housed in the landscape. `matrix` stores abundances as a
2-dimensional array (species × grid cells) for computational efficiency, and
`grid` is a 3-dimensional view of the same data (species × x × y). Random draws
during simulation use Julia's task-local default RNG, so no generator state is
stored here.
"""
mutable struct GridLandscape
    matrix::Matrix{Int64}
    grid::Array{Int64, 3}

    function GridLandscape(abun::Matrix{Int64}, dimension::Tuple)
        a = abun
        return new(a, reshape(a, dimension))
    end
end
import Base.copy
function copy(gl::GridLandscape)
    return GridLandscape(copy(gl.matrix), size(gl.grid))
end
"""
    GridLandscape(sl::SavedLandscape, dimension::Tuple)

Restore a `GridLandscape` from a [`SavedLandscape`](@ref), reshaping the
abundance matrix to `dimension`. The saved RNG snapshot is restored separately
at the load site (see [`loadfile`](@ref)).
"""
function GridLandscape(sl::SavedLandscape, dimension::Tuple)
    return GridLandscape(sl.matrix, dimension)
end

"""
    SavedLandscape(gl::GridLandscape, rngs::Vector{Random.Xoshiro})

Convert a [`GridLandscape`](@ref) to a `SavedLandscape` for serialisation,
preserving the abundance matrix and a snapshot of the per-species RNG streams
`rngs` so a cached run can be resumed with a reproducible random stream.
"""
function SavedLandscape(gl::GridLandscape, rngs::Vector{Random.Xoshiro})
    return SavedLandscape(gl.matrix, copy.(rngs))
end

"""
    CachedGridLandscape

Ecosystem abundances for a cached simulation. `matrix` is an `AxisArray` over
time where each slot holds either a [`GridLandscape`](@ref) or `missing`.
`outputfolder` is the path to the folder where JLD2 cache files are written.
`timestep` is the simulation step (the granularity of the time axis) and
`saveinterval` is the (possibly coarser) interval at which checkpoints are
written to disk; it must be a multiple of `timestep`. Because the simulation
always advances by `timestep`, the results are independent of `saveinterval`.
"""
mutable struct CachedGridLandscape
    matrix::AxisArray{Union{GridLandscape, Missing}, 1}
    outputfolder::String
    saveinterval::Unitful.Time
    timestep::Unitful.Time
end

"""
    CachedGridLandscape(file::String, times::StepRangeLen;
                        saveinterval::Unitful.Time = step(times))

Construct a `CachedGridLandscape` backed by the folder `file`, initialising all
timepoints in the range `times` to `missing`. The simulation timestep is the
step size of `times`; `saveinterval` sets how often checkpoints are written to
disk and must be a multiple of the timestep (it defaults to saving every step).
"""
function CachedGridLandscape(file::String, times::StepRangeLen;
                             saveinterval::Unitful.Time = step(times))
    timestep = step(times)
    iszero(mod(saveinterval, timestep)) ||
        error("saveinterval ($saveinterval) must be a multiple of the timestep ($timestep)")
    v = Vector{Union{GridLandscape, Missing}}(undef, length(times))
    fill!(v, missing)
    a = AxisArray(v, Axis{:time}(times))
    return CachedGridLandscape(a, file, saveinterval, timestep)
end

"""
    emptygridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)

Create an empty [`GridLandscape`](@ref) given a [`GridAbioticEnv`](@ref) and a
[`SpeciesList`](@ref).
"""
function emptygridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)
    mat = zeros(Int64, counttypes(spplist, true), countsubcommunities(gae))

    dimension = (counttypes(spplist, true), _getdimension(gae.habitat)...)
    return GridLandscape(mat, dimension)
end
