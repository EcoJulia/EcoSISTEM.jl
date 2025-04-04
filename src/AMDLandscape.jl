# SPDX-License-Identifier: LGPL-3.0-or-later

using Missings
using AxisArrays
using AMDGPU
using Random

"""
    AMDGridLandscape

Ecosystem abundances housed in the landscape. These are represented in both 2
dimensions (for computational efficiency in simulations) and 3 dimensions (to
represent species, their abundances and position in the grid).

"""
mutable struct AMDGridLandscape
    matrix::Matrix{Int64}
    grid::Array{Int64, 3}
    rocmatrix::ROCArray{Int64, 2, AMDGPU.Runtime.Mem.HIPBuffer}
    rngs::Vector{MersenneTwister}

    function AMDGridLandscape(abun::Matrix{Int64}, dimension::Tuple)
        a = abun
        return new(a, reshape(a, dimension), ROCArray(a),
                   [MersenneTwister(rand(UInt)) for _ in 1:Threads.nthreads()])
    end
    function AMDGridLandscape(abun::Matrix{Int64}, dimension::Tuple,
                           rngs::Vector{MersenneTwister})
        a = abun
        return new(a, reshape(a, dimension), ROCArray(a), rngs)
    end
end
import Base.copy
function copy(gl::AMDGridLandscape)
    return AMDGridLandscape(copy(gl.matrix), size(gl.grid))
end

"""
    emptyAMDgridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)

Function to create an empty AMDGridLandscape given a GridAbioticEnv and a
SpeciesList.
"""
function emptyAMDgridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)
    mat = zeros(Int64, counttypes(spplist, true), countsubcommunities(gae))

    dimension = (counttypes(spplist, true), _getdimension(gae.habitat)...)
    return AMDGridLandscape(mat, dimension)
end

function synchronise_to_gpu!(agl::AMDGridLandscape)
    copyto!(agl.rocmatrix, agl.matrix)
end

function synchronise_from_gpu!(agl::AMDGridLandscape)
    copyto!(agl.matrix, agl.rocmatrix)
end
