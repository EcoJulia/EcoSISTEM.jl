using Missings
using AxisArrays
using Random
using Compat

struct SavedLandscape
    matrix::Matrix{Int64}
    seed::Vector{MersenneTwister}
end

"""
    GridLandscape

Ecosystem abundances housed in the landscape. These are represented in both 2
dimensions (for computational efficiency in simulations) and 3 dimensions (to
represent species, their abundances and position in the grid).

"""
mutable struct GridLandscape
  matrix::Matrix{Int64}
  grid::Array{Int64, 3}
  seed::Vector{MersenneTwister}

  function GridLandscape(abun::Matrix{Int64}, dimension::Tuple)
    a = abun
    return new(a, reshape(a, dimension), [MersenneTwister(rand(UInt)) for _ in 1:Threads.nthreads()])
  end
  function GridLandscape(abun::Matrix{Int64}, dimension::Tuple, seed::Vector{MersenneTwister})
    a = abun
    return new(a, reshape(a, dimension), seed)
  end
end
import Base.copy
function copy(gl::GridLandscape)
    return GridLandscape(copy(gl.matrix), size(gl.grid))
end
function GridLandscape(sl::SavedLandscape, dimension::Tuple)
    GridLandscape(sl.matrix, dimension, sl.seed)
end

function SavedLandscape(gl::GridLandscape)
    SavedLandscape(gl.matrix, gl.seed)
end

"""
    CachedGridLandscape

Ecosystem abundances housed in the cached landscape. These are either stored in the matrix or output to a cache.

"""
mutable struct CachedGridLandscape
  matrix::AxisArray{Union{GridLandscape, Missing}, 1}
  outputfolder::String
  saveinterval::Unitful.Time
end

function CachedGridLandscape(file::String, rng::StepRangeLen)
  interval = step(rng)
  v = Vector{Union{GridLandscape, Missing}}(Compat.undef, length(rng))
  fill!(v, missing)
  a = AxisArray(v, Axis{:time}(rng))
  return CachedGridLandscape(a, file, interval)
end

"""
    emptygridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)

Function to create an empty GridLandscape given a GridAbioticEnv and a
SpeciesList.
"""
function emptygridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)
  mat = zeros(Int64, counttypes(spplist, true), countsubcommunities(gae))

  dimension = (counttypes(spplist, true), _getdimension(gae.habitat)...)
  return GridLandscape(mat, dimension)
end
