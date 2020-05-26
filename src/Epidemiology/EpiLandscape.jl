"""
    EpiLandscape

Disease class abundances housed in the landscape. These are represented in both 2 dimensions (for computational efficiency in simulations) and 3 dimensions (to represent disease classes, their abundances and position in the grid).

"""
mutable struct EpiLandscape
  matrix::Matrix{Int64}
  matrix_v::Matrix{Int64}
  grid::Array{Int64, 3}
  seed::Vector{MersenneTwister}
  function EpiLandscape(abun1::Matrix{Int64}, abun2::Matrix{Int64}, d1::Tuple)
    return new(abun1, abun2, reshape(abun1, d1), [MersenneTwister(rand(UInt)) for _ in 1:Threads.nthreads()])
  end
  function EpiLandscape(abun1::Matrix{Int64}, abun2::Matrix{Int64}, d1::Tuple, d2::Tuple, seed::Vector{MersenneTwister})
    return new(abun1, abun2, reshape(abun1, d1), seed)
  end
end

virus(x::EpiLandscape) = x.matrix_v
human(x::EpiLandscape) = x.matrix

import Base.copy
function copy(gl::EpiLandscape)
    return EpiLandscape(copy(gl.matrix), size(gl.grid))
end

function Base.isapprox(gl_1::EpiLandscape, gl_2::EpiLandscape; kwargs...)
    return isapprox(gl_1.matrix, gl_2.matrix; kwargs...)
end

"""
    emptyepilandscape(epienv::GridEpiEnv, epilist::EpiList)

Function to create an empty EpiLandscape given a GridEpiEnv and a
EpiList.
"""
function emptyepilandscape(epienv::GridEpiEnv, epilist::EpiList)
  mat = zeros(Int64, counttypes(epilist, true), countsubcommunities(epienv))
  dimension = (counttypes(epilist, true), _getdimension(epienv.habitat)...)
  return EpiLandscape(mat, mat, dimension)
end
