"""
    EpiLandscape

Disease class abundances housed in the landscape. These are represented in both 2 dimensions (for computational efficiency in simulations) and 3 dimensions (to represent disease classes, their abundances and position in the grid).

"""
mutable struct EpiLandscape
  matrix::Matrix{Int64}
  grid::Array{Int64, 3}
  seed::Vector{MersenneTwister}

  function EpiLandscape(abun::Matrix{Int64}, dimension::Tuple)
    a = abun
    return new(a, reshape(a, dimension), [MersenneTwister(rand(UInt)) for _ in 1:Threads.nthreads()])
  end
  function EpiLandscape(abun::Matrix{Int64}, dimension::Tuple, seed::Vector{MersenneTwister})
    a = abun
    return new(a, reshape(a, dimension), seed)
  end
end
import Base.copy
function copy(gl::EpiLandscape)
    return EpiLandscape(copy(gl.matrix), size(gl.grid))
end

"""
    emptyepilandscape(epienv::GridEpiEnv, epilist::EpiList)

Function to create an empty EpiLandscape given a GridEpiEnv and a
EpiList.
"""
function emptyepilandscape(epienv::GridEpiEnv, epilist::EpiList)
  mat = zeros(Int64, counttypes(epilist, true), countsubcommunities(epienv))

  dimension = (counttypes(epilist, true), _getdimension(epienv.habitat)...)
  return EpiLandscape(mat, dimension)
end
