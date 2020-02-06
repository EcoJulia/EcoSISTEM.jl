using Missings
using AxisArrays

"""
    MPIGridLandscape

Ecosystem abundances housed in the landscape. These are represented in both 2
dimensions (for computational efficiency in simulations) and 3 dimensions (to
represent species, their abundances and position in the grid).

"""
mutable struct MPIGridLandscape
  matrix::Matrix{Int64}
  grid::Array{Int64, 3}
  seed::Vector{MersenneTwister}

  function MPIGridLandscape(abun::Matrix{Int64}, dimension::Tuple)
    a = abun
    return new(a, reshape(a, dimension), [MersenneTwister(rand(UInt)) for _ in 1:Threads.nthreads()])
  end
  function MPIGridLandscape(abun::Matrix{Int64}, dimension::Tuple, seed::Vector{MersenneTwister})
    a = abun
    return new(a, reshape(a, dimension), seed)
  end
end
import Base.copy
function copy(gl::MPIGridLandscape)
    return MPIGridLandscape(copy(gl.matrix), size(gl.grid))
end


"""
    emptygridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)

Function to create an empty MPIGridLandscape given a GridAbioticEnv and a
SpeciesList.
"""
function emptyMPIgridlandscape(gae::GridAbioticEnv, spplist::SpeciesList)
  mat = zeros(Int64, counttypes(spplist, true), countsubcommunities(gae))

  dimension = (counttypes(spplist, true), _getdimension(gae.habitat)...)
  return MPIGridLandscape(mat, dimension)
end
