"""
    EpiLandscape

Disease class abundances housed in the landscape. These are represented in both 2 dimensions (for computational efficiency in simulations) and 3 dimensions (to represent disease classes, their abundances and position in the grid).

"""
mutable struct EpiLandscape{U <: Integer, VecRNGType <: AbstractVector{<:Random.AbstractRNG}}
  matrix::Matrix{U}
  matrix_v::Matrix{U}
  grid::Array{U, 3}
  seed::VecRNGType
end
function EpiLandscape(human_abun::Matrix{U}, virus_abun::Matrix{U}, d1::Tuple,
    rngtype::RNGType=Random.MersenneTwister) where {U <: Integer, RNGType}
  return EpiLandscape(human_abun, virus_abun, reshape(human_abun, d1),
    [rngtype(rand(UInt)) for _ in 1:Threads.nthreads()])
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
function emptyepilandscape(epienv::GridEpiEnv, epilist::EpiList, intnum::U,
    rngtype::RNGType=Random.MersenneTwister) where {U <: Integer, RNGType}
  mat_human = zeros(U, counttypes(epilist.human, true), countsubcommunities(epienv))
  mat_virus = zeros(U, counttypes(epilist.virus, true), countsubcommunities(epienv))
  dimension = (counttypes(epilist.human, true), _getdimension(epienv.habitat)...)
  return EpiLandscape(mat_human, mat_virus, dimension, rngtype)
end
