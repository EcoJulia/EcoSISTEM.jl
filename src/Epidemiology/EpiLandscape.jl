
"""
    EpiLandscape

Disease class abundances housed in the landscape. These are represented in both 2 dimensions (for computational efficiency in simulations) and 3 dimensions (to represent disease classes, their abundances and position in the grid).

"""
mutable struct EpiLandscape{U <: Integer, VecRNGType <: AbstractVector{<:Random.AbstractRNG}} <: AbstractLandscape
  matrix::Matrix{U}
  matrix_v::Matrix{U}
  grid::Array{U, 3}
  rngs::VecRNGType
end
function EpiLandscape(human_abun::Matrix{U}, virus_abun::Matrix{U}, d1::Tuple,
     Rngtype::Type{R} = Random.MersenneTwister
     ) where {U <: Integer, R <: Random.AbstractRNG}
  rngs = [Rngtype(rand(UInt)) for _ in 1:Threads.nthreads()]
  return EpiLandscape(human_abun, virus_abun, reshape(human_abun, d1), rngs)
end
function EpiLandscape(human_abun::Matrix{U}, virus_abun::Matrix{U}, d1::Tuple,
    d2::Tuple, rngs::VecRNGType) where {U <: Integer, VecRNGType <: AbstractVector{<:Random.AbstractRNG}}
  return EpiLandscape(human_abun, virus_abun, reshape(human_abun, d1), rngs)
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
    emptyepilandscape(epienv::GridEpiEnv, epilist::SpeciesList)

Function to create an empty EpiLandscape given a GridEpiEnv and a
SpeciesList.
"""
function emptyepilandscape(epienv::GridEpiEnv, epilist::SpeciesList, intnum::U,
    Rngtype::Type{R} = Random.MersenneTwister
    ) where {U <: Integer, R <: Random.AbstractRNG}
  mat_human = zeros(U, counttypes(epilist.species, true), countsubcommunities(epienv))
  mat_virus = zeros(U, counttypes(epilist.pathogens, true), countsubcommunities(epienv))
  dimension = (counttypes(epilist.species, true), _getdimension(epienv.habitat)...)
  return EpiLandscape(mat_human, mat_virus, dimension, Rngtype)
end
