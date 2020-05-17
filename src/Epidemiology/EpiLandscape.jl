abstract type AbstractAgent end # What to call this?

struct Virus <: AbstractAgent end
struct Human <: AbstractAgent end
struct Generic <: AbstractAgent end


"""
    EpiLandscape

Disease class abundances housed in the landscape. These are represented in both 2 dimensions (for computational efficiency in simulations) and 3 dimensions (to represent disease classes, their abundances and position in the grid).

"""
mutable struct EpiLandscape{T<:AbstractAgent}
  matrix::Matrix{Int64}
  grid::Array{Int64, 3}
  seed::Vector{MersenneTwister}

  EpiLandscape(A::Matrix{Int64}, args...) = EpiLandscape(Generic, A, args...)
  function EpiLandscape(Generic, abun::Matrix{Int64}, dimension::Tuple)
    a = abun
    return new{Generic}(a, reshape(a, dimension), [MersenneTwister(rand(UInt)) for _ in 1:Threads.nthreads()])
  end
  function EpiLandscape(Generic, abun::Matrix{Int64}, dimension::Tuple, seed::Vector{MersenneTwister})
    a = abun
    return new{Generic}(a, reshape(a, dimension), seed)
  end
end

import Base.copy
function copy(gl::EpiLandscape)
    return EpiLandscape(copy(gl.matrix), size(gl.grid))
end

function Base.isapprox(gl_1::EpiLandscape{T}, gl_2::EpiLandscape{T}; kwargs...) where {T<:AbstractAgent}
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
  return EpiLandscape(mat, dimension)
end
