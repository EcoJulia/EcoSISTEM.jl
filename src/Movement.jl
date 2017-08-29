"""
    AbstractMovement

Abstract supertype of movements
"""
abstract type AbstractMovement end


abstract type AbstractKernel end

"""
    GaussianKernel <: AbstractKernel

GaussianMovement holds parameters for a gaussian movement kernel; a vector of
dispersal variances per species, `var`, and a threshold, `thresh`, beyond which
dispersal cannot take place.
"""
mutable struct GaussianKernel <: AbstractKernel
  dist::Vector{Unitful.Length{Float64}}
  thresh::Float64
end

function GaussianKernel(dispersaldist::Unitful.Length{Float64}, numspecies::Int64,
  pthresh::Float64)
  dist = map(u-> uconvert(km, u), dispersaldist)
  GaussianKernel(repmat([dist], numspecies), pthresh)
end

mutable struct BirthOnlyMovement{K <: AbstractKernel} <: AbstractMovement
  kernel::K
end

mutable struct AlwaysMovement{K <: AbstractKernel} <: AbstractMovement
  kernel::K
end

mutable struct NoMovement{K <: AbstractKernel} <: AbstractMovement
  kernel::K
end

getkernel(m::BirthOnlyMovement) = m.kernel
getkernel(m::AlwaysMovement) = m.kernel
getkernel(m::NoMovement) = m.kernel
