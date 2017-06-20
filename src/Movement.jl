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
  var::Vector{Float64}
  thresh::Float64
end

function GaussianKernel(dispersalvar::Float64, numspecies::Int64,
  pthresh::Float64)
  GaussianKernel(repmat([dispersalvar], numspecies), pthresh)
end

mutable struct BirthOnlyMovement{K <: AbstractKernel} <: AbstractMovement
  kernel::K
end

mutable struct AlwaysMovement{K <: AbstractKernel} <: AbstractMovement
  kernel::K
end

getkernel(m::BirthOnlyMovement) = m.kernel
getkernel(m::AlwaysMovement) = m.kernel
