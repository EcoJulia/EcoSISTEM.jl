"""
    AbstractMovement

Abstract supertype of movements
"""
abstract AbstractMovement


abstract AbstractKernel

"""
    GaussianKernel <: AbstractKernel

GaussianMovement holds parameters for a gaussian movement kernel; a vector of
dispersal variances per species, `var`, and a threshold, `thresh`, beyond which
dispersal cannot take place.
"""
type GaussianKernel <: AbstractKernel
  var::Vector{Float64}
  thresh::Float64
end

function GaussianKernel(dispersalvar::Float64, numspecies::Int64,
  pthresh::Float64)
  GaussianKernel(repmat([dispersalvar], numspecies), pthresh)
end

type BirthOnlyMovement{K <: AbstractKernel} <: AbstractMovement
  kernel::K
end

type AlwaysMovement{K <: AbstractKernel} <: AbstractMovement
  kernel::K
end

getkernel(m::BirthOnlyMovement) = m.kernel
getkernel(m::AlwaysMovement) = m.kernel
