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
  dist::Unitful.Length{Float64}
  thresh::Float64

  function GaussianKernel(dispersaldist::Unitful.Length{Float64}, pthresh::Float64)
    dist = uconvert(km, dispersaldist)
    new(dist, pthresh)
  end
end


mutable struct LongTailKernel <: AbstractKernel
  dist::Unitful.Length{Float64}
  shape::Float64
  thresh::Float64
  
  function LongTailKernel(dispersaldist::Unitful.Length{Float64}, shape::Float64, pthresh::Float64)
    dist = uconvert(km, dispersaldist)
    new(dist, shape, pthresh)
  end
end

GLOBAL_typedict["GaussianKernel"] = GaussianKernel

abstract type BoundaryCondition end

mutable struct Cylinder <: BoundaryCondition end
mutable struct Torus <: BoundaryCondition end
mutable struct NoBoundary <: BoundaryCondition end

GLOBAL_typedict["Cylinder"] = Cylinder
GLOBAL_typedict["Torus"] = Torus
GLOBAL_typedict["NoBoundary"] = NoBoundary

mutable struct BirthOnlyMovement{K <: AbstractKernel, B <: BoundaryCondition} <: AbstractMovement
  kernels::Vector{K}
  boundary::B
end
function BirthOnlyMovement(kernels::Vector{K}) where K <: AbstractKernel
    return BirthOnlyMovement{K, NoBoundary}(kernel, NoBoundary())
end
GLOBAL_typedict["BirthOnlyMovement"] = BirthOnlyMovement
mutable struct AlwaysMovement{K <: AbstractKernel, B <: BoundaryCondition} <: AbstractMovement
  kernels::Vector{K}
  boundary::B
end
function AlwaysMovement(kernels::Vector{K}) where K <: AbstractKernel
    return AlwaysMovement{K, NoBoundary}(kernel, NoBoundary())
end
GLOBAL_typedict["AlwaysMovement"] = AlwaysMovement

mutable struct NoMovement{K <: AbstractKernel} <: AbstractMovement
  kernels::Vector{K}
end
GLOBAL_typedict["NoMovement"] = NoMovement
getkernels(m::BirthOnlyMovement) = m.kernels
getkernels(m::AlwaysMovement) = m.kernels
getkernels(m::NoMovement) = m.kernels

getboundary(m::BirthOnlyMovement) = m.boundary
getboundary(m::AlwaysMovement) = m.boundary
