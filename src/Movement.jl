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

function GaussianKernel(dispersaldist::Unitful.Length{Float64}, numspecies::Int64, pthresh::Float64)
  dist = map(u-> uconvert(km, u), dispersaldist)
  GaussianKernel(fill(dist, numspecies), pthresh)
end

mutable struct LongTailKernel <: AbstractKernel
  dist::Vector{Unitful.Length{Float64}}
  shape::Vector{Float64}
  thresh::Float64
end
function LongTailKernel(dispersaldist::Unitful.Length{Float64}, shape::Float64, numspecies::Int64, pthresh::Float64)
  dist = map(u-> uconvert(km, u), dispersaldist)
  LongTailKernel(fill(dist, numspecies), fill(shape, numspecies), pthresh)
end


abstract type BoundaryCondition end

mutable struct Cylinder <: BoundaryCondition end
mutable struct Torus <: BoundaryCondition end
mutable struct NoBoundary <: BoundaryCondition end

mutable struct BirthOnlyMovement{K <: AbstractKernel, B <: BoundaryCondition} <: AbstractMovement
  kernel::K
  boundary::B
end
function BirthOnlyMovement(kernel::K) where K <: AbstractKernel
    return BirthOnlyMovement{K, NoBoundary}(kernel, NoBoundary())
end

mutable struct AlwaysMovement{K <: AbstractKernel, B <: BoundaryCondition} <: AbstractMovement
  kernel::K
  boundary::B
end
function AlwaysMovement(kernel::K) where K <: AbstractKernel
    return AlwaysMovement{K, NoBoundary}(kernel, NoBoundary())
end

mutable struct NoMovement{K <: AbstractKernel} <: AbstractMovement
  kernel::K
end

getkernel(m::BirthOnlyMovement) = m.kernel
getkernel(m::AlwaysMovement) = m.kernel
getkernel(m::NoMovement) = m.kernel

getboundary(m::BirthOnlyMovement) = m.boundary
getboundary(m::AlwaysMovement) = m.boundary
