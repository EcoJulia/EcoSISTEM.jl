"""
    AbstractMovement

Abstract supertype of movements
"""
abstract type AbstractMovement end

abstract type AbstractKernel end

"""
    GaussianKernel <: AbstractKernel

GaussianMovement holds parameters for a gaussian movement kernel; a
dispersal variance for a species, `var`, and a threshold, `thresh`, beyond which
dispersal cannot take place.
"""
mutable struct GaussianKernel <: AbstractKernel
    dist::Unitful.Length{Float64}
    thresh::Float64

    function GaussianKernel(dispersaldist::Unitful.Length{Float64},
                            pthresh::Float64)
        dist = uconvert(km, dispersaldist)
        return new(dist, pthresh)
    end
end

"""
    LongTailKernel <: AbstractKernel

LongTailKernel holds parameters for a movement kernel; a
dispersal variance for a species, `var`, and a threshold, `thresh`, beyond which dispersal cannot take place.
"""
mutable struct LongTailKernel <: AbstractKernel
    dist::Unitful.Length{Float64}
    shape::Float64
    thresh::Float64

    function LongTailKernel(dispersaldist::Unitful.Length{Float64},
                            shape::Float64, pthresh::Float64)
        dist = uconvert(km, dispersaldist)
        return new(dist, shape, pthresh)
    end
end

"""
    BoundaryCondition

An abstract type for what should happen at the boundaries of an ecosystem.
"""
abstract type BoundaryCondition end

"""
    Cylinder <: BoundaryCondition

A cylindrical boundary where species can cross the x boundary but not the y.
"""
mutable struct Cylinder <: BoundaryCondition end
"""
    Torus <: BoundaryCondition

A toroidal boundary where species can cross both boundaries.
"""
mutable struct Torus <: BoundaryCondition end
"""
    NoBoundary <: BoundaryCondition

A hard boundary where no species can cross.
"""
mutable struct NoBoundary <: BoundaryCondition end

"""
    BirthOnlyMovement{K <: AbstractKernel, B <: BoundaryCondition} <: AbstractMovement

Movement can only happen to individuals that have just been born ("plant-like").
"""
mutable struct BirthOnlyMovement{K <: AbstractKernel, B <: BoundaryCondition} <:
               AbstractMovement
    kernels::Vector{K}
    boundary::B
end
function BirthOnlyMovement(kernels::Vector{K}) where {K <: AbstractKernel}
    return BirthOnlyMovement{K, NoBoundary}(kernels, NoBoundary())
end

"""
    AlwaysMovement{K <: AbstractKernel, B <: BoundaryCondition} <: AbstractMovement

Movement can happen to any individual ("animal-like").
"""
mutable struct AlwaysMovement{K <: AbstractKernel, B <: BoundaryCondition} <:
               AbstractMovement
    kernels::Vector{K}
    boundary::B
end
function AlwaysMovement(kernels::Vector{K}) where {K <: AbstractKernel}
    return AlwaysMovement{K, NoBoundary}(kernels, NoBoundary())
end

"""
    NoMovement{K <: AbstractKernel, B <: BoundaryCondition} <: AbstractMovement

No movement can take place.
"""
mutable struct NoMovement{K <: AbstractKernel} <: AbstractMovement
    kernels::Vector{K}
end

getkernels(m::BirthOnlyMovement) = m.kernels
getkernels(m::AlwaysMovement) = m.kernels
getkernels(m::NoMovement) = m.kernels

getboundary(m::BirthOnlyMovement) = m.boundary
getboundary(m::AlwaysMovement) = m.boundary
