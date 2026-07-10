# SPDX-License-Identifier: LGPL-3.0-or-later

"""
    AbstractMovement

Abstract supertype of movements
"""
abstract type AbstractMovement end

abstract type AbstractKernel end

"""
    GaussianKernel <: AbstractKernel

GaussianMovement holds parameters for a gaussian movement kernel; a dispersal
variance for a species, `var`, and a threshold, `thresh`, beyond which dispersal
cannot take place.
"""
struct GaussianKernel <: AbstractKernel
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

LongTailKernel holds parameters for a movement kernel; a dispersal variance for
a species, `var`, and a threshold, `thresh`, beyond which dispersal cannot take
place.
"""
struct LongTailKernel <: AbstractKernel
    dist::Unitful.Length{Float64}
    shape::Float64
    thresh::Float64

    function LongTailKernel(dispersaldist::Unitful.Length{Float64},
                            shape::Float64,
                            pthresh::Float64)
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
struct Cylinder <: BoundaryCondition end
"""
    Torus <: BoundaryCondition

A toroidal boundary where species can cross both boundaries.
"""
struct Torus <: BoundaryCondition end
"""
    NoBoundary <: BoundaryCondition

A hard boundary where no species can cross.
"""
struct NoBoundary <: BoundaryCondition end

"""
    BirthOnlyMovement{K <: AbstractKernel, B <: BoundaryCondition} <: AbstractMovement

Movement can only happen to individuals that have just been born ("plant-like").
"""
struct BirthOnlyMovement{K <: AbstractKernel, B <: BoundaryCondition} <:
       AbstractMovement
    kernels::Vector{K}
    boundary::B
end
function BirthOnlyMovement(kernels::Vector{K}) where {K <: AbstractKernel}
    return BirthOnlyMovement{K, NoBoundary}(kernels, NoBoundary())
end
@doc (@doc BirthOnlyMovement) BirthOnlyMovement(::Vector{K} where {K <:
                                                                   AbstractKernel})

"""
    AlwaysMovement{K <: AbstractKernel, B <: BoundaryCondition} <: AbstractMovement

Movement can happen to any individual ("animal-like").
"""
struct AlwaysMovement{K <: AbstractKernel, B <: BoundaryCondition} <:
       AbstractMovement
    kernels::Vector{K}
    boundary::B
end
function AlwaysMovement(kernels::Vector{K}) where {K <: AbstractKernel}
    return AlwaysMovement{K, NoBoundary}(kernels, NoBoundary())
end
@doc (@doc AlwaysMovement) AlwaysMovement(::Vector{K} where {K <:
                                                             AbstractKernel})

"""
    NoMovement{K <: AbstractKernel, B <: BoundaryCondition} <: AbstractMovement

No movement can take place.
"""
struct NoMovement{K <: AbstractKernel} <: AbstractMovement
    kernels::Vector{K}
end

"""
    getkernels(m::AbstractMovement)

Extract the vector of dispersal kernels from a movement type.
"""
getkernels(m::BirthOnlyMovement) = m.kernels
getkernels(m::AlwaysMovement) = m.kernels
getkernels(m::NoMovement) = m.kernels

"""
    getboundary(m::AbstractMovement)

Extract the boundary condition from a movement type.
"""
getboundary(m::BirthOnlyMovement) = m.boundary
getboundary(m::AlwaysMovement) = m.boundary
