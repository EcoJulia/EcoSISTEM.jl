"""
    EpiMovement{MO <: AbstractMovement} <: AbstractMovement

Movement can happen at several different levels, local gaussian processes, `localmoves`, and longer distance moves, `regionmoves`.
"""
mutable struct EpiMovement{MO1 <: AbstractMovement, MO2 <: AbstractMovement} <: AbstractMovement
    localmoves::MO1
    regionmoves::MO2
end

function EpiMovement(localkernels::Vector{K}, move_record::DataFrame) where K <: AbstractKernel
    localmoves = AlwaysMovement{K, NoBoundary}(localkernels, NoBoundary())
    regionmoves = LongDistance(move_record)
    return EpiMovement{typeof(localmoves), typeof(regionmoves)}(localmoves, regionmoves)
end

function EpiMovement(localkernels::Vector{K}) where K <: AbstractKernel
    localmoves = AlwaysMovement{K, NoBoundary}(localkernels, NoBoundary())
    regionmoves = LongDistance(DataFrame([(from=1.0, to=1.0, count=0.0)]))
    return EpiMovement{typeof(localmoves), typeof(regionmoves)}(localmoves, regionmoves)
end

mutable struct LongDistance <: AbstractMovement
    move_record::DataFrame
end


getdispersaldist(m::EpiMovement, sp::Int64) = m.localmoves.kernels[sp].dist
getdispersalvar(m::EpiMovement, sp::Int64) = (m.localmoves.kernels[sp].dist)^2 * pi / 4
