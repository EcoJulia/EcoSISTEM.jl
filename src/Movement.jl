# SPDX-License-Identifier: LGPL-3.0-or-later
#
# How offspring disperse: the dispersal kernels, the movement rules that draw from them, and the
# per-species lookup table of dispersal probabilities that the hot loop actually reads.

using Unitful

"""
    AbstractKernel

The shape of the probability distribution over where an individual's offspring lands, relative to its
own cell. [`GaussianKernel`](@ref) and [`LongTailKernel`](@ref) are the concrete kernels.

A kernel says only *how far and how likely*. What happens when a draw crosses an edge is the grid's
[`AbstractTopology`](@ref), and what becomes of one aimed at a dead cell is `disperse_safely` on the
movement - three separate questions, answered by three separate things.
"""
abstract type AbstractKernel end

"""
    GaussianKernel(dist::Unitful.Length, thresh::Float64)

A dispersal kernel whose distances are Rayleigh distributed - the two-dimensional case of dispersing
a Gaussian displacement in each direction independently.

# Fields

  - `dist`: the species' dispersal distance, converted to kilometres on construction.
  - `thresh`: the probability below which a destination is dropped, which is what bounds the lookup
    table to a finite neighbourhood.
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
    LongTailKernel(dist::Unitful.Length, shape::Float64, thresh::Float64)

A dispersal kernel with a heavier tail than [`GaussianKernel`](@ref)'s, so that rare long-distance
dispersal is not effectively impossible.

# Fields

  - `dist`: the species' dispersal distance, converted to kilometres on construction.
  - `shape`: how heavy the tail is. Smaller values reach further.
  - `thresh`: the probability below which a destination is dropped, which is what bounds the lookup
    table to a finite neighbourhood.
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
    AbstractMovement

Which individuals disperse at all: [`BirthOnlyMovement`](@ref) (only the newly born, plant-like),
[`AlwaysMovement`](@ref) (any individual, animal-like) or [`NoMovement`](@ref) (none).

A movement carries the per-species part of dispersal - one [`AbstractKernel`](@ref) per species, and
one `disperse_safely` flag per species. The grid's own edges are not part of it: they belong to the
[`GridHabitat`](@ref), because two species on one grid cannot be on different topologies.

`disperse_safely` says what becomes of an individual dispersing into a **dead cell** - one off the
grid, or an inactive one. `true` redistributes it among the reachable destinations as though it had
never aimed there; `false` loses it, so the share of the kernel pointing at dead cells becomes
mortality.

**It is per species because it is a fact about the disperser**, not about the map. A wind-dispersed
seed blown out to sea is gone; an animal-dispersed one is carried to somewhere the animal can
actually reach, so two species on one grid may legitimately differ. It cannot empty a cell outright:
staying put is one of the kernel's own destinations, and a cell only disperses while it is active, so
some share always survives.
"""
abstract type AbstractMovement end

"""
    BirthOnlyMovement{K <: AbstractKernel} <: AbstractMovement

Only individuals that have just been born disperse - plant-like.

# Fields

  - `kernels`: one [`AbstractKernel`](@ref) per species.
  - `disperse_safely`: one flag per species, as described on [`AbstractMovement`](@ref). Defaults to
    `true` for every species.
"""
struct BirthOnlyMovement{K <: AbstractKernel} <: AbstractMovement
    kernels::Vector{K}
    disperse_safely::Vector{Bool}

    function BirthOnlyMovement(kernels::Vector{K},
                               disperse_safely::AbstractVector{Bool} = fill(true,
                                                                            length(kernels))) where {K <:
                                                                                                     AbstractKernel}
        _checkdispersesafely(kernels, disperse_safely)
        return new{K}(kernels, collect(disperse_safely))
    end
end

"""
    AlwaysMovement{K <: AbstractKernel} <: AbstractMovement

Any individual may disperse, not only the newly born - animal-like.

# Fields

  - `kernels`: one [`AbstractKernel`](@ref) per species.
  - `disperse_safely`: one flag per species, as described on [`AbstractMovement`](@ref). Defaults to
    `true` for every species.
"""
struct AlwaysMovement{K <: AbstractKernel} <: AbstractMovement
    kernels::Vector{K}
    disperse_safely::Vector{Bool}

    function AlwaysMovement(kernels::Vector{K},
                            disperse_safely::AbstractVector{Bool} = fill(true,
                                                                         length(kernels))) where {K <:
                                                                                                  AbstractKernel}
        _checkdispersesafely(kernels, disperse_safely)
        return new{K}(kernels, collect(disperse_safely))
    end
end

"""
    NoMovement{K <: AbstractKernel} <: AbstractMovement

Nothing disperses: every individual stays in the cell it was born in.

# Fields

  - `kernels`: one [`AbstractKernel`](@ref) per species, carried but never drawn from, so that a
    model can be switched between movements without rebuilding its species.
"""
struct NoMovement{K <: AbstractKernel} <: AbstractMovement
    kernels::Vector{K}
end

"""
    Lookup

One species' dispersal neighbourhood, precomputed once and read by the hot loop: every destination
its kernel can reach, as an offset from the source cell, with the probability of landing there.

The last two fields are scratch, rewritten on every cell of every timestep. Holding them here rather
than allocating them per call is what keeps the movement step allocation-free.

# Fields

  - `y`, `x`: the destination offsets, relative to the cell dispersing. In that order, as everywhere
    else in the package.
  - `p`: the kernel's probability for each, fixed for the run.
  - `pnew`: those probabilities renormalised over the destinations actually available from this cell.
  - `moves`: how many individuals go to each destination this step.
"""
struct Lookup
    y::Vector{Int64}
    x::Vector{Int64}
    p::Vector{Float64}
    pnew::Vector{Float64}
    moves::Vector{Int64}
end
# == Functions ==================================================================================

# One line, because the default prints all five vectors - measured at 230 320 characters for a
# single species' neighbourhood. What identifies a lookup is how far it reaches, not the numbers in
# it: the last three fields are scratch rewritten on every cell of every timestep.
function Base.show(io::IO, l::Lookup)
    return print(io, "Lookup($(length(l.y)) destinations)")
end

# One flag per kernel, or the two vectors silently describe different species - `zip` would truncate
# and the last species would take a neighbour's setting.
function _checkdispersesafely(kernels, disperse_safely)
    length(disperse_safely) == length(kernels) ||
        error("`disperse_safely` has $(length(disperse_safely)) entries but there are " *
              "$(length(kernels)) dispersal kernels: it is a per-species setting, so it needs one " *
              "entry per species (or a single value applied to all).")
    return nothing
end

# Accepts (and ignores) a `disperse_safely` vector, so `NoMovement` can be constructed the same way
# as `BirthOnlyMovement`/`AlwaysMovement` when the movement type is chosen as a value and called
# uniformly (`movement(kernels, disperse_safely)`). Nothing disperses, so nothing can be lost.
function NoMovement(kernels::Vector{K},
                    ::AbstractVector{Bool}) where {K <: AbstractKernel}
    return NoMovement(kernels)
end
@doc (@doc NoMovement) NoMovement(::Vector{K} where {K <: AbstractKernel},
                                  ::AbstractVector{Bool})

"""
    getkernels(m::AbstractMovement)

Return the vector of dispersal kernels a movement holds, one per species.
"""
getkernels(m::BirthOnlyMovement) = m.kernels
getkernels(m::AlwaysMovement) = m.kernels
getkernels(m::NoMovement) = m.kernels

"""
    dispersesafely(m::AbstractMovement, sp::Integer)

Return whether species `sp` is redistributed (`true`) or lost (`false`) when its dispersal is aimed
at a dead cell. `NoMovement` answers `true`, since nothing disperses and so nothing can be lost.
"""
dispersesafely(m::BirthOnlyMovement, sp::Integer) = m.disperse_safely[sp]
dispersesafely(m::AlwaysMovement, sp::Integer) = m.disperse_safely[sp]
dispersesafely(::NoMovement, ::Integer) = true
