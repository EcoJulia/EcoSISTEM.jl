# SPDX-License-Identifier: LGPL-3.0-or-later
#
# *Where* an intervention acts. Operations after the first share one resolved region, which is the
# only way to act twice on the same randomly chosen cells.

using Unitful
using Unitful: 𝐓

# An integer number of cells, or a rate: a real number per unit time. The dimension is pinned, so a
# length or a bare float is refused by the signature rather than by a check further in.
const CELLCOUNT = Union{Integer, Unitful.Quantity{<:Real, 𝐓^-1}}

"""
    AbstractRegion

**Where** an [`Intervention`](@ref) applies - [`AllCells`](@ref), [`ActiveCells`](@ref),
[`CellMask`](@ref), [`RandomCells`](@ref) or [`SpreadingCells`](@ref).

A region resolves to the linear cell indices an operation should touch. The two random regions draw
from a **counter-based** stream keyed on the step (see [`Intervention`](@ref)), never from a species'
stream and never from the global RNG, so a selection is identical on every rank and thread and a run
replays exactly.
"""
abstract type AbstractRegion end

"""    AllCells() <: AbstractRegion - every cell in the grid, active or not. """
struct AllCells <: AbstractRegion end

"""    ActiveCells() <: AbstractRegion - only the cells currently marked active. """
struct ActiveCells <: AbstractRegion end

"""
    CellMask(mask::AbstractMatrix{Bool})

The cells where `mask` is true.

# Arguments

  - `mask`: a `(Y, X)` boolean grid, the same shape as the habitat's own.
"""
struct CellMask{M <: AbstractMatrix{Bool}} <: AbstractRegion
    mask::M
end

"""
    RandomCells(count::CELLCOUNT)

`count` cells drawn at random, without replacement, from the currently active ones.

# Arguments

  - `count`: how many cells, either exactly (`RandomCells(20)`) or as a per-time rate
    (`RandomCells(0.05/year)`). A rate is applied to each candidate cell independently over the
    step, so the number taken is a binomial draw and varies from one firing to the next.
"""
struct RandomCells{C <: CELLCOUNT} <: AbstractRegion
    count::C

    RandomCells(count::CELLCOUNT) = new{typeof(count)}(count)
end

"""
    SpreadingCells(count::CELLCOUNT)

`count` cells forming a contiguous cluster: a random active seed cell, then repeatedly a random
active neighbour of the cells already chosen. Falls back to a fresh seed if the cluster is boxed in,
so it always returns `count` cells where that many are available.

# Arguments

  - `count`: how many cells, exactly or as a per-time rate, as [`RandomCells`](@ref).
"""
struct SpreadingCells{C <: CELLCOUNT} <: AbstractRegion
    count::C

    SpreadingCells(count::CELLCOUNT) = new{typeof(count)}(count)
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# As for the schedules: these nest inside an `Intervention`, and the default prints the whole type
# signature where the caller wrote `RandomCells(20)`. A `CellMask` reports its shape rather than its
# contents - the mask is data, and printing it is what these methods exist to prevent.
Base.show(io::IO, r::RandomCells) = print(io, "RandomCells($(r.count))")
function Base.show(io::IO, r::SpreadingCells)
    return print(io, "SpreadingCells($(r.count))")
end
function Base.show(io::IO, r::CellMask)
    ny, nx = size(r.mask)
    return print(io, "CellMask($(ny) × $(nx), $(count(r.mask)) selected)")
end
