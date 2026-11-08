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
[`CellMask`](@ref), [`CellsAt`](@ref), [`RandomCells`](@ref) or [`SpreadingCells`](@ref). An
`Intervention` also takes a shape as its region, as a study area's `within` does: a
[`NaturalEarthSpec`](@ref) or any other [`AbstractShapeSpec`](@ref) means the cells it covers at
least half of, and a [`ShapeMaskSpec`](@ref) says which cells under another rule.

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

"""
    CellsAt(places...; snap = 0.0km)

The cells at `places`: a site named by where it is rather than by its position in the grid.

```julia
Intervention(AtTime(10year), CellsAt(LatLong(55.95°, -3.19°)), AddAbundance("Rhododendron", 500))
```

Each place is found in the grid by [`getcellat`](@ref), and a place off the grid is refused before
the run's first step. The cells are taken whether they are active or not, as a [`CellMask`](@ref)'s
are, so a [`Reactivate`](@ref) at a site reaches it. With `snap`, a place in an inactive cell is
moved to the nearest active cell whose centre lies within `snap` of its own, ties to the lower
`(y, x)`, and stays where it is when none is that near. Two places in one cell are that cell once.

# Arguments

  - `places`: one or more places, each a [`LatLong`](@ref), a [`SpatialLocation`](@ref) in the grid's
    own frame, or a grid key, `CartesianIndex(y, x)`.
  - `snap`: how far a place in an inactive cell may move to reach an active one, as a length; none
    by default.
"""
struct CellsAt{P <: AbstractVector, S <: Unitful.Length} <: AbstractRegion
    places::P
    snap::S

    function CellsAt(places::Union{SpatialLocation, CartesianIndex{2}}...;
                     snap::Unitful.Length = 0.0km)
        isempty(places) &&
            error("`CellsAt` needs at least one place: a `LatLong`, a `SpatialLocation` or a " *
                  "grid key, `CartesianIndex(y, x)`.")
        _checksnap(snap)
        given = collect(places)
        return new{typeof(given), typeof(snap)}(given, snap)
    end
end

# A shape as a region, as `within` takes one: `given` as written, the mask spec it means, and the
# cells it covers in the one ecosystem it was last resolved for, with the records of the shape's
# files. Resolving rasterises the shape on the habitat's grid, a collective build under MPI, so it is
# done once for each ecosystem rather than at every firing.
struct _ShapeCells{G, M <: ShapeMaskSpec} <: AbstractRegion
    given::G
    spec::M
    found::Base.RefValue{Any}
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
function Base.show(io::IO, r::CellsAt)
    print(io, "CellsAt(", join(map(p -> sprint(show, p), r.places), ", "))
    iszero(r.snap) || print(io, ", snap = ", r.snap)
    return print(io, ")")
end
Base.show(io::IO, r::_ShapeCells) = show(io, r.given)
