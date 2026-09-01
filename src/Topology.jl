# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What happens at the edge of the grid: a boundary condition per dimension, and the pairing of two
# of them that a habitat carries.

"""
    AbstractBoundaryCondition

What happens at **one** edge pair of a grid: [`Periodic`](@ref) (the edges join) or
[`Bounded`](@ref) (they do not). A condition on a single axis, not on the grid as a whole.
"""
abstract type AbstractBoundaryCondition end

"""    Bounded <: AbstractBoundaryCondition - this axis does not wrap: it has real edges. """
struct Bounded <: AbstractBoundaryCondition end

"""    Periodic <: AbstractBoundaryCondition - this axis wraps: its two edges join. """
struct Periodic <: AbstractBoundaryCondition end

"""
    AbstractTopology

The shape of the simulated world - how its edges join, if they join at all.

A property of the **grid**, held and set on [`GridHabitat`](@ref), so every species on one grid
shares it. The name promises no particular geometry: [`EdgeTopology`](@ref) is a rectangular lattice
and is the only subtype at present, but a world that is not one can subtype this.

What becomes of an individual dispersing into a *dead* cell is a separate question, answered per
species by `disperse_safely` on the movement rather than by the map.
"""
abstract type AbstractTopology end

"""
    EdgeTopology{BCY, BCX} <: AbstractTopology
    EdgeTopology(; y, x)

The topology of a rectangular grid: one [`AbstractBoundaryCondition`](@ref) per axis.

Three of the four combinations have names of their own - [`Torus`](@ref), [`Cylinder`](@ref) and
[`Island`](@ref) - and those are the usual way to write them. The keyword form is for the fourth,
which has no ecological name.

# Arguments

  - `y`: the boundary condition for the **rows** - northing on a projected grid - given as a type,
    `EdgeTopology(y = Periodic, x = Bounded)`.
  - `x`: the same for the **columns**, easting.

# Type parameters

  - `BCY`, `BCX`: those two conditions, in that order. `(y, x)` throughout the package, matching
    `DimensionalData`'s own `Y`/`X` dims.
"""
struct EdgeTopology{BCY <: AbstractBoundaryCondition,
                    BCX <: AbstractBoundaryCondition} <: AbstractTopology end

"""
    Torus = EdgeTopology{Periodic, Periodic}

A toroidal grid: both pairs of edges join, so nothing ever leaves.

True of a synthetic grid. On one with a real-world position it asserts something false - a study
area is a window on a globe, not a world that wraps - so it is accepted with a warning.
"""
const Torus = EdgeTopology{Periodic, Periodic}

"""
    Cylinder = EdgeTopology{Bounded, Periodic}

A cylindrical grid: the **east-west** (X) edges join, the north-south ones do not. The parameter
order is `{Y, X}`, so `Bounded` comes first - X is the axis that wraps.

On a grid with a real-world position a wrapping X is right only where the grid spans the whole
longitude sweep, and is warned about otherwise. Where it does span it, `Cylinder` is exact rather
than an approximation, and it is the only wrapping topology that can be: a wrapping Y is never
right, because latitude does not wrap however wide the grid.
"""
const Cylinder = EdgeTopology{Bounded, Periodic}

"""
    Island = EdgeTopology{Bounded, Bounded}

A grid with hard edges: neither pair joins. **The default on every grid.**

A step past an edge is treated exactly as a step into a dead cell, so what becomes of it is answered
per species by `disperse_safely`. With its default of `true` the blocked weight is redistributed
among the reachable destinations, making the edge **reflecting** rather than absorbing.
"""
const Island = EdgeTopology{Bounded, Bounded}
