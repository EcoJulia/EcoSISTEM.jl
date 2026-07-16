# SPDX-License-Identifier: LGPL-3.0-or-later

# ---------------------------------------------------------------------------
# Niche axes
# ---------------------------------------------------------------------------
# A `NicheAxis` names *what* a layer measures — an axis of the (Hutchinsonian) niche — in
# plain English, independent of its physical unit, and orthogonal to `Role` (which says how
# a layer is *used*). A layer and the trait matched to it share the same `NicheAxis`, so
# "same axis" is the matching rule rather than "same unit". Users add their own with
# `struct MyAxis <: NicheAxis end` plus, if needed, the interface overrides below.

"""
    NicheAxis

Abstract supertype of the niche-axis markers naming what a layer measures (e.g.
[`Temperature`](@ref), [`Rainfall`](@ref)), independent of its unit and orthogonal to
[`Role`](@ref). Extend it with your own `struct MyAxis <: NicheAxis end`.
"""
abstract type NicheAxis end

"""    Temperature <: NicheAxis """
struct Temperature <: NicheAxis end

"""    Rainfall <: NicheAxis """
struct Rainfall <: NicheAxis end

"""    SolarRadiation <: NicheAxis """
struct SolarRadiation <: NicheAxis end

"""    VolumetricWater <: NicheAxis """
struct VolumetricWater <: NicheAxis end

# --- The NicheAxis interface -----------------------------------------------
# The small set of hooks that replace unit-based inference of a layer's meaning. Each has a
# safe default; a built-in or user axis overrides only what applies. The budget/requirement/
# change types they name live in later-included files and are resolved at call time (nothing
# here runs at load). Not wired into the build path until stage B.

"""
    canonicalunit(::NicheAxis)

The unit a layer on this axis is normalised to when materialised (e.g. `Temperature → K`),
or `nothing` to leave values as-is.
"""
canonicalunit(::NicheAxis) = nothing
canonicalunit(::Temperature) = K
canonicalunit(::Rainfall) = mm

"""
    budgettype(::NicheAxis)

The `AbstractBudget` concrete type for this axis when used as a `Budget` resource, or
`nothing` if the axis is not a consumable resource (so it errors clearly rather than
silently guessing one if a budget is nonetheless requested).
"""
budgettype(::NicheAxis) = nothing
budgettype(::Rainfall) = WaterBudget
budgettype(::SolarRadiation) = SolarBudget
budgettype(::VolumetricWater) = VolWaterBudget

"""
    requirementtype(::NicheAxis)

The `AbstractRequirement` concrete type a species uses to consume this axis' resource, or
`nothing` if the axis is not a resource.
"""
requirementtype(::NicheAxis) = nothing
requirementtype(::Rainfall) = WaterRequirement
requirementtype(::SolarRadiation) = SolarRequirement
requirementtype(::VolumetricWater) = VolWaterRequirement

"""
    dynamics(::NicheAxis)

The default per-timestep change function for a layer on this axis (default
[`NoChange`](@ref)); overridden per layer at build time.
"""
dynamics(::NicheAxis) = NoChange
