# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Where a decided grid property came from — the provenance a `StudyAreaReport` records beside each
# number it reports, so that "0.5°" and "0.5° because you asked for it" are different statements.

"""
    AbstractDecisionSource

Where a study area's `crs` or `cellsize` came from — [`GivenByUser`](@ref) if you supplied it,
otherwise how it was derived. [`StudyAreaReport`](@ref) records one beside each of those two values,
which is what lets the area announce what it **guessed** without also announcing what you gave it.
"""
abstract type AbstractDecisionSource end

"""    GivenByUser <: AbstractDecisionSource — you supplied it. """
struct GivenByUser <: AbstractDecisionSource end

"""    AdoptedFromLayers <: AbstractDecisionSource — a CRS adopted from the layers. """
struct AdoptedFromLayers <: AbstractDecisionSource end

"""    NoRealWorldPosition <: AbstractDecisionSource — a synthetic area, with no real-world position for a CRS to describe. """
struct NoRealWorldPosition <: AbstractDecisionSource end

"""    TakenFromAlignedLayer <: AbstractDecisionSource — a cell size taken from the layer being aligned to. """
struct TakenFromAlignedLayer <: AbstractDecisionSource end

"""    AgreedByAllLayers <: AbstractDecisionSource — a cell size every layer agrees on. """
struct AgreedByAllLayers <: AbstractDecisionSource end

"""    MeasuredAcrossProjection <: AbstractDecisionSource — a cell size measured across the projection, no layer being in the target CRS. """
struct MeasuredAcrossProjection <: AbstractDecisionSource end

"""    RoundedFromMeasurement <: AbstractDecisionSource — a [`MeasuredAcrossProjection`](@ref) cell size, floored to whole units. """
struct RoundedFromMeasurement <: AbstractDecisionSource end

# --- Saying it in prose ------------------------------------------------------

# How a decision is described, given the source itself: the phrase completes the parenthesis in
# `StudyAreaReport`'s `show` and in the `@info` lines an investigated area emits.
#
# Dispatch is on the source's own type, and must stay that way. Reaching these through a `Symbol`
# promoted into a `Val` would be dynamic dispatch with none of `Val`'s compile-time payoff, because
# the compiler cannot see a tag that only arrives at runtime.
_sourcephrase(::GivenByUser) = "given"

_sourcephrase(::AdoptedFromLayers) = "adopted from the layers"

_sourcephrase(::TakenFromAlignedLayer) = "the aligned layer's own"

_sourcephrase(::AgreedByAllLayers) = "the layers agree on it"

function _sourcephrase(::MeasuredAcrossProjection)
    return "measured across the projection, since no layer is in the target CRS — approximate, pass `cellsize` for an exact one"
end

function _sourcephrase(::RoundedFromMeasurement)
    return "measured across the projection and rounded down to whole units — pass `cellsize` for an exact one"
end

_sourcephrase(::NoRealWorldPosition) = "no real-world position"
