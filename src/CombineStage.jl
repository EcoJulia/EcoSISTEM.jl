# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Which grid a `ConstructedSpec`'s combine runs on — the source's own, or the target the study
# area decided.

"""
    AbstractCombineStage

When a [`ConstructedSpec`](@ref)'s `combine` runs, relative to putting its layers on the study
grid: on the target grid ([`CombineOnTargetGrid`](@ref), the default) or on the layers' own
([`CombineOnSourceGrid`](@ref)). A type rather than a symbol, so an unrecognised stage is refused
by the signature where it is written rather than by a check inside the constructor.
"""
abstract type AbstractCombineStage end

"""
    CombineOnSourceGrid <: AbstractCombineStage

Combine on the layers' own shared grid first, and sample the *result* onto the study grid.

Needed whenever the combine does not commute with regridding: one that looks beyond its own cell
(a neighbourhood, a crop, an aggregation), and equally one that is cell-wise but **nonlinear**,
since the ratio of two interpolations is not the interpolation of the ratio. Growing-season
rainfall `gsp / gsl` is the case that matters: divided early it is a mean daily rate, divided late
a total over a mean season length, and those are different quantities.

The layers must share one native grid, checked when they are read — combining across datasets of
different resolutions has no meaning before either is resampled.
"""
struct CombineOnSourceGrid <: AbstractCombineStage end

"""
    CombineOnTargetGrid <: AbstractCombineStage

Sample every layer onto the study grid first, and combine there — the default, and right whenever
the combine commutes with regridding. It is what keeps `compress_landcover` interpolating the
per-class *percentages* it is given rather than argmaxing them first and interpolating between the
resulting class codes, which is meaningless.
"""
struct CombineOnTargetGrid <: AbstractCombineStage end
