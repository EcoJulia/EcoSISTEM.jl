# SPDX-License-Identifier: LGPL-3.0-or-later
#
# How a change value relates to the layer it changes: a position, an interval, or a rate.

using Unitful

"""
    AbstractChangeMode

How a change's values are to be read, and so which unit they must be in. The three that carry
values differ in exactly that: [`AbsoluteChange`](@ref) is a *position* in the layer's own unit,
[`RelativeChange`](@ref) an *interval* from the layer's captured values, and [`RateChange`](@ref)
an interval per unit time. [`NoChange`](@ref) carries no values at all.

A mode is a type parameter of [`AbstractLayerChange`](@ref), so every change declares one and
`changeunit` always has an answer.
"""
abstract type AbstractChangeMode end

"""    AbsoluteChange <: AbstractChangeMode — a position: the layer's value *is* this, in the layer's unit. """
struct AbsoluteChange <: AbstractChangeMode end

"""    RelativeChange <: AbstractChangeMode — an interval from the layer's captured values, so `K` rather than `°C`. """
struct RelativeChange <: AbstractChangeMode end

"""    RateChange <: AbstractChangeMode — an interval per unit time, accumulated each step. """
struct RateChange <: AbstractChangeMode end

"""    NoChange <: AbstractChangeMode — no values, and so no `changeunit`. """
struct NoChange <: AbstractChangeMode end
