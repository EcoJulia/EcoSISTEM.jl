# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Over what period a layer's values accumulate - a fixed span, one per slice of a stack, or a span
# that varies per cell (a growing season). What turns a declared amount into a rate.

using Unitful

"""
    AbstractAccumulationPeriod

The interval a layer's value accumulated over - what turns a total into an honest rate. One of
[`ConstantAccumulationPeriod`](@ref), [`PerSliceAccumulationPeriod`](@ref) or
[`PerCellAccumulationPeriod`](@ref); `nothing` where no period applies.

A type family rather than a bare unit because the three cases are not interchangeable: a constant
period is a fixed scale factor that a layer and a species tolerance both convert by, while a period
that *varies* makes the two readings different biological claims, and the code must be able to tell
which it is holding.
"""
abstract type AbstractAccumulationPeriod end

"""
    ConstantAccumulationPeriod(duration)

One fixed interval for the whole layer - `year` for a heat sum, `month_mean_duration` for a layer
whose month is unknowable. The stock->rate conversion is then a constant scale factor.

# Fields

  - `duration`: the interval itself, as a `Unitful.Units` - what an accumulated value is divided by
    to become a rate.
"""
struct ConstantAccumulationPeriod <: AbstractAccumulationPeriod
    duration::Unitful.Units
end

"""
    PerSliceAccumulationPeriod()

Each time slice accumulated over **its own calendar month**, so slice 1 is January's 31 days and
slice 2 February's. Written `perslice=calendar_month` in the shipped tables.

This is the case a single unit cannot express, and the reason the column is not simply a unit: a
12-slice monthly read has twelve different divisors, not one.
"""
struct PerSliceAccumulationPeriod <: AbstractAccumulationPeriod end

"""
    PerCellAccumulationPeriod(code)

The period is **another layer**, varying by cell - `gsp` (growing-season precipitation) accumulates
over `gsl` (growing-season length). Written `percell=gsl`.

A trait cannot be converted without reading that layer, so a stock and a rate reading of such a
value are genuinely different hypotheses rather than two spellings of one.

# Fields

  - `code`: the [`CODE_TYPE`](@ref) of the layer holding the period - `:gsl` for `gsp`. It names a
    layer of the *same* dataset, which is what makes the second read resolvable from the first.
"""
struct PerCellAccumulationPeriod <: AbstractAccumulationPeriod
    code::CODE_TYPE
end
