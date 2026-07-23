# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using EcoSISTEM.Units

import Base.eltype
"""
    AbstractNicheFit{TR}

The abstract supertype of relationships between a trait and its environment,
parameterised on any TR.

"""
abstract type AbstractNicheFit{TR} end

"""
    NicheSuitability{TR} <: AbstractNicheFit{TR}

The nichefit between a [`NicheTolerance`](@ref) continuous trait and its environment: the density of the
trait's response distribution evaluated at the current regime value, parameterised on any `TR`.
Works for any `Distributions.ContinuousUnivariateDistribution` (e.g. [`Trapezoid`](@ref) or `Uniform`).
"""
mutable struct NicheSuitability{TR} <: AbstractNicheFit{TR}
end

function (::NicheSuitability{TR})(dist::ContinuousUnivariateDistribution,
                                  current::TR) where {TR}
    # The distribution's parameters are bare numbers in the axis's canonical unit, and regimes are
    # stored in that same unit, so strip `current` to a bare number before evaluating the pdf.
    return pdf(dist, ustrip(current))
end
iscontinuous(tr::NicheSuitability) = true
function Base.eltype(tr::NicheSuitability{TR}) where {TR}
    return TR
end

# The deprecated `Gauss`/`Trapeze`/`Unif` nichefit shims (all `NicheSuitability` now) live in
# `src/deprecations.jl`.

"""
    MatchSuitability{TR} <: AbstractNicheFit{TR}

The nichefit between a discrete trait and its environment,
paramaterised on any TR. Current conditions are matched to a trait preference
and checked for a match.

"""
mutable struct MatchSuitability{TR} <: AbstractNicheFit{TR}
end
function (::MatchSuitability{TR})(niche::TR, pref::TR) where {TR}
    if niche == pref
        return 1.0
    else
        return 0.5
    end
end
iscontinuous(tr::MatchSuitability) = false
function Base.eltype(tr::MatchSuitability{TR}) where {TR}
    return TR
end

mutable struct LandCoverSuitability{TR} <: AbstractNicheFit{TR}
end
function (::LandCoverSuitability{TR})(niche::TR, pref::Vector{TR}) where {TR}
    if niche in pref
        return 1.0
    else
        return 0.0
    end
end
iscontinuous(tr::LandCoverSuitability) = false
function Base.eltype(::LandCoverSuitability{TR}) where {TR}
    return TR
end

"""
    NoFitContinuous{TR} <: AbstractNicheFit{TR}

The absence of a nichefit between a continuous trait and its environment,
paramaterised on any TR. Returns the value 1.

"""
mutable struct NoFitContinuous{TR} <: AbstractNicheFit{TR}
end
function (::NoFitContinuous{TR})(::TR, ::TR, ::TR) where {TR}
    return 1.0
end
iscontinuous(tr::NoFitContinuous) = true
function Base.eltype(tr::NoFitContinuous{TR}) where {TR}
    return TR
end

"""
    NoFitDiscrete{TR} <: AbstractNicheFit{TR}

The absence of a nichefit between a discrete trait and its environment,
paramaterised on any TR. Returns the value 1.

"""
mutable struct NoFitDiscrete{TR} <: AbstractNicheFit{TR}
end
function (::NoFitDiscrete{TR})(niche::TR, pref::TR) where {TR}
    return 1.0
end
iscontinuous(tr::NoFitDiscrete) = false
function Base.eltype(tr::NoFitDiscrete{TR}) where {TR}
    return TR
end

"""
    multiplicativeFit2{TR1, TR2} <: AbstractNicheFit{Tuple{TR1, TR2}}

Type that houses multiple AbstractNicheFits for two trait and
regime levels.

"""
mutable struct multiplicativeFit2{TR1, TR2} <:
               AbstractNicheFit{Tuple{TR1, TR2}}
    one::TR1
    two::TR2
end
function iscontinuous(tr::multiplicativeFit2)
    return [iscontinuous(tr.one), iscontinuous(tr.two)]
end
function Base.eltype(mtr::multiplicativeFit2)
    return [eltype(mtr.one), eltype(mtr.two)]
end

"""
    multiplicativeFit3{TR1, TR2, TR3} <: AbstractNicheFit{Tuple{TR1, TR2, TR3}}

Type that houses multiple AbstractNicheFits for three trait and
regime levels.

"""
mutable struct multiplicativeFit3{TR1, TR2, TR3} <:
               AbstractNicheFit{Tuple{TR1, TR2, TR3}}
    one::TR1
    two::TR2
    three::TR3
end
function iscontinuous(tr::multiplicativeFit3)
    return [iscontinuous(tr.one), iscontinuous(tr.two), iscontinuous(tr.three)]
end
function Base.eltype(mtr::multiplicativeFit3)
    return [eltype(mtr.one), eltype(mtr.two), eltype(mtr.three)]
end

"""
    additiveFit2{TR1, TR2} <: AbstractNicheFit{Tuple{TR1, TR2}}

Type that houses multiple AbstractNicheFits for two trait and
regime levels.

"""
mutable struct additiveFit2{TR1, TR2} <:
               AbstractNicheFit{Tuple{TR1, TR2}}
    one::TR1
    two::TR2
end
function iscontinuous(tr::additiveFit2)
    return [iscontinuous(tr.one), iscontinuous(tr.two)]
end
function Base.eltype(mtr::additiveFit2)
    return [eltype(mtr.one), eltype(mtr.two)]
end

"""
    multiplicativeFit3{TR1, TR2, TR3} <: AbstractNicheFit{Tuple{TR1, TR2, TR3}}

Type that houses multiple AbstractNicheFits for three trait and
regime levels.

"""
mutable struct additiveFit3{TR1, TR2, TR3} <:
               AbstractNicheFit{Tuple{TR1, TR2, TR3}}
    one::TR1
    two::TR2
    three::TR3
end
function iscontinuous(tr::additiveFit3)
    return [iscontinuous(tr.one), iscontinuous(tr.two), iscontinuous(tr.three)]
end
function Base.eltype(mtr::additiveFit3)
    return [eltype(mtr.one), eltype(mtr.two), eltype(mtr.three)]
end

"""
    combinefit

Function that combines the output of multiple niche fits, which varies
depending on whether multiplicative, additive etc.

"""
function combinefit(tr::Union{multiplicativeFit2, multiplicativeFit3})
    return *
end
function combinefit(tr::Union{additiveFit2, additiveFit3})
    return +
end
