# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using EcoSISTEM.Units

import Base.eltype
"""
    AbstractNicheFit{NF}

The abstract supertype of relationships between a trait and its environment,
parameterised on any NF.

"""
abstract type AbstractNicheFit{NF} end

"""
    NicheSuitability{NF} <: AbstractNicheFit{NF}

The nichefit between a [`NicheTolerance`](@ref) continuous trait and its environment: the density of the
trait's response distribution evaluated at the current regime value, parameterised on any `NF`.
Works for any `Distributions.ContinuousUnivariateDistribution` (e.g. [`Trapezoid`](@ref) or `Uniform`).
"""
mutable struct NicheSuitability{NF} <: AbstractNicheFit{NF}
end

function (::NicheSuitability{NF})(dist::ContinuousUnivariateDistribution,
                                  current) where {NF}
    return pdf(dist, _toframe(NF, current))
end

# Convert `current` to the frame `NF` expects before stripping. When `NF` is a concrete `Quantity` type
# (the real case, via `_default_suitability`'s tolerance-derived `NF`), this is a real `uconvert` — so a
# dimensionally-compatible-but-differently-scaled regime is corrected rather than silently wrong, and a
# genuinely incompatible dimension still throws, just via `uconvert` instead of dispatch. `NF` can also be
# an abstract dimension (e.g. `Unitful.Temperature`) or a bare non-quantity type — neither names a single
# fixed target scale to convert to, so just strip as before.
_toframe(::Type{NF}, current) where {NF <: Quantity} = ustrip(unit(NF), current)
_toframe(::Type{NF}, current) where {NF} = ustrip(current)
iscontinuous(nichefit::NicheSuitability) = true
function Base.eltype(nichefit::NicheSuitability{NF}) where {NF}
    return NF
end

# The deprecated `Gauss`/`Trapeze`/`Unif` nichefit shims (all `NicheSuitability` now) live in
# `src/deprecations.jl`.

"""
    MatchSuitability{NF} <: AbstractNicheFit{NF}

The nichefit between a discrete trait and its environment,
paramaterised on any NF. Current conditions are matched to a trait preference
and checked for a match.

"""
mutable struct MatchSuitability{NF} <: AbstractNicheFit{NF}
end
function (::MatchSuitability{NF})(niche::NF, pref::NF) where {NF}
    if niche == pref
        return 1.0
    else
        return 0.5
    end
end
iscontinuous(nichefit::MatchSuitability) = false
function Base.eltype(nichefit::MatchSuitability{NF}) where {NF}
    return NF
end

mutable struct LandCoverSuitability{NF} <: AbstractNicheFit{NF}
end
function (::LandCoverSuitability{NF})(niche::NF, pref::Vector{NF}) where {NF}
    if niche in pref
        return 1.0
    else
        return 0.0
    end
end
iscontinuous(nichefit::LandCoverSuitability) = false
function Base.eltype(::LandCoverSuitability{NF}) where {NF}
    return NF
end

"""
    NoFitContinuous{NF} <: AbstractNicheFit{NF}

The absence of a nichefit between a continuous trait and its environment,
paramaterised on any NF. Returns the value 1.

"""
mutable struct NoFitContinuous{NF} <: AbstractNicheFit{NF}
end
function (::NoFitContinuous{NF})(::NF, ::NF, ::NF) where {NF}
    return 1.0
end
iscontinuous(nichefit::NoFitContinuous) = true
function Base.eltype(nichefit::NoFitContinuous{NF}) where {NF}
    return NF
end

"""
    NoFitDiscrete{NF} <: AbstractNicheFit{NF}

The absence of a nichefit between a discrete trait and its environment,
paramaterised on any NF. Returns the value 1.

"""
mutable struct NoFitDiscrete{NF} <: AbstractNicheFit{NF}
end
function (::NoFitDiscrete{NF})(niche::NF, pref::NF) where {NF}
    return 1.0
end
iscontinuous(nichefit::NoFitDiscrete) = false
function Base.eltype(nichefit::NoFitDiscrete{NF}) where {NF}
    return NF
end

"""
    multiplicativeFit2{NF1, NF2} <: AbstractNicheFit{Tuple{NF1, NF2}}

Type that houses multiple AbstractNicheFits for two trait and
regime levels.

"""
mutable struct multiplicativeFit2{NF1, NF2} <:
               AbstractNicheFit{Tuple{NF1, NF2}}
    one::NF1
    two::NF2
end
function iscontinuous(nichefit::multiplicativeFit2)
    return [iscontinuous(nichefit.one), iscontinuous(nichefit.two)]
end
function Base.eltype(mnichefit::multiplicativeFit2)
    return [eltype(mnichefit.one), eltype(mnichefit.two)]
end

"""
    multiplicativeFit3{NF1, NF2, NF3} <: AbstractNicheFit{Tuple{NF1, NF2, NF3}}

Type that houses multiple AbstractNicheFits for three trait and
regime levels.

"""
mutable struct multiplicativeFit3{NF1, NF2, NF3} <:
               AbstractNicheFit{Tuple{NF1, NF2, NF3}}
    one::NF1
    two::NF2
    three::NF3
end
function iscontinuous(nichefit::multiplicativeFit3)
    return [
        iscontinuous(nichefit.one),
        iscontinuous(nichefit.two),
        iscontinuous(nichefit.three)
    ]
end
function Base.eltype(mnichefit::multiplicativeFit3)
    return [
        eltype(mnichefit.one),
        eltype(mnichefit.two),
        eltype(mnichefit.three)
    ]
end

"""
    additiveFit2{NF1, NF2} <: AbstractNicheFit{Tuple{NF1, NF2}}

Type that houses multiple AbstractNicheFits for two trait and
regime levels.

"""
mutable struct additiveFit2{NF1, NF2} <:
               AbstractNicheFit{Tuple{NF1, NF2}}
    one::NF1
    two::NF2
end
function iscontinuous(nichefit::additiveFit2)
    return [iscontinuous(nichefit.one), iscontinuous(nichefit.two)]
end
function Base.eltype(mnichefit::additiveFit2)
    return [eltype(mnichefit.one), eltype(mnichefit.two)]
end

"""
    multiplicativeFit3{NF1, NF2, NF3} <: AbstractNicheFit{Tuple{NF1, NF2, NF3}}

Type that houses multiple AbstractNicheFits for three trait and
regime levels.

"""
mutable struct additiveFit3{NF1, NF2, NF3} <:
               AbstractNicheFit{Tuple{NF1, NF2, NF3}}
    one::NF1
    two::NF2
    three::NF3
end
function iscontinuous(nichefit::additiveFit3)
    return [
        iscontinuous(nichefit.one),
        iscontinuous(nichefit.two),
        iscontinuous(nichefit.three)
    ]
end
function Base.eltype(mnichefit::additiveFit3)
    return [
        eltype(mnichefit.one),
        eltype(mnichefit.two),
        eltype(mnichefit.three)
    ]
end

"""
    combinefit

Function that combines the output of multiple niche fits, which varies
depending on whether multiplicative, additive etc.

"""
function combinefit(nichefit::Union{multiplicativeFit2, multiplicativeFit3})
    return *
end
function combinefit(nichefit::Union{additiveFit2, additiveFit3})
    return +
end
