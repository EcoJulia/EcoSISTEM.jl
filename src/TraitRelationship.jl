# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using EcoSISTEM.Units

import Base.eltype
"""
    AbstractTraitRelationship{TR}

The abstract supertype of relationships between a trait and its environment,
parameterised on any TR.

"""
abstract type AbstractTraitRelationship{TR} end

"""
    DistRel{TR} <: AbstractTraitRelationship{TR}

The relationship between a [`NicheTolerance`](@ref) continuous trait and its environment: the density of the
trait's response distribution evaluated at the current regime value, parameterised on any `TR`.
Works for any `Distributions.ContinuousUnivariateDistribution` (e.g. [`Trapezoid`](@ref) or `Uniform`).
"""
mutable struct DistRel{TR} <: AbstractTraitRelationship{TR}
end

function (::DistRel{TR})(dist::ContinuousUnivariateDistribution,
                         current::TR) where {TR}
    # The distribution's parameters are bare numbers in the axis's canonical unit, and regimes are
    # stored in that same unit, so strip `current` to a bare number before evaluating the pdf.
    return pdf(dist, ustrip(current))
end
iscontinuous(tr::DistRel) = true
function Base.eltype(tr::DistRel{TR}) where {TR}
    return TR
end

# The deprecated `Gauss`/`Trapeze`/`Unif` relationship shims (all `DistRel` now) live in
# `src/deprecations.jl`.

"""
    Match{TR} <: AbstractTraitRelationship{TR}

The relationship between a discrete trait and its environment,
paramaterised on any TR. Current conditions are matched to a trait preference
and checked for a match.

"""
mutable struct Match{TR} <: AbstractTraitRelationship{TR}
end
function (::Match{TR})(niche::TR, pref::TR) where {TR}
    if niche == pref
        return 1.0
    else
        return 0.5
    end
end
iscontinuous(tr::Match) = false
function Base.eltype(tr::Match{TR}) where {TR}
    return TR
end

mutable struct LCmatch{TR} <: AbstractTraitRelationship{TR}
end
function (::LCmatch{TR})(niche::TR, pref::Vector{TR}) where {TR}
    if niche in pref
        return 1.0
    else
        return 0.0
    end
end
iscontinuous(tr::LCmatch) = false
function Base.eltype(::LCmatch{TR}) where {TR}
    return TR
end

"""
    NoRelContinuous{TR} <: AbstractTraitRelationship{TR}

The absence of a relationship between a continuous trait and its environment,
paramaterised on any TR. Returns the value 1.

"""
mutable struct NoRelContinuous{TR} <: AbstractTraitRelationship{TR}
end
function (::NoRelContinuous{TR})(::TR, ::TR, ::TR) where {TR}
    return 1.0
end
iscontinuous(tr::NoRelContinuous) = true
function Base.eltype(tr::NoRelContinuous{TR}) where {TR}
    return TR
end

"""
    NoRelDiscrete{TR} <: AbstractTraitRelationship{TR}

The absence of a relationship between a discrete trait and its environment,
paramaterised on any TR. Returns the value 1.

"""
mutable struct NoRelDiscrete{TR} <: AbstractTraitRelationship{TR}
end
function (::NoRelDiscrete{TR})(niche::TR, pref::TR) where {TR}
    return 1.0
end
iscontinuous(tr::NoRelDiscrete) = false
function Base.eltype(tr::NoRelDiscrete{TR}) where {TR}
    return TR
end

"""
    multiplicativeTR2{TR1, TR2} <: AbstractTraitRelationship{Tuple{TR1, TR2}}

Type that houses multiple AbstractTraitRelationships for two trait and
regime levels.

"""
mutable struct multiplicativeTR2{TR1, TR2} <:
               AbstractTraitRelationship{Tuple{TR1, TR2}}
    one::TR1
    two::TR2
end
function iscontinuous(tr::multiplicativeTR2)
    return [iscontinuous(tr.one), iscontinuous(tr.two)]
end
function Base.eltype(mtr::multiplicativeTR2)
    return [eltype(mtr.one), eltype(mtr.two)]
end

"""
    multiplicativeTR3{TR1, TR2, TR3} <: AbstractTraitRelationship{Tuple{TR1, TR2, TR3}}

Type that houses multiple AbstractTraitRelationships for three trait and
regime levels.

"""
mutable struct multiplicativeTR3{TR1, TR2, TR3} <:
               AbstractTraitRelationship{Tuple{TR1, TR2, TR3}}
    one::TR1
    two::TR2
    three::TR3
end
function iscontinuous(tr::multiplicativeTR3)
    return [iscontinuous(tr.one), iscontinuous(tr.two), iscontinuous(tr.three)]
end
function Base.eltype(mtr::multiplicativeTR3)
    return [eltype(mtr.one), eltype(mtr.two), eltype(mtr.three)]
end

"""
    additiveTR2{TR1, TR2} <: AbstractTraitRelationship{Tuple{TR1, TR2}}

Type that houses multiple AbstractTraitRelationships for two trait and
regime levels.

"""
mutable struct additiveTR2{TR1, TR2} <:
               AbstractTraitRelationship{Tuple{TR1, TR2}}
    one::TR1
    two::TR2
end
function iscontinuous(tr::additiveTR2)
    return [iscontinuous(tr.one), iscontinuous(tr.two)]
end
function Base.eltype(mtr::additiveTR2)
    return [eltype(mtr.one), eltype(mtr.two)]
end

"""
    multiplicativeTR3{TR1, TR2, TR3} <: AbstractTraitRelationship{Tuple{TR1, TR2, TR3}}

Type that houses multiple AbstractTraitRelationships for three trait and
regime levels.

"""
mutable struct additiveTR3{TR1, TR2, TR3} <:
               AbstractTraitRelationship{Tuple{TR1, TR2, TR3}}
    one::TR1
    two::TR2
    three::TR3
end
function iscontinuous(tr::additiveTR3)
    return [iscontinuous(tr.one), iscontinuous(tr.two), iscontinuous(tr.three)]
end
function Base.eltype(mtr::additiveTR3)
    return [eltype(mtr.one), eltype(mtr.two), eltype(mtr.three)]
end

"""
    combineTR

Function that combines the output of multiple trait relationships, which varies
depending on whether multiplicative, additive etc.

"""
function combineTR(tr::Union{multiplicativeTR2, multiplicativeTR3})
    return *
end
function combineTR(tr::Union{additiveTR2, additiveTR3})
    return +
end
