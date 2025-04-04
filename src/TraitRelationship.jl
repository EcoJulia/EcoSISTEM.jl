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
    Gauss{TR} <: AbstractTraitRelationship{TR}

The Gaussian relationship between a continuous trait and its environment,
paramaterised on any TR.

"""
mutable struct Gauss{TR} <: AbstractTraitRelationship{TR}
end

function (::Gauss{TR})(current::TR, opt::TR, var::TR) where {TR}
    pref = (1.0) / sqrt(2 * π * var^2) *
           exp(-abs(current - opt)^2 / (2 * var^2))
    return pref * unit(current)
end
iscontinuous(tr::Gauss) = true
function Base.eltype(tr::Gauss{TR}) where {TR}
    return TR
end

"""
    Trapeze{TR} <: AbstractTraitRelationship{TR}

The relationship between a continuous trait and its environment,
paramaterised on any TR.

"""
mutable struct Trapeze{TR} <: AbstractTraitRelationship{TR}
end

function (::Trapeze{TR})(dist::Trapezoid, current::TR) where {TR}
    return pdf(dist, ustrip(current))
end
iscontinuous(tr::Trapeze) = true
function Base.eltype(tr::Trapeze{TR}) where {TR}
    return TR
end

"""
    Trapeze{TR} <: AbstractTraitRelationship{TR}

The relationship between a continuous trait and its environment,
paramaterised on any TR.

"""
mutable struct Unif{TR} <: AbstractTraitRelationship{TR}
end

function (::Unif{TR})(dist::Uniform, current::TR) where {TR}
    return pdf(dist, uconvert(NoUnits, current / mm))
end
iscontinuous(tr::Unif) = true
function Base.eltype(tr::Unif{TR}) where {TR}
    return TR
end

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

The absense of a relationship between a continuous trait and its environment,
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

The absense of a relationship between a discrete trait and its environment,
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
habitat levels.

"""
mutable struct multiplicativeTR2{TR1, TR2} <:
               AbstractTraitRelationship{Tuple{TR1, TR2}}
    tr1::TR1
    tr2::TR2
end
function iscontinuous(tr::multiplicativeTR2)
    return [iscontinuous(tr.tr1), iscontinuous(tr.tr2)]
end
function Base.eltype(mtr::multiplicativeTR2)
    return [eltype(mtr.tr1), eltype(mtr.tr2)]
end

"""
    multiplicativeTR3{TR1, TR2, TR3} <: AbstractTraitRelationship{Tuple{TR1, TR2, TR3}}

Type that houses multiple AbstractTraitRelationships for three trait and
habitat levels.

"""
mutable struct multiplicativeTR3{TR1, TR2, TR3} <:
               AbstractTraitRelationship{Tuple{TR1, TR2, TR3}}
    tr1::TR1
    tr2::TR2
    tr3::TR3
end
function iscontinuous(tr::multiplicativeTR3)
    return [iscontinuous(tr.tr1), iscontinuous(tr.tr2), iscontinuous(tr.tr3)]
end
function Base.eltype(mtr::multiplicativeTR3)
    return [eltype(mtr.tr1), eltype(mtr.tr2), eltype(mtr.tr3)]
end

"""
    additiveTR2{TR1, TR2} <: AbstractTraitRelationship{Tuple{TR1, TR2}}

Type that houses multiple AbstractTraitRelationships for two trait and
habitat levels.

"""
mutable struct additiveTR2{TR1, TR2} <:
               AbstractTraitRelationship{Tuple{TR1, TR2}}
    tr1::TR1
    tr2::TR2
end
function iscontinuous(tr::additiveTR2)
    return [iscontinuous(tr.tr1), iscontinuous(tr.tr2)]
end
function Base.eltype(mtr::additiveTR2)
    return [eltype(mtr.tr1), eltype(mtr.tr2)]
end

"""
    multiplicativeTR3{TR1, TR2, TR3} <: AbstractTraitRelationship{Tuple{TR1, TR2, TR3}}

Type that houses multiple AbstractTraitRelationships for three trait and
habitat levels.

"""
mutable struct additiveTR3{TR1, TR2, TR3} <:
               AbstractTraitRelationship{Tuple{TR1, TR2, TR3}}
    tr1::TR1
    tr2::TR2
    tr3::TR3
end
function iscontinuous(tr::additiveTR3)
    return [iscontinuous(tr.tr1), iscontinuous(tr.tr2), iscontinuous(tr.tr3)]
end
function Base.eltype(mtr::additiveTR3)
    return [eltype(mtr.tr1), eltype(mtr.tr2), eltype(mtr.tr3)]
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
