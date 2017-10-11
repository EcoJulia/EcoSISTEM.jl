using Unitful
import Base.eltype
"""
    AbstractTraitRelationship{TR}

The abstract supertype of relationships between a trait and its environment,
parameterised on any TR.

"""

abstract type AbstractTraitRelationship{TR}
end
"""
    Gauss{TR} <: AbstractTraitRelationship{TR}

The Gaussian relationship between a continuous trait and its environment,
paramaterised on any TR.

"""
mutable struct Gauss{TR} <: AbstractTraitRelationship{TR}
end

function (::Gauss{TR})(current::TR, opt::TR, var::TR) where TR
    return (1.0 * unit(current))/sqrt(2 * π * var^2) * exp(-abs(current - opt)^2/(2 * var^2))
end
iscontinuous(tr::Gauss{TR}) where TR = true
function eltype{TR}(tr::Gauss{TR})
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
function (::Match{TR})(niche::TR, pref::TR) where TR
  if niche == pref
    return 1.0
  else
    return 0.5
  end
end
iscontinuous(tr::Match{TR}) where TR = false
function eltype{TR}(tr::Match{TR})
    return TR
end
"""
    NoRelContinuous{TR} <: AbstractTraitRelationship{TR}

The absense of a relationship between a continuous trait and its environment,
paramaterised on any TR. Returns the value 1.

"""
mutable struct NoRelContinuous{TR} <: AbstractTraitRelationship{TR}
end
function (::NoRelContinuous{TR})(niche::TR, pref::TR) where TR
    return 1.0
end
iscontinuous(tr::NoRelContinuous{TR}) where TR = true
function eltype{TR}(tr::NoRelContinuous{TR})
    return TR
end
"""
    NoRelDiscrete{TR} <: AbstractTraitRelationship{TR}

The absense of a relationship between a discrete trait and its environment,
paramaterised on any TR. Returns the value 1.

"""
mutable struct NoRelDiscrete{TR} <: AbstractTraitRelationship{TR}
end
function (::NoRelDiscrete{TR})(niche::TR, pref::TR) where TR
    return 1.0
end
iscontinuous(tr::NoRelDiscrete{TR}) where TR = true
function eltype{TR}(tr::NoRelDiscrete{TR})
    return TR
end
"""
    multiplicativeTR2{TR1, TR2} <: AbstractTraitRelationship{Tuple{TR1, TR2}}

Type that houses multiple AbstractTraitRelationships for two trait and
habitat levels.

"""
mutable struct multiplicativeTR2{TR1, TR2} <: AbstractTraitRelationship{Tuple{TR1, TR2}}
  tr1::TR1
  tr2::TR2
end
function eltype(mtr::multiplicativeTR2)
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

function eltype(mtr::multiplicativeTR3)
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
