using Unitful
import Base.eltype
"""
    TraitRelationship

The relationship between a trait and its environment, represented as a Matrix
of Floats.

"""
abstract type AbstractTraitRelationship{TR}
end

mutable struct TraitRelationship{TR} <: AbstractTraitRelationship{TR}
    func::Function
end
function eltype{TR}(trait::TraitRelationship{TR})
    return TR
end

mutable struct multiplicativeTR2{TR1, TR2} <: AbstractTraitRelationship{Tuple{TR1, TR2}}
  tr1::TR1
  tr2::TR2
end
function eltype(mtr::multiplicativeTR2)
    return [eltype(mtr.tr1), eltype(mtr.tr2)]
end
mutable struct multiplicativeTR3{TR1, TR2, TR3} <:
    AbstractTraitRelationship{Tuple{TR1, TR2, TR3}}
  tr1::TR1
  tr2::TR2
  tr3::TR3
end
function eltype(mtr::multiplicativeTR3)
    return [eltype(mtr.tr1), eltype(mtr.tr2), eltype(mtr.tr3)]
end
function combineTR(tr::Union{multiplicativeTR2, multiplicativeTR3})
    return *
end

function GaussTemp(temp::Unitful.Temperature{Float64}, opttemp::Unitful.Temperature{Float64},
  var::Unitful.Temperature{Float64})
  1°C/sqrt(2 * π * var^2) * exp(-abs(temp - opttemp)^2/(2 * var^2))
end

function SimpleNiche(niche::String, pref::String)
  if niche == pref
    return 0.5
  else
    return 1.0
  end
end
