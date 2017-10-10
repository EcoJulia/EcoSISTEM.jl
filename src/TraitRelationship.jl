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

function GaussFunc(current::Float64, opt::Float64, var::Float64)
    1.0 /sqrt(2 * π * var^2) * exp(-abs(current - opt)^2/(2 * var^2))
end

function GaussFunc{T <: Unitful.Quantity{Float64}}(current::T,
      opt::T, var::T)
    (1.0 * unit(current))/sqrt(2 * π * var^2) * exp(-abs(current - opt)^2/(2 * var^2))
end

function SimpleDiscrete(niche::Int64, pref::Int64)
  if niche == pref
    return 1.0
  else
    return 0.5
  end
end
function NoRel(niche::Int64, pref::Int64)
    return 1.0
end
function NoRel{T <: Union{Float64, Unitful.Quantity{Float64}}}(val::T,
     optval::T, var::T)
    return 1.0
end
