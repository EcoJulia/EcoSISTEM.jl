"""
    TraitRelationship

The relationship between a trait and its environment, represented as a Matrix
of Floats.

"""
mutable struct TraitRelationship
  traitfun::Function
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
