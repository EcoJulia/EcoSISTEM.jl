"""
    TraitRelationship

The relationship between a trait and its environment, represented as a Matrix
of Floats.

"""
mutable struct TraitRelationship
  traitfun::Function
end

function GaussTemp(temp::Float64, opttemp::Float64, var::Float64)
  1/sqrt(2 * π * var^2) * exp(-abs(temp - opttemp)^2/(2 * var^2))
end

temp = 5.0
opttemp = 5.0
var = 0.1
GaussTemp(temp, opttemp, var)
exp((-2 * (abs(opttemp - temp))^2) / (var^2)) / (var)
