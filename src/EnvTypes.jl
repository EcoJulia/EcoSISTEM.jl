
abstract type Envtype end

mutable struct Temp <: Envtype
end

mutable struct Niches <: Envtype
end

mutable struct None <: Envtype
end
