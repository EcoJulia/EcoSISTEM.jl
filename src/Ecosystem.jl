using Diversity
abstract AbstractHabitat

type Habitat <: AbstractHabitat
  environment::Vector{Float64}
end

type Location{H<:AbstractHabitat}
  habitat::H
  metacommunity::Metacommunity
end

abstract AbstractLandscape

type MatrixLandscape{Loc <: Location} <: AbstractLandscape
  matrix::Matrix{Loc}
end

type Ecosystem{Land<:AbstractLandscape}
  landscape::Land
end

function get_landscape(eco::Ecosystem)
  eco.landscape
end

import Base.start, Base.next, Base.done, Base.eltype, Base.length
function start(land::MatrixLandscape)
  start(land.matrix)
end
function next(land::MatrixLandscape, state)
  next(land.matrix, state)
end
function done(land::MatrixLandscape, state)
  done(land.matrix, state)
end
function eltype(land::MatrixLandscape)
  eltype(land.matrix)
end
function length(land::MatrixLandscape)
  length(land.matrix)
end
