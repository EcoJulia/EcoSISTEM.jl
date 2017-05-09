"""
    AbstractMovement

Abstract supertype of movements
"""
abstract AbstractMovement
"""
    GaussianMovement <: AbstractMovement

GaussianMovement holds parameters for a gaussian movement kernel; a vector of
dispersal variances per species, `var`, and a threshold, `thresh`, beyond which
dispersal cannot take place.
"""
type GaussianMovement <: AbstractMovement
  var::Vector{Float64}
  thresh::Float64
end

function GaussianMovement(dispersalvar::Float64, numspecies::Int64,
  pthresh::Float64)
  GaussianMovement(repmat([dispersalvar], numspecies), pthresh)
end
