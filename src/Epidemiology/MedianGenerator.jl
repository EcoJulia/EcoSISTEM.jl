struct MedianGenerator <: Random.AbstractRNG end

MedianGenerator(args...) = MedianGenerator()

import Random: rand, rand!
import Distributions: rand, rand!, median

function median(dist::Multinomial)
  v = dist.n * dist.p
  output = Int64.(floor.(v))
  sum(output) == dist.n && return output
  @assert dist.n - length(dist.p) <= sum(output) < dist.n
  p = sortperm(v .- output)
  for rp in reverse(p)
    sum(output) == dist.n && break
    output[rp] += 1
  end
  @assert sum(output) == dist.n
  return output
end

# Distributions/src/multivariates.jl exports rand(::AbstractRNG, _) which shouldn't be called
# with rand(::MedianGenerator, _)
function rand(::MedianGenerator, dist)
  throw(ArgumentError("rand(::MedianGenerator, dist) not implemented for dist=$dist."))
end
rand(::MedianGenerator, dist::Distribution{Univariate,S}) where {S<:ValueSupport} = median(dist)
rand(::MedianGenerator, dist::Binomial) = median(dist)
rand(::MedianGenerator, dist::Distribution{Multivariate,S}) where {S<:ValueSupport} = median(dist)

function rand!(::MedianGenerator, dist::Distribution{Multivariate,S},
        container::AbstractVector{T}) where {S<:ValueSupport, T<:Integer}
    container .= median(dist)
    return nothing
end
