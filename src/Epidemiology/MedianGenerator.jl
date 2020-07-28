struct MedianGenerator <: Random.AbstractRNG end

MedianGenerator(args...) = MedianGenerator()

import Random: rand, rand!
import Distributions: rand, rand!, median

function median(dist::Multinomial)
    values = dist.n .* dist.p
    rounded = round.(Int64, values)
    values .-= rounded
    order = sortperm(values)
    count = sum(rounded) - dist.n

    n = 1
    while count > 0
        if rounded[order[n]] > 0
            rounded[order[n]] -= 1
            count -= 1
        end
        n += 1
    end

    if count < 0
        rounded[order[(length(order)+count-1):length(order)]] .+= 1
    end

    return rounded
end

rand(::MedianGenerator, dist::Distribution{Univariate,S}) where {S<:ValueSupport} = median(dist)
rand(::MedianGenerator, dist::Binomial) = median(dist)

function rand!(::MedianGenerator, dist::Distribution{Multivariate,S},
        container::AbstractVector{T}) where {S<:ValueSupport, T<:Integer}
    container .= T.(median(dist))
    return nothing
end

