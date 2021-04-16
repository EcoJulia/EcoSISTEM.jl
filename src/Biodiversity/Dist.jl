using Distributions
import Random: rand
using Random

import Distributions: @check_args, ContinuousUnivariateDistribution,
rand, params, GLOBAL_RNG, pdf

"""
    Trapezoid{T<:Real} <: ContinuousUnivariateDistribution

Trapezoidal distribution as described at https://en.wikipedia.org/wiki/Trapezoidal_distribution.
"""
struct Trapezoid{T<:Real} <: ContinuousUnivariateDistribution
    a::T
    b::T
    c::T
    d::T

    Trapezoid{T}(a::T, b::T, c::T, d::T) where {T} = (@check_args(Trapezoid, a < d); new{T}(a, b, c, d))
end
Trapezoid(a::T, b::T, c::T, d::T) where {T<:Real} = Trapezoid{T}(a, b, c, d)
Trapezoid(a::Real, b::Real, c::Real, d::Real) = Uniform(promote(a, b, c, d)...)
Trapezoid(a::Integer, b::Integer, c::Integer, d::Integer) = Trapezoid(Float64(a), Float64(b),
    Float64(c), Float64(d))
Trapezoid() = Trapezoid(0.0, 0.0, 1.0, 1.0)

params(d::Trapezoid) = (d.a, d.b, d.c, d.d)
rand(d::Trapezoid) = rand(GLOBAL_RNG, d)
function rand(rng::AbstractRNG, T::Trapezoid)
    (a, b, c, d) = params(T)
    b_m_a = b - a
    c_m_b = c - b
    d_m_c = d - c
    Cϕ = 4 / ((b_m_a * 2) + (c_m_b * 4) + (d_m_c * 2))
    pi1 = Cϕ * b_m_a/2
    pi2 = Cϕ * c_m_b
    pi3 = Cϕ * d_m_c/2
    u = rand(rng)
    if (u >= 0 && u <= pi1)
        return a + (u/pi1)^(1/2) * b_m_a
    elseif (u > pi1 && u <= (1 - pi3))
        return b + ((u - pi1)/pi2) * c_m_b
    else
        return d - ((1 - u)/pi3)^(1/2) * d_m_c
    end
end

function pdf(T::Trapezoid, x::Real)
    (a, b, c, d) = params(T)
    d_p_c = d + c
    a_p_b = a + b
    b_m_a = b - a
    d_m_c = d - c
    if (x >= a && x < b)
        return (2 / (d_p_c- a_p_b)) * ((x-a)/b_m_a)
    elseif (x >= b && x < c)
        return (2 / (d_p_c- a_p_b))
    elseif (x >=c && x <= d)
        return (2 / (d_p_c- a_p_b)) * ((d-x)/d_m_c)
    else
        return 0.0
    end
end
