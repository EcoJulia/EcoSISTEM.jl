using Distributions

"""
    jexp(theta, n::Int64=1)

Function to sample randomly from the Exponential distribution with rate theta.
Number sampled is one unless otherwise specified.
"""
function jexp(theta, n::Int64=1)
  rand(Exponential(theta), n)
end

"""
    jpois(lambda, n::Int64=1)

Function to sample randomly from the Poisson distribution with mean lambda.
Number sampled is one unless otherwise specified.
"""
function jpois(lambda, n::Int64=1)
  rand(Poisson(lambda), n)
end

"""
    jbinom(n::Int64, size::Int64, p::Real)

Function to sample randomly from the Binomial distribution with size `size`,
probability `p`, a specified number of times (`n`).
"""
function jbinom(n::Int64, size::Int64, p::Real)
  rand(Binomial(size, p), n)
end

"""
    jmulti(n::Int64, p::AbstractArray)

Function to sample randomly from the Multinomial distribution with size `size`,
probability `p`, a specified number of times (`n`).
"""
function jmulti(n::Int64, p::AbstractArray)
  rand(Multinomial(n, p))
end
function jmulti(n::Int64,size::Int64, p::Real)
  rand(Multinomial(size, repmat([p], n)))
end


"""
    junif(a,b)

Function to sample randomly from the Uniform distribution within bounds `a` and
`b`.
"""
function junif(a::Float64, b::Float64)
  rand(Uniform(a, b))
end

function junif(a::Float64, b::Float64, n::Int64)
  rand(Uniform(a,b), n)
end

function junif(a::Int64, b::Int64)
  rand(Uniform(a, b))
end

function junif(a::Int64, b::Int64, n::Int64)
  rand(Uniform(a,b), n)
end
"""
    jdir(k,a)

Function to sample randomly from the Dirichlet distribution with number of
categories k and positive scalar `a`.
"""
function jdir(k::Int64, a::Int64)
  rand(Dirichlet(k,a))
end

"""
    jnorm(μ::Float64, σ::Float64, n::Int=1)

Function to sample randomly from the Normal distribution with mean μ and
variance σ. Number sampled (`n`) is one unless otherwise specified.
"""
function jnorm(μ,σ,n::Int=1)
  rand(Normal(μ,σ),n)
end
"""
    tnorm(μ, σ, l = 0, u = 1, n::Int = 1)

Function to sample randomly from a truncated Normal distribution with mean μ and
variance σ. `l` and `u` are the upper and lower bounds of truncation, respectively.
Number sampled (`n`) is one unless otherwise specified.
"""
function tnorm(μ, σ, l = 0, u = 1, n::Int = 1)
  rand(Truncated(Normal(μ, σ), l, u), n)
end
"""
    jgamma(shape, scale)

Function to sample randomly from a Gamma distribution with shape and scale
parameters (`shape` and `scale`, respectively).
"""
function jgamma(shape::Float64, scale::Float64)
  rand(Gamma(shape, scale))
end

import Distributions: @check_args, ContinuousUnivariateDistribution,
rand, params, GLOBAL_RNG, pdf
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
