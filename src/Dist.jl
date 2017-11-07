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
