using Distributions

function jexp(theta, n::Int64=1)
  rand(Exponential(theta), n)
end
function jpois(gamma, n::Int64=1)
  rand(Poisson(gamma), n)
end

function jbinom(n::Int64, size::Int64, p::Real)
  rand(Binomial(size,p), n)
end

function jmulti(n::Int64, p::AbstractArray)
  rand(Multinomial(n, p))
end
function jmulti(n::Int64, p::Real)
  rand(Multinomial(n, repmat([p], n)))
end


# Function to sample randomly from the Uniform distribution
function junif(a,b)
  rand(Uniform(a,b))
end

# Function to sample randomly from the Dirichlet distribution
function jdir(k,a)
  rand(Dirichlet(k,a))
end

function jnorm(μ,σ,n::Int=1)
  rand(Normal(μ,σ),n)
end
