using Diversity
"""
    meta_simpson(eco::Ecosystem, qs::Vector{Float64})
Function to calculate the Simpson diversity for the entire ecosystem.
"""
function meta_simpson(eco::Ecosystem, qs::Vector{Float64})
    div = meta_gamma(eco, 2.0)
    div[:diversity] = 1 ./ div[:diversity]
    return div
end

function meta_simpson(eco::Ecosystem, qs::Float64)
    div = meta_gamma(eco, 2.0)
    div[:diversity] = 1 ./ div[:diversity]
    return div
end
"""
    meta_shannon(eco::Ecosystem, qs::Vector{Float64})
Function to calculate the Shannon entropy for the entire ecosystem.
"""
function meta_shannon(eco::Ecosystem, qs::Vector{Float64})
    div = meta_gamma(eco, 1.0)
    div[:diversity] = log.(div[:diversity])
    return div
end
function meta_shannon(eco::Ecosystem, qs::Float64)
    div = meta_gamma(eco, 1.0)
    div[:diversity] = log(div[:diversity])
    return div
end
"""
    meta_speciesrichness(eco::Ecosystem, qs::Vector{Float64})
Function to calculate the species richness for the entire ecosystem.
"""
function meta_speciesrichness(eco::Ecosystem, qs::Vector{Float64})
    return meta_gamma(eco, 0.0)
end
function meta_speciesrichness(eco::Ecosystem, qs::Float64)
    return meta_gamma(eco, 0.0)
end
