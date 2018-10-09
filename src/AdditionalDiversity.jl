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


function mean_abun(eco::Ecosystem, qs::Vector{Float64})
    SR = meta_speciesrichness(eco, 0.0)
    SR[:diversity] = sum(eco.abundances.matrix) ./ SR[:diversity]
    SR[:measure] = "Mean abundance"
    return SR
end
function geom_mean_abun(eco::Ecosystem, qs::Float64)
    SR = meta_speciesrichness(eco, 0.0)
    SR[:diversity] = exp(sum(log.(mapslices(sum, eco.abundances.matrix, 1))) ./
                        SR[:diversity])
    SR[:measure] = "Geometric mean abundance"
    return SR
end
