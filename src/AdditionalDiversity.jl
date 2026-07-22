# SPDX-License-Identifier: LGPL-3.0-or-later

using Diversity

"""
    makeunique(eco::Ecosystem)
Convert type of similarity in [`SpeciesList`](@ref) to UniqueTypes, i.e. an
identity matrix.
"""
function makeunique(eco::Ecosystem)
    sppl = eco.spplist
    spp = length(sppl.names)
    EcoSISTEM.invalidatecaches!(eco)
    newsppl = SpeciesList{typeof(sppl.tolerance),
                          typeof(sppl.demand),
                          typeof(sppl.movement),
                          UniqueTypes,
                          typeof(sppl.params)}(sppl.names,
                                               sppl.tolerance,
                                               sppl.abun,
                                               sppl.demand,
                                               UniqueTypes(spp),
                                               sppl.movement,
                                               sppl.params,
                                               sppl.native)
    newsppl.susceptible = sppl.susceptible
    return Ecosystem{typeof(eco.habitat), typeof(newsppl),
                     typeof(eco.nichefit)}(eco.abundances,
                                           newsppl,
                                           eco.habitat,
                                           eco.ordinariness,
                                           eco.nichefit,
                                           eco.lookup,
                                           eco.cache,
                                           eco.rngs)
end

"""
    meta_simpson(eco::Ecosystem, qs::Vector{Float64})
Calculate the Simpson diversity for the entire ecosystem.
"""
function meta_simpson(eco::Ecosystem, qs::Vector{Float64})
    eco = makeunique(eco)
    div = meta_gamma(eco, 2.0)
    div[!, :diversity] = 1 ./ div[!, :diversity]
    return div
end

function meta_simpson(eco::Ecosystem, qs::Float64)
    eco = makeunique(eco)
    div = meta_gamma(eco, 2.0)
    div[!, :diversity] = 1 ./ div[!, :diversity]
    return div
end

"""
    meta_shannon(eco::Ecosystem, qs::Vector{Float64})
Calculate the Shannon entropy for the entire ecosystem.
"""
function meta_shannon(eco::Ecosystem, qs::Vector{Float64})
    eco = makeunique(eco)
    div = meta_gamma(eco, 1.0)
    div[!, :diversity] = log.(div[!, :diversity])
    return div
end

function meta_shannon(eco::Ecosystem, qs::Float64)
    eco = makeunique(eco)
    div = meta_gamma(eco, 1.0)
    div[!, :diversity] = log.(div[!, :diversity])
    return div
end

"""
    meta_speciesrichness(eco::Ecosystem, qs::Vector{Float64})
Calculate the species richness for the entire ecosystem.
"""
function meta_speciesrichness(eco::Ecosystem, qs::Vector{Float64})
    eco = makeunique(eco)
    return meta_gamma(eco, 0.0)
end

function meta_speciesrichness(eco::Ecosystem, qs::Float64)
    eco = makeunique(eco)
    return meta_gamma(eco, 0.0)
end

"""
    mean_abun(eco::Ecosystem, qs::Vector{Float64})
Calculate the mean arithmetic abundance for the entire ecosystem.
"""
function mean_abun(eco::Ecosystem, qs::Vector{Float64})
    eco = makeunique(eco)
    SR = meta_speciesrichness(eco, 0.0)
    SR[:, :diversity] .= sum(eco.abundances.matrix) ./
                         size(eco.abundances.matrix, 1)
    SR[:, :measure] .= "Mean abundance"
    return SR
end

function mean_abun(eco::Ecosystem, qs::Float64)
    return mean_abun(eco, [qs])
end

"""
    geom_mean_abun(eco::Ecosystem, qs::Vector{Float64})
Calculate the geometric mean abundance for the entire ecosystem.
"""
function geom_mean_abun(eco::Ecosystem, qs::Vector{Float64})
    eco = makeunique(eco)
    SR = meta_speciesrichness(eco, 0.0)
    SR[:, :diversity] .= exp.(sum(log.(mapslices(sum, eco.abundances.matrix,
                                                 dims = 2) .+ 1)) ./
                              size(eco.abundances.matrix, 1)) .- 1
    SR[:, :measure] .= "Geometric mean abundance"
    return SR
end

function geom_mean_abun(eco::Ecosystem, qs::Float64)
    return geom_mean_abun(eco, [qs])
end

"""
    sorenson(eco::Ecosystem, qs::Vector{Float64})
Calculate the Sorenson similarity for the entire ecosystem.
"""
function sorenson(eco::Ecosystem, qs::Vector{Float64})
    eco = makeunique(eco)
    SR = meta_speciesrichness(eco, 0.0)
    ab1 = eco.spplist.abun
    ab2 = mapslices(sum, eco.abundances.matrix, dims = 2)
    SR[:, :diversity] .= 1 - abs(sum(ab1 .- ab2)) / sum(ab1 .+ ab2)
    SR[:, :measure] .= "Sorenson"
    return SR
end

function sorenson(eco::Ecosystem, qs::Float64)
    return sorenson(eco, [qs])
end

"""
    pd(eco::Ecosystem, qs::Vector{Float64})
Calculate Faith's phylogenetic diversity (PD) for the entire ecosystem.
"""
function pd(eco::Ecosystem, qs::Vector{Float64})
    PD = meta_gamma(eco, 0.0)
    PD[:, :diversity] .= PD[:, :diversity] /
                         mean(heightstoroot(eco.spplist.types.tree))
    return PD
end

function pd(eco::Ecosystem, qs::Float64)
    PD = meta_gamma(eco, 0.0)
    PD[:, :diversity] .= PD[:, :diversity] /
                         mean(heightstoroot(eco.spplist.types.tree))
    return PD
end
