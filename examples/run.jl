using Simulation
using Diversity

habitat=Habitat([1.0, 2.0])
meta=Metacommunity([0.2 0.2 0.3; 0.3 0.0 0.0]')
loc=Location(habitat, meta)

mat=reshape([loc, loc, loc, loc],(2,2))
landscape=MatrixLandscape(mat)
ecosystem=Ecosystem(landscape)


function calc_diversity(eco::Ecosystem, diversityFn, qs)
  ls = get_landscape(eco)
  map(elt -> diversityFn(elt.metacommunity, qs), ls)
end

function calc_diversity{Eco <: Ecosystem, V <: AbstractVector{Float64}}(eco::Eco, diversityFn, qs::V)
  ls = get_landscape(eco)
  map(elt -> diversityFn(elt.metacommunity, qs)::V, ls)
end

function calc_diversity{Eco <: Ecosystem, FP <: AbstractFloat}(eco::Eco, diversityFn, q::FP)
  ls = get_landscape(eco)
  map(elt -> diversityFn(elt.metacommunity, q)::FP, ls)
end


@code_warntype calc_diversity(ecosystem, metacommunityAbar, [0.0, 1.0, 2.0])
@code_warntype calc_diversity(ecosystem, metacommunityAbar, [0, 1, 2])
@code_warntype calc_diversity(ecosystem, metacommunityAbar, 0.0)
@code_warntype calc_diversity(ecosystem, metacommunityAbar, 0)
calc_diversity(ecosystem, metacommunityAbar, [0.0, 1.0, 2.0])
calc_diversity(ecosystem, metacommunityAbar, [0, 1, 2])
