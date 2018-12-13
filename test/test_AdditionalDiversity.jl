using Simulation
using Diversity
using Phylo
using Statistics
using Test

eco = TestEcosystem()

Ab = mapslices(sum, eco.abundances.matrix, dims = 2)
RelAb = Ab / sum(eco.abundances.matrix)
## TEST makeunique
@test_nowarn makeunique(eco)
@test typeof(makeunique(eco).spplist.types) <: UniqueTypes

## TEST meta_simpson
@test meta_simpson(eco, 1.0)[:diversity] == meta_simpson(eco, 2.0)[:diversity]
@test meta_simpson(eco, 1.0)[:diversity][1] ≈
    sum(RelAb .^ 2)
@test meta_simpson(eco, 1.0)[:diversity] == 1 ./ meta_gamma(makeunique(eco), 2.0)[:diversity]
@test meta_simpson(eco, 2.0) == meta_simpson(makeunique(eco), 2.0)

## TEST meta_simpson
@test meta_shannon(eco, 1.0)[:diversity] == meta_shannon(eco, 2.0)[:diversity]
@test meta_shannon(eco, 1.0)[:diversity][1] ≈ -sum(RelAb .* log.(RelAb))
@test meta_shannon(eco, 1.0)[:diversity] == log.(meta_gamma(makeunique(eco), 1.0)[:diversity])
@test meta_shannon(eco, 2.0) == meta_shannon(makeunique(eco), 2.0)

## TEST meta_speciesrichness
@test meta_speciesrichness(eco, 1.0)[:diversity] == meta_speciesrichness(eco, 0.0)[:diversity]
@test meta_speciesrichness(eco, 1.0)[:diversity][1] ≈ sum(Ab .> 0)
@test meta_speciesrichness(eco, 1.0)[:diversity] == meta_gamma(makeunique(eco), 0.0)[:diversity]
@test meta_speciesrichness(eco, 0.0) == meta_speciesrichness(makeunique(eco), 0.0)

## TEST mean_abun
@test mean_abun(eco, 1.0)[:diversity] == mean_abun(eco, 0.0)[:diversity]
@test mean_abun(eco, 1.0)[:diversity][1] ≈sum(eco.abundances.matrix) ./
    size(eco.abundances.matrix, 1)
@test mean_abun(eco, 0.0) == mean_abun(makeunique(eco), 0.0)
@test mean_abun(eco, 0.0) == mean_abun(eco, [0.0, 1, 2])

## TEST geom_mean_abun
@test geom_mean_abun(eco, 1.0)[:diversity] == geom_mean_abun(eco, 0.0)[:diversity]
@test geom_mean_abun(eco, 1.0)[:diversity][1] ≈ exp.(sum(log.(mapslices(sum, eco.abundances.matrix, dims = 2) .+ 1)) ./
                    size(eco.abundances.matrix, 1)) .- 1
@test geom_mean_abun(eco, 0.0) == geom_mean_abun(makeunique(eco), 0.0)
@test geom_mean_abun(eco, 0.0) == geom_mean_abun(eco, [0.0, 1, 2])

## TEST sorenson
@test sorenson(eco, 1.0)[:diversity] == sorenson(eco, 0.0)[:diversity]
@test sorenson(eco, 0.0) == sorenson(makeunique(eco), 0.0)
@test sorenson(eco, 0.0) == sorenson(eco, [0.0, 1, 2])


## TEST pd (against R and julia)
@test pd(eco, 1.0)[:diversity] == pd(eco, 0.0)[:diversity]
tree = eco.spplist.types.tree
mat = reshape(mapslices(sum, eco.abundances.matrix, dims = 2), 1, 150)
@rput tree
@rput mat
R"library(picante)
       mat = data.frame(mat, row.names ='eco1')
       colnames(mat) = tree$tip.label
       p = pd(mat, tree, include.root=F)"
@test pd(eco, 1.0)[:diversity] ≈ @rget(p)[:PD] / mean(heightstoroot(eco.spplist.types.tree))
@test pd(eco, 1.0)[:diversity][1] ≈ sum(map(b -> getlength(eco.spplist.types.tree,b),
           getbranchnames(eco.spplist.types.tree))) / mean(heightstoroot(eco.spplist.types.tree))
@test pd(eco, 0.0) == pd(eco, [0.0, 1, 2])
