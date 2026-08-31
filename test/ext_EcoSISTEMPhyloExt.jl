# SPDX-License-Identifier: LGPL-3.0-or-later

module TestEcoSISTEMPhyloExt

using Phylo
using DataFrames
using EcoSISTEM
# `[C7-VIS]` B1/B2/B3: these are `public` rather than exported, so they must be named.
using EcoSISTEM: gettraits, reroot!, assigntraits!, resettraits!
using EcoSISTEM: varcovar, fitbrownian
using LinearAlgebra: diag, issymmetric
using Random
using Test

@testset "Phylo traits" begin
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
    reroot!(tree, "tip 2")
    @test getroot(tree) == "NewRoot"

    # Assign continuous trait
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
    traits = DataFrame([[1.0], [0.5]], [:start, :σ²])
    assigntraits!(tree, traits)
    @test nrow(gettraits(tree)) == 10
    resettraits!(tree)
    @test_throws ErrorException gettraits(tree)

    # Assign discrete trait
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
    switch = [0.5, 0.5]
    traits = DataFrame(trait1 = collect(1:2))
    assigntraits!(tree, switch, traits)
    @test nrow(gettraits(tree)) == 10
    resettraits!(tree)
    @test_throws ErrorException gettraits(tree)
end

# The phylogenetic trait models, which had NO test of any kind until 2026-08-28. That is not a
# hypothetical gap: `varcovar` called `indmax`, removed in Julia 0.7, so it threw an `UndefVarError`
# on **every** call for years, and it was found by moving the file rather than by running it.
#
# `varcovar` is the one worth pinning properly, because the other two are built on it and because it
# has an exact, checkable definition: entry (i, j) is the root-to-tip distance shared by tips i and
# j - the distance from the root to their most recent common ancestor.
@testset "varcovar is the shared-ancestry matrix" begin
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(6))
    V = varcovar(tree)
    tips = collect(nodenamefilter(isleaf, tree))
    root = first(collect(nodenamefilter(isroot, tree)))

    @test size(V) == (length(tips), length(tips))
    @test issymmetric(V)

    # On an ultrametric tree every tip is the same distance from the root, so the diagonal - each
    # tip's shared ancestry with itself, i.e. its whole path - is constant. This is what makes the
    # fixture worth being `Ultrametric` rather than any random tree.
    @test all(≈(distance(tree, root, first(tips))), diag(V))

    # And the defining property, asserted against the tree rather than against stored numbers: an
    # off-diagonal entry is the root-to-MRCA distance for that pair.
    for i in 1:length(tips), j in 1:length(tips)
        i < j || continue
        mrca = argmax(map(x -> distance(tree, root, x),
                          getancestors(tree, tips[i]) ∩
                          getancestors(tree, tips[j])))
        shared = getancestors(tree, tips[i]) ∩ getancestors(tree, tips[j])
        @test V[i, j] ≈ distance(tree, root, shared[mrca])
    end

    # Shared ancestry cannot exceed a tip's own path to the root, which is what an off-diagonal
    # entry larger than the diagonal would mean.
    @test all(V[i, j] <= V[i, i] + 1e-9 for i in axes(V, 1), j in axes(V, 2))
end

# The two fitted models. Their numbers come from an optimiser, so what is pinned is the shape of the
# answer and the relationship between the two - not coefficients that would re-bless on any change
# to `Optim`.
@testset "fitbrownian returns a usable fitted model" begin
    Random.seed!(20260828)
    tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(8))
    tips = collect(nodenamefilter(isleaf, tree))
    # Random rather than a linear ramp: an evenly spaced trait is close to degenerate against an
    # ultrametric tree's covariance.
    traits = randn(length(tips))

    bm = fitbrownian(tree, traits)
    @test bm isa EcoSISTEM.Brownian
    @test isfinite(bm.LL)
    @test length(bm.optimum) == 2               # σ² and the root state
    @test all(isfinite, bm.se)

    # Prints as itself rather than as a parameter dump - it defines its own `show`.
    @test occursin("Brownian", sprint(show, bm))
end

# `fitlambda`, its `Lambda` model and the released `fitLambda` were all DELETED in v0.5.0 - with no
# shim, because the function never worked on Julia 1.x and so no result of it can be depended on.
#
# Two independent faults: `varcovar` threw `UndefVarError` on every call from Julia 0.7 (`indmax`),
# and the fit itself threw `ArgumentError: matrix contains Infs or NaNs` - bounded optimisation on
# lambda meant `Fminbox`, which needs a gradient, and the objective's `log(abs(det(x[1] * V)))`
# underflows for a scaled n x n covariance. `fitbrownian` shares the expression and survives only
# because its optimisation is unbounded, so Nelder-Mead never probes where it blows up.
@testset "the lambda model is gone, in every namespace it lived in" begin
    for n in (:fitlambda, :Lambda, :fitLambda)
        @test !isdefined(EcoSISTEM, n)
        @test !isdefined(EcoSISTEM.ClimatePref, n)
    end
end

end
