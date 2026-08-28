# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDiversityInterface

using EcoSISTEM
using EcoSISTEM: makeunique
using Diversity
# `faith_pd`/`generalisedfaith_pd` live in `Diversity.Ecology` and are NOT re-exported into
# `Diversity` itself, so plain `using Diversity` does not reach them.
using Diversity.Ecology
using Phylo
using Statistics
using Test

include("TestCases.jl")

@testset "diversity" begin
    eco = Test1Ecosystem()

    Ab = mapslices(sum, eco.abundances.matrix, dims = 2)
    RelAb = Ab / sum(eco.abundances.matrix)
    ## TEST makeunique
    @test_nowarn makeunique(eco)
    @test typeof(makeunique(eco).spplist.types) <: UniqueTypes

    ## TEST faith_pd (against R and julia)
    # **`faith_pd` is DIVERSITY's now, not ours**, supplied for an ecosystem by
    # `EcoSISTEMPhyloExt` — so this is `Diversity.Ecology.faith_pd`, reached through
    # `using Diversity.Ecology` above, and it carries Diversity's meaning: **per subcommunity**,
    # not one number for the
    # metacommunity, which is what a hand-rolled copy of it would be tempted to return.
    tree = eco.spplist.types.tree
    total_length = sum(map(b -> getlength(tree, b),
                           collect(getbranchnames(tree))))
    # Faith's PD is the **total branch length** of the tree spanned by the types present. Every
    # subcommunity of this fixture holds every species, so each one spans the whole tree — which is
    # why they all equal the total, and is exactly the case that made the metacommunity and
    # subcommunity readings look interchangeable.
    @test length(faith_pd(eco)[!, :diversity]) == countsubcommunities(eco)
    @test all(≈(total_length), faith_pd(eco)[!, :diversity])
    # …and the metacommunity reading is still reachable, now by asking for it explicitly.
    @test generalisedfaith_pd(metacommunityDiversity, eco)[!, :diversity][1] ≈
          total_length
    # No division by `mean(heightstoroot(tree))` is applied, and this is what says it would be a
    # no-op: every tree EcoSISTEM builds is ultrametric with height 1, so the divisor is exactly 1.
    # It is asserted rather than assumed, because leaving it out is only safe while that holds — a
    # tree of any other height would make such a `pd` return something
    # that was not Faith's PD at all.
    @test mean(heightstoroot(tree)) ≈ 1.0
    # There is no `q` to pass. Faith's PD is q = 0 by definition, and the old signature's `qs` was
    # accepted and silently ignored — so `pd(eco, 1.0)` looked like a q = 1 measure and was not.
    @test_throws MethodError faith_pd(eco, 1.0)
end

@testset "makeunique keeps the species names" begin
    # **The names must be non-numeric or this proves nothing.** The bug was `UniqueTypes(count)`,
    # which generates `"1"`, `"2"`, … — indistinguishable from the real names on any list that already
    # has numeric ones, which is every list `build_species` produces *and* `Test1Ecosystem`'s.
    # Only `names` is renamed: `types` is a **type parameter** of `SpeciesList`, so it cannot be
    # reassigned in place — which is exactly why `makeunique` rebuilds the whole species list.
    eco = Test1Ecosystem()
    named = ["sp_" * n for n in eco.spplist.names]
    eco.spplist.names .= named

    unique = makeunique(eco)
    @test typeof(unique.spplist.types) <: UniqueTypes
    @test gettypenames(unique.spplist.types, true) == named
    @test unique.spplist.names == named
end

end
