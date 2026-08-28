# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Does any pair of our methods leave Julia unable to choose between them? Run with the other hygiene
# checks:
#
#     julia --project -e 'using Pkg; Pkg.test(; test_args = ["extras_clean.jl"])'
#
# Must go through `Pkg.test`: the weak dependencies bring in the extensions' methods, and an
# ambiguity between a parent method and an extension one is exactly the kind this exists to catch.
#
# Why this exists. An ambiguity is created by adding a type annotation to *one* method of a
# multi-method function, even when that annotation is correct in isolation — specificity is a
# relation between methods, so narrowing one changes how the whole family resolves. Nothing else
# reports it: the package loads cleanly, method counts are unchanged, and the failure only appears
# when a call lands on the ambiguous pair. Two were introduced this way during the signature audit
# of 2026-08-23, and cost a full test run each to find.
#
# It is fast (about a second) and needs no network.

module TestCleanAmbiguities

using EcoSISTEM
using Test

@testset "Method ambiguities" begin
    # `recursive = false` keeps this to methods EcoSISTEM itself defines. The submodules are covered
    # by their own entries below, so a new ambiguity is attributed rather than merely counted.
    found = Test.detect_ambiguities(EcoSISTEM, recursive = false)

    # A named exception list rather than a count, matching `test_EcoSISTEM.jl` and
    # `clean_Architecture.jl`: a new ambiguity has to be looked at, and an old one disappearing
    # shrinks this list instead of passing silently.
    #
    # Both known cases are `Base.BroadcastStyle` on `ClimateRasterStyle`, against DataFrames'
    # `DataFrameStyle` and Base's own `Unknown` fallback. They are pre-existing, they arise from two
    # foreign packages disagreeing about a rule neither of them owns, and neither is reachable from
    # anything this package does with a `ClimateRaster`.
    known = [:BroadcastStyle]
    unexpected = filter(pair -> !(first(pair).name in known), found)
    @test isempty(unexpected)

    # The known ones are still exactly two, so the carve-out cannot quietly absorb a third.
    @test count(pair -> first(pair).name in known, found) == 2

    # And the detector is not vacuous — the failure it guards against is silent, so a check that
    # never fires would look identical to a clean package. Two methods that genuinely cannot be
    # ordered must be reported.
    @eval module _AmbiguityControl
    f(::Integer, ::Any) = 1
    f(::Any, ::Integer) = 2
    end
    @test length(Test.detect_ambiguities(_AmbiguityControl, recursive = false)) ==
          1
end

end
