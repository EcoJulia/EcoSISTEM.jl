# SPDX-License-Identifier: LGPL-3.0-or-later

module TestTolerance

using EcoSISTEM
# `[C7-VIS]` B1/B2/B3: these are `public` rather than exported, so they must be named.
using EcoSISTEM: getdist, gettraits, discrete_evolve, continuous_evolve
using Test
using Distributions
using Unitful, Unitful.DefaultSymbols
using EcoSISTEM.Units
using Phylo
using DataFrames

import EcoSISTEM: SimpleCategoricalTolerance

@testset "Traits" begin
    numSpecies = 10
    opts = fill(5.0°C, numSpecies)
    vars = rand(Uniform(0, 25 / 9), numSpecies) * °C

    @testset "gaussian trait (Normal NicheTolerance)" begin
        # A Gaussian preference is a `NicheTolerance` with a `Normal` response on a named axis; the °C input vectors
        # are read per role and built in the default (canonical K) frame.
        bin = NicheTolerance(Temperature, Normal, opts, vars)
        @test bin isa NicheTolerance{Temperature}
        @test EcoSISTEM.iscontinuous(bin) == true
        @test eltype(bin) == typeof(1.0K)
        # Params are stored bare in the canonical (K) frame. `σ` is a standard deviation (a temperature
        # *interval*): 9°F, 5°C and 5K all give the same 5 K width; `μ` keeps the affine offset
        # (5°C → 278.15 K).
        @test params(getdist(NicheTolerance(Temperature, Normal,
                                            [300.0u"°F"],
                                            [9.0u"°F"]),
                             1))[2] ≈ 5.0
        @test params(getdist(NicheTolerance(Temperature, Normal, [5.0°C],
                                            [5.0°C]), 1))[2] ==
              5.0
        @test params(getdist(NicheTolerance(Temperature, Normal, [278.15K],
                                            [5.0K]),
                             1))[2] ==
              5.0
        @test params(getdist(NicheTolerance(Temperature, Normal, [5.0°C],
                                            [5.0°C]), 1))[1] ==
              278.15
    end
    @testset "vector constructor: imputation, bare + mixed units" begin
        # reads the vectors' °C unit, builds in the default (canonical K) frame → identical to the
        # explicit-matrix K form
        @test params(getdist(NicheTolerance(Temperature, Normal, [5.0°C],
                                            [2.0°C]), 1)) ==
              params(getdist(NicheTolerance(Temperature, Normal,
                                            [278.15 2.0],
                                            support = K),
                             1))
        # a bare (dimensionless) shape param alongside a unitful scale (Gamma: α, θ) imputes θ's unit
        @test params(getdist(NicheTolerance(Precipitation, Gamma, [2.0],
                                            [3.0mm / day]), 1)) ==
              (2.0, 3.0)
        # bare vectors on a dimensioned axis are read in its canonical unit
        @test params(getdist(NicheTolerance(Temperature, Normal, [274.0],
                                            [2.0]), 1)) ==
              (274.0, 2.0)
        # differing units across parameter vectors error
        @test_throws ErrorException NicheTolerance(Temperature, Normal,
                                                   [5.0°C],
                                                   [2.0K])
    end
    # The deprecated `GaussTrait` shim is covered in `test/test_deprecations.jl`.
    @testset "categorical trait" begin
        # Categorical trait
        @test_nowarn SimpleCategoricalTolerance(fill(1, 10),
                                                axis = EcoSISTEM.TypologyAxis)
        @test EcoSISTEM.iscontinuous(SimpleCategoricalTolerance(fill(1, 10),
                                                                axis = EcoSISTEM.TypologyAxis)) ==
              false
        @test eltype(SimpleCategoricalTolerance(fill(1, 10),
                                                axis = EcoSISTEM.TypologyAxis)) <:
              Int
    end
    # The merge of `CategoricalTolerance` and `LandCoverTolerance`: one type taking either
    # spelling, with the soft/hard distinction demoted from a *type* to a **number**.
    @testset "one categorical tolerance takes either spelling" begin
        # The two spellings must give the *same* stored shape, or the merge has not happened —
        # a single preferred class is a set of size one, which is the whole argument for merging.
        single = SimpleCategoricalTolerance([1, 2, 3],
                                            axis = EcoSISTEM.TypologyAxis)
        sets = SimpleCategoricalTolerance([[1], [2], [3]],
                                          axis = EcoSISTEM.TypologyAxis)
        @test single.vals == sets.vals == [[1], [2], [3]]
        @test EcoSISTEM.getpref(single, 2) == [2]
        # A genuine set per species works, and is the same type: one class or several is a value,
        # not a type distinction.
        multi = SimpleCategoricalTolerance([[1, 2], [3]],
                                           axis = EcoSISTEM.TypologyAxis)
        @test typeof(multi) === typeof(single)
        @test EcoSISTEM.getpref(multi, 1) == [1, 2]

        # The penalty is the whole soft/hard axis, and `_categoryweight` is what the hot loop
        # asks. Inside the set it is always 1; outside it is whatever was asked for.
        soft = SimpleCategoricalTolerance([[1, 2], [3]],
                                          axis = EcoSISTEM.TypologyAxis,
                                          penalty = 0.5)
        @test EcoSISTEM._categoryweight(soft, 1, 2) == 1.0
        @test EcoSISTEM._categoryweight(soft, 1, 9) == 0.5
        @test EcoSISTEM._categoryweight(multi, 1, 9) == 0.0    # the default: hard exclusion
        @test EcoSISTEM._categoryweight(SimpleCategoricalTolerance([1],
                                                                   axis = EcoSISTEM.TypologyAxis,
                                                                   penalty = 0.25),
                                        1, 9) == 0.25

        # A weight outside [0, 1] is not a suitability, and is refused at the constructor rather
        # than producing a birth/death factor nothing else would question.
        @test_throws ErrorException SimpleCategoricalTolerance([1],
                                                               axis = EcoSISTEM.TypologyAxis,
                                                               penalty = 1.5)
        @test_throws ErrorException SimpleCategoricalTolerance([1],
                                                               axis = EcoSISTEM.TypologyAxis,
                                                               penalty = -0.1)
        @test_nowarn SimpleCategoricalTolerance([1],
                                                axis = EcoSISTEM.TypologyAxis,
                                                penalty = 1.0)
        @test_nowarn SimpleCategoricalTolerance([1],
                                                axis = EcoSISTEM.TypologyAxis,
                                                penalty = 0.0)

        # The merged type carries a niche axis, which neither of its predecessors did — this is
        # what lets a categorical tolerance be checked against its regime at all.
        @test EcoSISTEM.axisof(SimpleCategoricalTolerance([1],
                                                          axis = EcoSISTEM.NicheAxis)) ===
              EcoSISTEM.NicheAxis
        @test EcoSISTEM.axisof(SimpleCategoricalTolerance([1],
                                                          axis = LandCoverTypology)) ===
              LandCoverTypology
    end
    @testset "temperature bin" begin
        # Temperature bin
        @test_nowarn TempTolerance(repeat([1 2 3 4], 10))
        @test EcoSISTEM.iscontinuous(TempTolerance(repeat([1 2 3 4], 10))) ==
              true
        @test eltype(TempTolerance(repeat([1 2 3 4], 10))) <:
              Unitful.Temperature
    end
    @testset "rainfall bin" begin
        # Rainfall bin
        @test_nowarn RainTolerance(repeat([1 2], 10))
        @test EcoSISTEM.iscontinuous(RainTolerance(repeat([1 2], 10))) == true
        @test eltype(RainTolerance(repeat([1 2], 10))) == typeof(1.0mm / day)
        @test_nowarn SpeciesRequirementCollection((TempTolerance(repeat([1 2 3 4
                                                                         ], 10)),
                                                   RainTolerance(repeat([1 2],
                                                                        10))))
    end
    @testset "support = distribution frame (matrix + vector)" begin
        # `support` is the frame the distribution is built in (= its eltype = the required regime unit),
        # default the axis's canonical unit; bare matrix numbers are taken as-is *in that frame*.
        @test params(getdist(NicheTolerance(Precipitation, Gamma, [2.0 3.0]),
                             1)) ==
              (2.0, 3.0)                                          # default mm frame
        b = NicheTolerance(Temperature, Normal, [0.0 2.0], support = °C)
        @test collect(params(getdist(b, 1))) ≈ [0.0, 2.0]        # °C frame, bare as-is
        @test eltype(b) == typeof(1.0°C)
        @test collect(params(getdist(NicheTolerance(Temperature, Uniform,
                                                    [0.0 10.0],
                                                    support = °C), 1))) ≈
              [0.0, 10.0]
        # the vector constructor reads the inputs' unit and converts to the `support` frame per role:
        # a location properly (0°C → 273.15 K on a K frame), a scale as an interval (2°C → 2 K).
        @test collect(params(getdist(NicheTolerance(Temperature, Normal,
                                                    [0.0°C],
                                                    [2.0°C],
                                                    support = K), 1))) ≈
              [273.15, 2.0]
        @test collect(params(getdist(NicheTolerance(Temperature, Uniform,
                                                    [0.0°C],
                                                    [10.0°C], support = K), 1))) ≈
              [273.15, 283.15]
        # a support (frame) of the wrong dimension errors
        @test_throws ErrorException NicheTolerance(Precipitation, Gamma,
                                                   [2.0 3.0],
                                                   support = °C)
        # a shape-only distribution (Beta) needs offset + scale to be placed on a dimensioned frame
        @test_throws ErrorException NicheTolerance(Temperature, Beta,
                                                   [2.0 3.0])
        bb = NicheTolerance(Temperature, Beta, [2.0 3.0], support = K,
                            offset = 270.0,
                            scale = 30.0)
        @test getdist(bb, 1) isa Distributions.LocationScale
        # the K and °C frames encode the *same* preference: the density at a physical value is
        # frame-independent when evaluated with `ustrip(support, x)` in each frame (the zero-cost contract,
        # and what the hot path's `ustrip(current)` does when the regime is in the frame's unit).
        bK = NicheTolerance(Temperature, Normal, [5.0°C], [2.0°C],
                            support = K)
        bC = NicheTolerance(Temperature, Normal, [5.0°C], [2.0°C],
                            support = °C)
        @test eltype(bK) == typeof(1.0K)
        @test eltype(bC) == typeof(1.0°C)
        @test pdf(getdist(bK, 1), ustrip(K, 6.0°C)) ≈
              pdf(getdist(bC, 1), ustrip(°C, 6.0°C))
    end
    @testset "multiple tolerance" begin
        # Multiple tolerance
        # **Named**: both deprecated builders put their tolerance on the ROOT axis, so derived
        # names cannot tell them apart. That is the rule working — two members claiming nothing are
        # exactly as indistinguishable as two claiming the same thing.
        warmth = TempTolerance(repeat([1 2 3 4], 10))
        wet = RainTolerance(repeat([1 2], 10))
        tr2 = SpeciesRequirementCollection((; warmth, wet))
        @test map(EcoSISTEM.iscontinuous, values(tr2)) == (true, true)
        @test map(eltype, values(tr2)) == (typeof(1.0K), typeof(1.0mm / day))
        gbin = NicheTolerance(Temperature, Normal, opts, vars)
        # **Named** throughout: `gbin` and the deprecated `TempTolerance` are both on a
        # temperature axis, and `warmth`/`wet` are both on the root — neither pair can be told
        # apart by a derived name.
        niche = gbin
        @test_nowarn SpeciesRequirementCollection((; niche, warmth, wet))
        tr3 = SpeciesRequirementCollection((; niche, warmth, wet))
        @test map(EcoSISTEM.iscontinuous, values(tr3)) == (true, true, true)
        @test map(eltype, values(tr3)) ==
              (typeof(1.0K), typeof(1.0K), typeof(1.0mm / day))
    end
    @testset "evolution" begin
        # Test evolution of tolerance
        tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
        @test_nowarn discrete_evolve(2, tree, 0.5)
        @test typeof(gettraits(tree)) <: DataFrame
        @test maximum(gettraits(tree)[!, :trait1]) == 2

        tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
        @test_nowarn continuous_evolve(1.0, 0.1, tree)
        @test typeof(gettraits(tree)) <: DataFrame
        @test all(gettraits(tree)[!, :σ²] .== 0.1)
        @test length(unique(gettraits(tree)[!, :start])) > 1
    end
end

end
