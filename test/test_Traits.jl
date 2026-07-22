# SPDX-License-Identifier: LGPL-3.0-or-later

module TestTraits

using EcoSISTEM
using Test
using Distributions
using Unitful, Unitful.DefaultSymbols
using EcoSISTEM.Units
using Phylo
using DataFrames

import EcoSISTEM: DiscreteTolerance

@testset "Traits" begin
    numSpecies = 10
    opts = fill(5.0°C, numSpecies)
    vars = rand(Uniform(0, 25 / 9), numSpecies) * °C

    @testset "gaussian trait (Normal NicheTolerance)" begin
        # A Gaussian preference is a `NicheTolerance` with a `Normal` response on a named axis; the °C input vectors
        # are read per role and built in the default (canonical K) frame.
        bin = NicheTolerance(MeanTemperature, Normal, opts, vars)
        @test bin isa NicheTolerance{MeanTemperature}
        @test EcoSISTEM.iscontinuous(bin) == true
        @test eltype(bin) == typeof(1.0K)
        # Params are stored bare in the canonical (K) frame. `σ` is a standard deviation (a temperature
        # *interval*): 9°F, 5°C and 5K all give the same 5 K width; `μ` keeps the affine offset
        # (5°C → 278.15 K).
        @test params(getdist(NicheTolerance(MeanTemperature, Normal,
                                            [300.0u"°F"],
                                            [9.0u"°F"]),
                             1))[2] ≈ 5.0
        @test params(getdist(NicheTolerance(MeanTemperature, Normal, [5.0°C],
                                            [5.0°C]), 1))[2] ==
              5.0
        @test params(getdist(NicheTolerance(MeanTemperature, Normal, [278.15K],
                                            [5.0K]),
                             1))[2] ==
              5.0
        @test params(getdist(NicheTolerance(MeanTemperature, Normal, [5.0°C],
                                            [5.0°C]), 1))[1] ==
              278.15
    end
    @testset "vector constructor: imputation, bare + mixed units" begin
        # reads the vectors' °C unit, builds in the default (canonical K) frame → identical to the
        # explicit-matrix K form
        @test params(getdist(NicheTolerance(MeanTemperature, Normal, [5.0°C],
                                            [2.0°C]), 1)) ==
              params(getdist(NicheTolerance(MeanTemperature, Normal,
                                            [278.15 2.0];
                                            support = K),
                             1))
        # a bare (dimensionless) shape param alongside a unitful scale (Gamma: α, θ) imputes θ's unit
        @test params(getdist(NicheTolerance(Precipitation, Gamma, [2.0],
                                            [3.0mm]), 1)) ==
              (2.0, 3.0)
        # bare vectors on a dimensioned axis are read in its canonical unit
        @test params(getdist(NicheTolerance(MeanTemperature, Normal, [274.0],
                                            [2.0]), 1)) ==
              (274.0, 2.0)
        # differing units across parameter vectors error
        @test_throws ErrorException NicheTolerance(MeanTemperature, Normal,
                                                   [5.0°C],
                                                   [2.0K])
    end
    # The deprecated `GaussTrait` shim is covered in `test/test_deprecations.jl`.
    @testset "discrete trait" begin
        # Discrete trait
        @test_nowarn DiscreteTolerance(fill(1, 10))
        @test EcoSISTEM.iscontinuous(DiscreteTolerance(fill(1, 10))) == false
        @test eltype(DiscreteTolerance(fill(1, 10))) <: Int
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
        @test eltype(RainTolerance(repeat([1 2], 10))) <: Unitful.Length
        @test_nowarn ToleranceCollection2(TempTolerance(repeat([1 2 3 4], 10)),
                                          RainTolerance(repeat([1 2], 10)))
    end
    @testset "support = distribution frame (matrix + vector)" begin
        # `support` is the frame the distribution is built in (= its eltype = the required regime unit),
        # default the axis's canonical unit; bare matrix numbers are taken as-is *in that frame*.
        @test params(getdist(NicheTolerance(Precipitation, Gamma, [2.0 3.0]),
                             1)) ==
              (2.0, 3.0)                                          # default mm frame
        b = NicheTolerance(MeanTemperature, Normal, [0.0 2.0]; support = u"°C")
        @test collect(params(getdist(b, 1))) ≈ [0.0, 2.0]        # °C frame, bare as-is
        @test eltype(b) == typeof(1.0u"°C")
        @test collect(params(getdist(NicheTolerance(MeanTemperature, Uniform,
                                                    [0.0 10.0];
                                                    support = u"°C"), 1))) ≈
              [0.0, 10.0]
        # the vector constructor reads the inputs' unit and converts to the `support` frame per role:
        # a location properly (0°C → 273.15 K on a K frame), a scale as an interval (2°C → 2 K).
        @test collect(params(getdist(NicheTolerance(MeanTemperature, Normal,
                                                    [0.0°C],
                                                    [2.0°C];
                                                    support = K), 1))) ≈
              [273.15, 2.0]
        @test collect(params(getdist(NicheTolerance(MeanTemperature, Uniform,
                                                    [0.0°C],
                                                    [10.0°C]; support = K), 1))) ≈
              [273.15, 283.15]
        # a support (frame) of the wrong dimension errors
        @test_throws ErrorException NicheTolerance(Precipitation, Gamma,
                                                   [2.0 3.0];
                                                   support = u"°C")
        # a shape-only distribution (Beta) needs offset + scale to be placed on a dimensioned frame
        @test_throws ErrorException NicheTolerance(MeanTemperature, Beta,
                                                   [2.0 3.0])
        bb = NicheTolerance(MeanTemperature, Beta, [2.0 3.0]; support = K,
                            offset = 270.0,
                            scale = 30.0)
        @test getdist(bb, 1) isa Distributions.LocationScale
        # the K and °C frames encode the *same* preference: the density at a physical value is
        # frame-independent when evaluated with `ustrip(support, x)` in each frame (the zero-cost contract,
        # and what the hot path's `ustrip(current)` does when the regime is in the frame's unit).
        bK = NicheTolerance(MeanTemperature, Normal, [5.0°C], [2.0°C];
                            support = K)
        bC = NicheTolerance(MeanTemperature, Normal, [5.0°C], [2.0°C];
                            support = °C)
        @test eltype(bK) == typeof(1.0K)
        @test eltype(bC) == typeof(1.0°C)
        @test pdf(getdist(bK, 1), ustrip(K, 6.0°C)) ≈
              pdf(getdist(bC, 1), ustrip(°C, 6.0°C))
    end
    @testset "multiple tolerance" begin
        # Multiple tolerance
        tr2 = ToleranceCollection2(TempTolerance(repeat([1 2 3 4], 10)),
                                   RainTolerance(repeat([1 2], 10)))
        @test EcoSISTEM.iscontinuous(tr2) == [true, true]
        @test eltype(tr2) == [typeof(1.0K), typeof(1.0mm)]
        gbin = NicheTolerance(MeanTemperature, Normal, opts, vars)
        @test_nowarn ToleranceCollection3(gbin,
                                          TempTolerance(repeat([1 2 3 4], 10)),
                                          RainTolerance(repeat([1 2], 10)))
        tr3 = ToleranceCollection3(gbin,
                                   TempTolerance(repeat([1 2 3 4], 10)),
                                   RainTolerance(repeat([1 2], 10)))
        @test EcoSISTEM.iscontinuous(tr3) == [true, true, true]
        @test eltype(tr3) == [typeof(1.0K), typeof(1.0K), typeof(1.0mm)]
    end
    @testset "evolution" begin
        # Test evolution of tolerance
        tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
        @test_nowarn DiscreteEvolve(2, tree, 0.5)
        @test typeof(get_traits(tree)) <: DataFrame
        @test maximum(get_traits(tree)[!, :trait1]) == 2

        tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
        @test_nowarn ContinuousEvolve(1.0, 0.1, tree)
        @test typeof(get_traits(tree)) <: DataFrame
        @test all(get_traits(tree)[!, :σ²] .== 0.1)
        @test length(unique(get_traits(tree)[!, :start])) > 1
    end
end

end
