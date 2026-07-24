# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDist

using EcoSISTEM
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using Unitful
using Distributions
using Random

import EcoSISTEM: read_distribution, param_roles_resolved, all_positions,
                  role_units, param_roles, guess_params
# `param_units` is public (unexported); reach it through the module.
const param_units = EcoSISTEM.param_units

@testset "Dist" begin
    @testset "Trapezoid distribution" begin
        a, b, c, d = 1, 2, 3, 4
        dist = Trapezoid(a, b, c, d)
        @test params(dist) == (dist.a, dist.b, dist.c, dist.d)
        # `Integer` arguments are widened to `Float64`
        @test eltype(params(dist)) == Float64
        @test_nowarn rand(dist)
        @test_nowarn rand(MersenneTwister(1), dist)
        # piecewise-linear density: rises on [a,b], flat 0.5 on [b,c], falls on [c,d], 0 outside
        @test pdf(dist, 1.0) == 0
        @test pdf(dist, 4.0) == 0
        @test pdf(dist, 2.0) == 0.5
        @test pdf(dist, 3.0) == 0.5
        @test pdf(dist, 1.5) == 0.25          # halfway up the rising edge
        @test pdf(dist, 3.5) == 0.25          # halfway down the falling edge
        @test pdf(dist, 0.5) == 0.0           # below the support
        # support bounds
        @test minimum(dist) == 1.0
        @test maximum(dist) == 4.0

        # the argument guard: `a ≤ b ≤ c ≤ d` and `a < d` (non-degenerate)
        @test_throws DomainError Trapezoid(4, 3, 2, 1)   # reversed
        @test_throws DomainError Trapezoid(1, 1, 1, 1)   # a == d (point mass)
        @test_nowarn Trapezoid(0, 0, 1, 1)               # flat top/bottom is allowed

        # zero-argument default is the uniform distribution on [0, 1]
        u = Trapezoid()
        @test params(u) == (0.0, 0.0, 1.0, 1.0)
        @test pdf(u, 0.25) == pdf(u, 0.5) == pdf(u, 0.9) == 1.0
    end

    @testset "parameter role resolution" begin
        # location/scale, all-position (bounds), shape/scale, and shape-only families
        @test param_roles_resolved(Normal) == [:location, :scale]
        @test param_roles_resolved(Uniform) == [:location, :location]
        @test param_roles_resolved(Trapezoid) ==
              [:location, :location, :location, :location]
        @test param_roles_resolved(Gamma) == [:shape, :scale]
        @test param_roles_resolved(Beta) == [:shape, :shape]
        # LogNormal/LogitNormal are special-cased to shape-only (log/logit-space params)
        @test param_roles_resolved(LogNormal) == [:shape, :shape]
        @test param_roles_resolved(LogitNormal) == [:shape, :shape]

        # `all_positions`: a bounds family's support moves with its parameters; a fixed-support
        # (Beta) or infinite-support (Normal) family's does not
        @test all_positions(Uniform, param_roles(Uniform)) == true
        @test all_positions(Beta, param_roles(Beta)) == false
        @test all_positions(Normal, param_roles(Normal)) == false

        # `guess_params` reconstructs a representative instance of the family
        @test guess_params(Normal) isa Normal
        @test guess_params(Gamma) isa Gamma
        # `param_roles` returns the introspection named tuple
        info = param_roles(Normal)
        @test info.roles == [[:location], [:scale]]
        @test info.dist isa Normal
    end

    @testset "role_units + param_units (absolute-unit introspection)" begin
        # location/scale carry the ABSOLUTE support unit (K for °C), rate its inverse, shape none
        @test role_units(:location, K) == K
        @test role_units(:scale, °C) == K          # absolute unit of an affine one
        @test role_units(:rate, mm) == inv(mm)
        @test role_units(:shape, K) == NoUnits
        @test_throws ErrorException role_units(:bogus, K)

        @test param_units(Normal, K) == [K, K]
        @test param_units(Normal, °C) == [K, K]     # affine → absolute
        @test param_units(Uniform, mm) == [mm, mm]
        @test param_units(Gamma, mm) == [NoUnits, mm]  # shape dimensionless, scale carries it
        # shape-only family: every parameter is dimensionless (regression — Beta must not be
        # mistaken for a bounds family and given [K, K])
        @test param_units(Beta, K) == [NoUnits, NoUnits]
        @test param_units(Normal, NoUnits) == [NoUnits, NoUnits]
    end

    @testset "read_distribution (affine-aware value conversion)" begin
        # no unit ⇒ bare passthrough
        @test params(read_distribution(Normal, NoUnits, [1.0, 2.0])) ==
              (1.0, 2.0)

        # location is a position (proper affine), scale is an interval (width): on the default K
        # frame, 0 °C → 273.15 K but a 2 °C width → 2 K
        @test collect(params(read_distribution(Normal, °C, [0.0, 2.0]))) ≈
              [273.15, 2.0]
        # in the °C frame the location stays 0, the width still 2
        @test collect(params(read_distribution(Normal, °C, [0.0, 2.0];
                                               canonical = °C))) ≈ [0.0, 2.0]
        # both Uniform bounds are positions → both shift
        @test collect(params(read_distribution(Uniform, °C, [0.0, 10.0]))) ≈
              [273.15, 283.15]
        # a K frame leaves K inputs unchanged
        @test collect(params(read_distribution(Normal, K, [274.0, 2.0]))) ≈
              [274.0, 2.0]
        # Gamma: shape bare, scale a width in mm
        @test collect(params(read_distribution(Gamma, mm, [2.0, 3.0]))) ≈
              [2.0, 3.0]

        # a shape-only family is placed on the axis via a LocationScale from offset + scale
        bb = read_distribution(Beta, K, [2.0, 3.0]; offset = 270.0,
                               scale = 30.0)
        @test bb isa Distributions.LocationScale
        @test bb.μ == 270.0
        @test bb.σ == 30.0
        @test params(bb.ρ) == (2.0, 3.0)

        # error paths
        @test_throws ErrorException read_distribution(Normal, K, [1.0])   # wrong param count
        @test_throws ErrorException read_distribution(Beta, K, [2.0, 3.0]) # shape-only, no offset/scale
    end
end

end
