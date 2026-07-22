# SPDX-License-Identifier: LGPL-3.0-or-later

module TestTraitfuns

using EcoSISTEM
using Test
using Distributions
using Unitful, Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using AxisArrays

grid = (10, 10)
area = 25.0km^2
totalK = 10000.0kJ / km^2
active = fill(true, grid)
@testset "Trait functions" begin
    # TEST simplehabitatAE
    fillval = 1.0
    habitat1 = simplehabitatAE(fillval, grid, totalK, area)
    habitat2 = tempgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K / month)

    regime = RegimeCollection2(habitat1.regime, habitat2.regime)
    tolerance = ToleranceCollection2(NicheTolerance(Unclassified, Normal,
                                                    fill(1.0, 10),
                                                    fill(0.1, 10)),
                                     NicheTolerance(MeanTemperature, Normal,
                                                    fill(1.0K, 10),
                                                    fill(0.1K, 10)))
    rel = multiplicativeTR2(DistRel{Float64}(), DistRel{Unitful.Temperature}())
    @test_nowarn EcoSISTEM._traitfun(regime, tolerance, rel, 1, 1)
    @test getpref(tolerance, :one) == tolerance.one
    @test getpref(tolerance, :two) == tolerance.two
    @test EcoSISTEM.getrelationship(rel, :one) == rel.one
    @test EcoSISTEM.getrelationship(rel, :two) == rel.two

    temp = AxisArray(fill(1.0K, 10, 10, 3),
                     Axis{:latitude}((1:10) .* °),
                     Axis{:longitude}((1:10) .* °),
                     Axis{:time}(collect(1:3) .* s))
    eratemp = ERA(temp)
    active = fill(true, 10, 10)
    solar = SolarTimeSupply(fill(10.0kJ, 10, 10, 3), 1)
    ea = eraAE(eratemp, solar, active)
    regime = ea.regime
    @test_nowarn EcoSISTEM._traitfun(regime, tolerance.two, rel.two, 1, 1)
    tolerance = TempTolerance(Array(hcat(fill(collect(1:4), 10)...)'))
    rel = DistRel{Unitful.Temperature}()
    @test_nowarn EcoSISTEM._traitfun(regime, tolerance, rel, 1, 1)
    @test getpref(tolerance, 1) == params(getdist(tolerance, 1)) ==
          (1.0, 2.0, 3.0, 4.0)

    rain = AxisArray(fill(1.0mm, 10, 10, 3),
                     Axis{:latitude}((1:10) .* °),
                     Axis{:longitude}((1:10) .* °),
                     Axis{:time}(collect(1:3) .* s))
    erarain = ERA(rain)
    active = fill(true, 10, 10)
    solar = SolarTimeSupply(fill(10.0kJ, 10, 10, 3), 1)
    ea = eraAE(erarain, solar, active)
    regime = ea.regime
    # a Uniform response takes 2 parameters, so each species' row has 2 entries
    tolerance = RainTolerance(Array(hcat(fill(collect(1:2), 10)...)'))
    rel = DistRel{typeof(1.0mm)}()
    @test_nowarn EcoSISTEM._traitfun(regime, tolerance, rel, 1, 1)
    @test getpref(tolerance, 1) == params(getdist(tolerance, 1)) == (1.0, 2.0)
end

@testset "NicheTolerance trait function is distribution-generic + type-stable" begin
    temp = AxisArray(fill(1.0K, 10, 10, 3),
                     Axis{:latitude}((1:10) .* °),
                     Axis{:longitude}((1:10) .* °),
                     Axis{:time}(collect(1:3) .* s))
    solar = SolarTimeSupply(fill(10.0kJ, 10, 10, 3), 1)
    regime = eraAE(ERA(temp), solar, fill(true, 10, 10)).regime

    # A Normal response (neither Trapezoid nor Uniform) proves any continuous distribution works.
    bin = NicheTolerance(MeanTemperature, Normal,
                         Array(hcat(fill([1.0, 2.0], 10)...)'))
    rel = DistRel{typeof(1.0K)}()
    @test EcoSISTEM._traitfun(regime, bin, rel, 1, 1) ==
          pdf(Normal(1.0, 2.0), 1.0)
    # hot-loop type stability (hard constraint): the `D(row...)` splat must not leak Any/Union.
    @test @inferred(EcoSISTEM._traitfun(regime, bin, rel, 1, 1)) isa Float64

    # the relationship's unit `TR` is imputed from the trait's axis / the regime — not typed by hand
    @test typeof(DistRel(bin)) == DistRel{typeof(1.0K)}
    @test eltype(DistRel(bin)) == eltype(bin)
    @test EcoSISTEM._default_relationship(bin, regime) isa DistRel{typeof(1.0K)}
end

end
