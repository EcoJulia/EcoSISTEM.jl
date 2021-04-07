using EcoSISTEM
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units

grid = (10, 10)
area = 25.0km^2
totalK = 10000.0kJ/km^2
active = fill(true, grid)
@testset "Trait functions" begin
    # TEST simplehabitatAE
    fillval = 1.0
    abenv1 = simplehabitatAE(fillval, grid, totalK, area)
    abenv2 = tempgradAE(-10.0K, 10.0K, grid, totalK, area, 0.01K/month)

    hab = HabitatCollection2(abenv1.habitat, abenv2.habitat)
    trts = TraitCollection2(GaussTrait(fill(1.0, 10), fill(0.1, 10)), GaussTrait(fill(1.0K, 10), fill(0.1K, 10)))
    rel = multiplicativeTR2(Gauss{Float64}(), Gauss{Unitful.Temperature}())
    @test_nowarn EcoSISTEM._traitfun(hab, trts, rel, 1, 1)
    @test getpref(trts, :t1) == trts.t1
    @test getpref(trts, :t2) == trts.t2
    @test EcoSISTEM.getrelationship(rel, :tr1) == rel.tr1
    @test EcoSISTEM.getrelationship(rel, :tr2) == rel.tr2

    temp = AxisArray(fill(1.0K, 10, 10, 3), Axis{:latitude}(1:10), Axis{:longitude}(1:10), Axis{:time}(collect(1:3) .* s))
    eratemp = ERA(temp)
    active = fill(true, 10, 10)
    solar = SolarTimeBudget(fill(10.0kJ, 10, 10, 3), 1)
    ea = eraAE(eratemp, solar, active)
    hab = ea.habitat
    @test_nowarn EcoSISTEM._traitfun(hab, trts.t2, rel.tr2, 1, 1)
    trts = TempBin(Array(hcat(fill(collect(1:4), 10)...)'))
    rel = Trapeze{Unitful.Temperature}()
    @test_nowarn EcoSISTEM._traitfun(hab, trts, rel, 1, 1)
    @test getpref(trts, 1) == trts.dist[1, :]


    rain = AxisArray(fill(1.0mm, 10, 10, 3), Axis{:latitude}(1:10), Axis{:longitude}(1:10), Axis{:time}(collect(1:3) .* s))
    erarain = ERA(rain)
    active = fill(true, 10, 10)
    solar = SolarTimeBudget(fill(10.0kJ, 10, 10, 3), 1)
    ea = eraAE(erarain, solar, active)
    hab = ea.habitat
    trts = RainBin(Array(hcat(fill(collect(1:4), 10)...)'))
    rel = Unif{typeof(1.0mm)}()
    @test_nowarn EcoSISTEM._traitfun(hab, trts, rel, 1, 1)
    @test getpref(trts, 1) == trts.dist[1, :]

end
