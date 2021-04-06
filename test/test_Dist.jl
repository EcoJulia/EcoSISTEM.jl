using EcoSISTEM
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using Unitful
using Distributions
using Random

@testset "Trapeze distribution" begin
    a = 1; b = 2; c = 3; d = 4
    dist = Trapezoid(a, b, c, d)
    @test params(dist) == (dist.a, dist.b, dist.c, dist.d)
    @test_nowarn rand(dist)
    @test_nowarn rand(MersenneTwister(1), dist)
    @test pdf(dist, 1.0) == 0
    @test pdf(dist, 4.0) == 0
    @test pdf(dist, 2.0) == 0.5
    @test pdf(dist, 3.0) == 0.5
end
