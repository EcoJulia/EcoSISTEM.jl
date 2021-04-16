using EcoSISTEM
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units
import EcoSISTEM: traitfun

include("TestCases.jl")

@testset "Epi traits" begin
    epi = TestEpiSystem()
    @test_nowarn traitfun(epi, 1, 1)
    @test traitfun(epi, 1, 1) == traitfun(epi, 2, 1)
end
