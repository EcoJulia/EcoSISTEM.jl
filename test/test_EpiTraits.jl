using EcoSISTEM
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units
import EcoSISTEM: traitfun

include("TestCases.jl")

@testset "Epi traits" begin
    epi = TestEpiSystem()
    @test_nowarn traitfun(epi, 1, 1, epi.spplist.pathogens)
    @test traitfun(epi, 1, 1, epi.spplist.pathogens) == traitfun(epi, 2, 1, epi.spplist.pathogens)
end
