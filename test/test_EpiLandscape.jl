module TestEpiLandscape

using EcoSISTEM
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units
import EcoSISTEM: emptyepilandscape

include("TestCases.jl")

@testset "EpiLandscape" begin
    epi = Test1EpiSystem()

    @test_nowarn emptyepilandscape(epi.abenv, epi.spplist, Int64(1))
    el = emptyepilandscape(epi.abenv, epi.spplist, Int64(1))
    @test size(el.grid, 2) * size(el.grid, 3) == size(el.matrix, 2)
    @test size(el.grid, 1) == size(el.matrix, 1)
    @test sum(el.grid) == sum(el.matrix)
    el.matrix[1, 1] += 1
    @test sum(el.grid) == sum(el.matrix)
    el.grid[1, 1, 1] += 1
    @test sum(el.grid) == sum(el.matrix)
end

end
