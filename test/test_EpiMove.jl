module TestEpiMove

using EcoSISTEM
using Test
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using DataFrames

@testset "Epi movement" begin
    kernels = dispersal_dists = fill(2.0km, 10 * 10)
    kernel = GaussianKernel.(dispersal_dists, 1e-10)

    @test_nowarn EpiMovement(kernel)
    move = EpiMovement(kernel)
    @test typeof(move.localmoves) <: AlwaysMovement
    @test typeof(move.regionmoves) <: LongDistance
    @test move.localmoves.kernels == kernel

    move_record = DataFrame([[1], [1], [1.0]], [:from, :to, :count])
    @test_nowarn EpiMovement(kernel, move_record)
    move = EpiMovement(kernel, move_record)
    @test typeof(move.localmoves) <: AlwaysMovement
    @test typeof(move.regionmoves) <: LongDistance
    @test move.localmoves.kernels == kernel
    @test move.regionmoves.move_record == move_record
end

end
