using Simulation
using Test
using Unitful.DefaultSymbols
using Simulation.Units
using DataFrames

@testset "Epi movement" begin
    kernels = dispersal_dists = fill(2.0km, 10 * 10)
    kernel = GaussianKernel.(dispersal_dists, 1e-10)

    @test_nowarn EpiMovement(kernel)
    move = EpiMovement(kernel)
    @test typeof(move.home) <: AlwaysMovement
    @test typeof(move.work) <: Commuting
    @test move.home.kernels == kernel

    work = DataFrame([[1], [1], [1.0]], [:from, :to, :count])
    @test_nowarn EpiMovement(kernel, work)
    move = EpiMovement(kernel, work)
    @test typeof(move.home) <: AlwaysMovement
    @test typeof(move.work) <: Commuting
    @test move.home.kernels == kernel
    @test move.work.home_to_work == work
end
