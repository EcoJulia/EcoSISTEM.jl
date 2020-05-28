@testset "Data Utils" begin
    path = Simulation.path("test", "examples", "scrc_demographics.h5")
    scotpop_10k = parse_scotpop(path)
    @test size(scotpop_10k) == (10000, 10000, 19)
    scotpop_1k = parse_scotpop(path, grid="1k")
    @test size(scotpop_1k) == (1000, 1000, 19)
    @test sum(scotpop_1k) == sum(scotpop_10k)
    @test_throws ArgumentError parse_scotpop(path, grid="5k")
end
