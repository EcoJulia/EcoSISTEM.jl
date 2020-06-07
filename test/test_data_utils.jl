@testset "Data Utils" begin
    path = Simulation.path("test", "examples", "scrc_demographics.h5")
    scotpop_10k = parse_hdf5(path)
    @test size(scotpop_10k) == (70, 47, 19)
    scotpop_1k = parse_hdf5(path, grid="1k")
    @test size(scotpop_1k) == (691, 465, 19)
    @test sum(scotpop_1k) == sum(scotpop_10k)
    @test_throws ArgumentError parse_hdf5(path, grid="5k")
end
