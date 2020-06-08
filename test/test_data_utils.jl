using HTTP

@testset "Data Utils" begin
    # Download and read in population sizes for Scotland
    file, io = mktemp()
    r = HTTP.request("GET", "https://raw.githubusercontent.com/ScottishCovidResponse/temporary_data/master/human/demographics/scotland/data/demographics.h5")
    write(io, r.body)
    close(io)
    scotpop_10k = parse_hdf5(file, component = "grid10km/5year/persons")
    @test size(scotpop_10k) == (47, 70, 19)
    scotpop_1k = parse_hdf5(file, grid="1km",
                            component = "grid1km/10year/persons")
    @test size(scotpop_1k) == (465, 691, 10)
    @test sum(scotpop_1k) == sum(scotpop_10k)
    @test_throws ArgumentError parse_hdf5(file, grid="4km")
end
