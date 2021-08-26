using HTTP
using BritishNationalGrid
import EcoSISTEM: get_en, get_bng, create_BNG_grid
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

@testset "BNG grid" begin
    bng_point = "NS2468"
    en = (224000, 668000)
    @test get_en(bng_point) == en
    @test get_en([bng_point]) == ([en[1]], [en[2]])
    @test get_bng(en[1], en[2]) == "NS 24 68"
    easts = norths = [1000, 2000]
    vals = fill(1.0, 4)
    grid = create_BNG_grid(easts, norths, vals)
    @test size(grid) == (2,2)
    ages = [1, 2]
    grid = create_BNG_grid(easts, norths, vals, ages)
    @test size(grid) == (2,2,2)
end