# using SimulationData
#
# @testset "simulationdata" begin
#     config = joinpath("simulationdata", "config.yaml")
#     dataconfig = joinpath("simulationdata/demographics", "data_config.yaml")
#     accessfile = joinpath("simulationdata", "access-example.yaml")
#     remove_accessfile() = rm(accessfile, force=true)
#     remove_data() = rm("simulationdata/demographics/human", recursive=true)
#
#     # Basic tests to check integration
#
#     @testset "Basic usage" begin
#         remove_accessfile()
#         @test !isfile(accessfile)
#         api = StandardAPI(config, "test_uri", "test_git_sha")
#         @test read_estimate(api, "parameter", "example-estimate") == 1.0
#         @test !isfile(accessfile)
#         close(api)
#         @test isfile(accessfile)
#         remove_accessfile()
#     end
#
#     @testset "Do-block usage" begin
#         remove_accessfile()
#         @test !isfile(accessfile)
#         StandardAPI(config, "test_uri", "test_git_sha") do api
#             @test read_estimate(api, "parameter", "example-estimate") == 1.0
#         end
#         @test isfile(accessfile)
#         remove_accessfile()
#     end
#     @testset "Population data" begin
#         try
#             download_data_registry(dataconfig)
#             api = StandardAPI(dataconfig, "test_uri", "test_git_sha")
#             scotpop = parse_scottish_population(api)
#         catch e
#             println("Can't download from boydorr.gla.ac.uk ftp server")
#             @test_broken isfile("simulationdata/demographics/human/demographics/population/scotland/1.0.0/1.0.0.h5")
#             @test_broken size(scotpop) == (465, 691, 10)
#             @test_broken sum(scotpop .< 0) == 0
#             @test_broken sum(scotpop) > 5e6
#         finally
#             remove_data()
#         end
#     end
# end
