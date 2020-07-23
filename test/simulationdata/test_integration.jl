using SimulationData

@testset "simulationdata" begin
    config = joinpath("simulationdata", "config.yaml")
    accessfile = joinpath("simulationdata", "access-example.yaml")
    remove_accessfile() = rm(accessfile, force=true)

    # Basic tests to check integration

    @testset "Basic usage" begin
        remove_accessfile()
        @test !isfile(accessfile)
        api = StandardAPI(config, "test_uri", "test_git_sha")
        @test read_estimate(api, "parameter", "example-estimate") == 1.0
        @test !isfile(accessfile)
        close(api)
        @test isfile(accessfile)
        remove_accessfile()
    end

    @testset "Do-block usage" begin
        remove_accessfile()
        @test !isfile(accessfile)
        StandardAPI(config, "test_uri", "test_git_sha") do api
            @test read_estimate(api, "parameter", "example-estimate") == 1.0
        end
        @test isfile(accessfile)
        remove_accessfile()
    end
end
