module TestInference

using EcoSISTEM
using Test
using Unitful.DefaultSymbols
using EcoSISTEM.Units

@testset "Inference" begin
    @testset "SIR wrapper" begin
        # Allocating version:
        param = (beta_env = 1.0/day, beta_force = 1.0/day, sigma = 0.02/day, virus_growth = 1e-3/day, virus_decay = 1e-3/day, mean_dispersal_dist = 10.0km)
        runparams = (times = 2years, interval = 1day, timestep = 1day)
        grid_size = (4,4)
        area = 100.0km^2
        @test_nowarn abuns = SIR_wrapper(grid_size, area, param, runparams)

        # Non-allocating version:
        times = runparams.times; interval = runparams.interval
        Ncells = grid_size[1] * grid_size[2]
        numclasses = 4
        abuns = zeros(Int64, numclasses, Ncells, convert(Int64, floor(times / interval)) + 1)
        @test_nowarn SIR_wrapper!(grid_size, area, param, runparams, abuns)

        @test abuns == SIR_wrapper!(grid_size, area, param, runparams, abuns)
    end
    @testset "SEI3HRD wrapper" begin
        # Allocating version:
        param = (beta_env = 1.0/day, beta_force = 1.0/day, virus_growth_symp = 1e-3/day, virus_growth_asymp = 1e-3/day, virus_growth_presymp = 1e-3/day, virus_decay = 1e-3/day, mean_dispersal_dist = 10.0km)
        runparams = (times = 2years, interval = 1day, timestep = 1day)
        grid_size = (4,4)
        area = 100.0km^2
        @test_nowarn abuns = SEI3HRD_wrapper(grid_size, area, param, runparams)

        # Non-allocating version:
        times = runparams.times; interval = runparams.interval
        Ncells = grid_size[1] * grid_size[2]
        numclasses = 8
        abuns = zeros(Int64, numclasses, Ncells, convert(Int64, floor(times / interval)) + 1)
        @test_nowarn SEI3HRD_wrapper!(grid_size, area, param, runparams, abuns)

        @test abuns == SEI3HRD_wrapper!(grid_size, area, param, runparams, abuns)


        # Multiple age categories:
        age_categories = 10
        param = (beta_env = fill(1.0/day, age_categories), beta_force = fill(1.0/day, age_categories), virus_growth_symp = fill(1e-3/day, age_categories), virus_growth_asymp = fill(1e-3/day, age_categories), virus_growth_presymp = fill(1e-3/day, age_categories), virus_decay = 1e-3/day, mean_dispersal_dist = 10.0km)
        runparams = (times = 2years, interval = 1day, timestep = 1day)
        grid_size = (4,4)
        area = 100.0km^2
        @test_nowarn abuns = SEI3HRD_wrapper(grid_size, area, param, runparams, age_categories)
    end
end

end
