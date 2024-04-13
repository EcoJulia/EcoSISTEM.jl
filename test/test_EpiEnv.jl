module TestEpiEnv

using EcoSISTEM
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using AxisArrays

@testset "Epi environments" begin
    grid = (5, 5)
    area = 25.0km^2
    active = fill(true, grid)
    control = NoControl()
    @testset "simple habitat environments" begin
        # TEST simplehabitatAE
        fillval = 1.0
        abenv = simplehabitatAE(fillval, grid, area, control)
        @test_nowarn simplehabitatAE(fillval, grid, area, control)
        @test_nowarn simplehabitatAE(fillval, grid, area, active, control)
        @test all(abenv.habitat.matrix .== fillval)
        @test size(abenv.habitat.matrix) == grid
        @test abenv.active == active
        @test all(abenv.active)
        @test abenv.habitat.size^2 * grid[1] * grid[2] == area
    end

    @testset "Shrink grid" begin
        grid = (10, 10)
        area = 100.0km^2
        fillval = 1.0

        active = fill(true, grid)
        # Set inactive cells on the outer layers
        # Make this asymmetric to be more general
        # These should be trimmed off
        for row in (1, 2, 3, 10)
            active[row, :] .= false
        end
        for col in (1, 2, 10)
            active[:, col] .= false
        end
        # Set some more inactive cells, these shouldn't result in any further trimming
        active[4, 4] = false
        active[5, 4] = false
        active[5, 5] = false
        active[8, 8] = false
        # This should stop col 2 being trimmed
        active[5, 2] = true
        # Should be trimmed from 10x10 to 6x8
        expected_grid = (6, 8)
        expected_active = active[4:9, 2:9]
        # Shouldn't change
        expected_gridlength = 1km
        control = NoControl()
        expected_matrix = fill(fillval, expected_grid)

        @testset "Construct directly" begin
            epienv = simplehabitatAE(fillval, grid, area, active, control)

            @test epienv.active == expected_active
            @test size(epienv.habitat.matrix) == expected_grid
            @test epienv.habitat.size == expected_gridlength
            @test epienv.habitat.matrix == expected_matrix
        end
    end
end

end
