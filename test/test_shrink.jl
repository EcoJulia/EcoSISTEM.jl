using EcoSISTEM
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using AxisArrays

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
expected_grid= (6, 8)
expected_active = active[4:9, 2:9]
# Shouldn't change
expected_gridlength = 1km
control = NoControl()
expected_matrix = fill(fillval, expected_grid)

@testset "_shrink_to_active" begin
    @testset "shrink a normal matrix" begin
        M = rand(grid...)
        M_ref = copy(M)
        M_shrunk = shrink_to_active(M, active)
        @test M_shrunk isa AxisArray
        @test axisvalues(M_shrunk) == (4:9, 2:9)
        @test size(M_shrunk) == expected_grid
        @test M_shrunk == M[4:9, 2:9]
        @test M == M_ref
    end
    @testset "automatically identify inactive regions" begin
        M = Matrix{Union{Float64, Missing}}(rand(grid...))
        M[.!active] .= NaN
        M[1, 1] = missing
        M_ref = copy(M)
        M_shrunk = shrink_to_active(M, active)
        @test M_shrunk isa AxisArray
        @test axisvalues(M_shrunk) == (4:9, 2:9)
        @test size(M_shrunk) == expected_grid
        # Should have preserved the NaNs
        @test any(isnan.(M_shrunk))
        @test isequal(M_shrunk, M[4:9, 2:9])
        @test isequal(M, M_ref)
    end
    @testset "shrink an AxisArray matrix" begin
        M = AxisArray(
            rand(grid...);
            x=Symbol.("x_", 1:grid[1]),
            y=Symbol.("y_", 1:grid[2]),
        )
        M_ref = copy(M)
        M_shrunk = shrink_to_active(M, active)
        @test M_shrunk isa AxisArray
        @test axisvalues(M_shrunk) == (Symbol.("x_", 4:9), Symbol.("y_", 2:9))
        @test size(M_shrunk) == expected_grid
        @test M_shrunk == M[4:9, 2:9]
    end
end

@testset "convert_population" begin
    M = Matrix{Union{Float64, Missing}}(rand(grid...))
    M[.!active] .= NaN
    M[1, 1] = missing
    M_ref = copy(M)

    @testset "Integer type $U" for U in (Int64, UInt64)
        M_expected = U.(round.(replace(M, NaN => 0, missing => 0)))

        M_converted = convert_population(M, U(1))

        @test M_converted isa Matrix{U}
        @test M_converted == M_expected
        # Should not modify M
        @test isequal(M, M_ref)
    end
end
