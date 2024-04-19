module TestDataCleaning

using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using AxisArrays
using Test

@testset "up- and down-resolution testing" begin
    ar2 = AxisArray(collect(reshape(1.0:81.0, 9, 9)))
    @test all(upresolution(downresolution(ar2, 2), 2) .≈ ar2)
    @test all(downresolution(upresolution(ar2, 3), 3) .≈ ar2)

    ar3 = AxisArray(collect(reshape(1.0:25.0, 5, 5, 1)))
    @test all(upresolution(downresolution(ar3, 2), 2) .≈ ar3)
    @test all(downresolution(upresolution(ar3, 3), 3) .≈ ar3)

    ar2b = AxisArray(collect(reshape(1.0:81.0, 9, 9)))
    art = upresolution(ar2b, 2)
    downresolution!(ar2b.data, art.data, 2)
    @test all(ar2b .≈ ar2)
end

end
