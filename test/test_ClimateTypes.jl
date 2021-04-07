using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

@testset "ECMWF types" begin
    temp = AxisArray(fill(1.0K, 10, 10, 3), Axis{:latitude}(1:10), Axis{:longitude}(1:10), Axis{:time}(collect(1:3) .* s))
    @test_nowarn ERA(temp)
    @test_nowarn CERA(temp)
    @test_nowarn CRUTS(temp)
end

@testset "Worldclim types" begin
    temp = AxisArray(fill(1.0K, 10, 10, 12), Axis{:latitude}(1:10), Axis{:longitude}(1:10), Axis{:time}(collect(1:12) .* s))
    @test_nowarn Worldclim(temp)
    @test_nowarn CHELSA(temp)
end

@testset "Bioclim types" begin
    temp = AxisArray(fill(1.0K, 10, 10, 19), Axis{:latitude}(1:10), Axis{:longitude}(1:10), Axis{:vars}(collect(1:19)))
    @test_nowarn Bioclim(temp)
end

@testset "Reference types" begin
    ref = AxisArray(fill(1, 10, 10), Axis{:latitude}(1:10), Axis{:longitude}(1:10))
    @test_nowarn Reference(ref)
end
