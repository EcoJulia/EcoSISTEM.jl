module ExtEcoSISTEMPlotsExt

using Unitful, Unitful.DefaultSymbols
using EcoSISTEM, EcoSISTEM.Units
using AxisArrays
using Plots
using Test

@testset "Plotting EcoSISTEM items" begin
    temp = AxisArray(fill(1.0K, 11, 11, 21),
                     Axis{:latitude}((0°):(1°):(10°)),
                     Axis{:longitude}((0°):(1°):(10°)),
                     Axis{:time}((2000year):(1year):(2020year)))
    era = EcoSISTEM.ERA(temp)
    @test plot(era, 2001year).n == 1
    @test plot(era, 2002year, 1° .. 4°, 5° .. 10°).n == 1
end

end
