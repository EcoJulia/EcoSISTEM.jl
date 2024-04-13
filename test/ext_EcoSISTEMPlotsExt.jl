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
    @test plot(era, 2000year).n == 1
    @test plot(era, 2002year, 1° .. 4°, 5° .. 10°).n == 1
    cera = EcoSISTEM.CERA(temp)
    @test plot(cera, 2001year).n == 1
    @test plot(cera, 2020year, 1° .. 4°, 5° .. 10°).n == 1

    temp = AxisArray(fill(2.0K, 11, 11, 12),
                     Axis{:latitude}((0°):(1°):(10°)),
                     Axis{:longitude}((0°):(1°):(10°)),
                     Axis{:time}(January:(1month):December))
    wm = EcoSISTEM.Worldclim_monthly(temp)
    @test plot(wm, January).n == 1
    @test plot(wm, September, 0° .. 2°, 5° .. 10°).n == 1
    chelsa = EcoSISTEM.CHELSA_monthly(temp)
    @test plot(chelsa, February).n == 1
    @test plot(chelsa, December, 3° .. 10°, 0° .. 3°).n == 1
end

end
