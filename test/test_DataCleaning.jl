# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDataCleaning

using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using AxisArrays
using RasterDataSources
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

@testset "wrapper up/down-resolution (ClimateRaster, ERA)" begin
    arr3 = AxisArray(collect(reshape(1.0:75.0, 5, 5, 3)),
                     Axis{:latitude}(collect(1:5) .* °),
                     Axis{:longitude}(collect(1:5) .* °),
                     Axis{:time}(collect(1:3) .* month))
    cr3 = ClimateRaster(WorldClim{Climate}, arr3)
    @test downresolution(cr3, 2) isa ClimateRaster            # default fn = mean
    @test downresolution(cr3, 2; fn = maximum) isa ClimateRaster
    @test upresolution(cr3, 2) isa ClimateRaster

    # a 2-D ClimateRaster (bioclim) routes to the keyword 2-D method
    cr2 = ClimateRaster(WorldClim{BioClim},
                        AxisArray(collect(reshape(1.0:81.0, 9, 9)),
                                  Axis{:latitude}(collect(1:9) .* °),
                                  Axis{:longitude}(collect(1:9) .* °)))
    @test downresolution(cr2, 2; fn = maximum) isa ClimateRaster

    @test downresolution(ERA(arr3), 2) isa ERA
end

@testset "up/down-resolution preserves axis names" begin
    # 2-D lat×long: the result must keep (:latitude, :longitude), not swap them.
    aa2 = AxisArray(collect(reshape(1.0:81.0, 9, 9)),
                    Axis{:latitude}(collect(1:9) .* °),
                    Axis{:longitude}(collect(1:9) .* °))
    @test AxisArrays.axisnames(downresolution(aa2, 3)) ==
          (:latitude, :longitude)
    @test AxisArrays.axisnames(upresolution(aa2, 3)) == (:latitude, :longitude)

    # 3-D with a :time third axis.
    aa3 = AxisArray(collect(reshape(1.0:75.0, 5, 5, 3)),
                    Axis{:latitude}(collect(1:5) .* °),
                    Axis{:longitude}(collect(1:5) .* °),
                    Axis{:time}(collect(1:3) .* month))
    @test AxisArrays.axisnames(downresolution(aa3, 2)) ==
          (:latitude, :longitude, :time)
    @test AxisArrays.axisnames(upresolution(aa3, 2)) ==
          (:latitude, :longitude, :time)

    # 3-D with a :var third axis (bioclim): the third axis name is preserved, not forced to :time.
    aav = AxisArray(collect(reshape(1.0:75.0, 5, 5, 3)),
                    Axis{:latitude}(collect(1:5) .* °),
                    Axis{:longitude}(collect(1:5) .* °),
                    Axis{:var}(1:3))
    @test AxisArrays.axisnames(downresolution(aav, 2)) ==
          (:latitude, :longitude, :var)
end

end
