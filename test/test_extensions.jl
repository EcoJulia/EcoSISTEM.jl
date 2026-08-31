# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests for `src/extensions.jl` - the hooks the parent declares for an extension to implement.
#
# This file must go through `Pkg.test`, not `julia --project=.. test_extensions.jl`: it exercises
# `compress_landcover`, whose sole method lives in `EcoSISTEMRasterDataSourcesExt`, and that trigger
# package is a weak dependency and so is absent from the package environment.
#
# What is asserted here is the hook's *behaviour* once the extension supplies it. Whether the hook
# gains a method at all is asserted in `ext_EcoSISTEMRasterDataSourcesExt.jl`, which tests the seam.

module TestExtensions

using EcoSISTEM
using Unitful, Unitful.DefaultSymbols
using DimensionalData
using DimensionalData: Dim
using RasterDataSources
using Test

@testset "compress_landcover" begin
    # A tiny synthetic 2×2 grid, 3 layers - the winning layer *position* (by construction, the
    # same thing `compress_landcover` reports as the class code) is known per cell.
    vals = zeros(2, 2, 3)
    vals[1, 1, :] = [1.0, 5.0, 2.0]   # layer 2 wins
    vals[1, 2, :] = [5.0, 1.0, 2.0]   # layer 1 wins
    vals[2, 1, :] = [1.0, 2.0, 5.0]   # layer 3 wins
    vals[2, 2, :] = [1.0, 5.0, 2.0]   # layer 2 wins
    raster = DimArray(vals, (Y(1:2), X(1:2), Dim{:layer}([:a, :b, :c])))
    landcover = ClimateRaster(EarthEnv{LandCover}, raster)
    compressed = compress_landcover(landcover)
    # **Derived, not `EarthEnv{LandCover}`** - an intended change, not a regression. A
    # dominant-class layer is computed *from* land cover and is not land cover: its values are class
    # codes where the inputs were cover fractions, so claiming the dataset would attach that
    # dataset's catalogue row to a quantity it does not describe. The lineage survives in the
    # parameter, which is what `DerivedData` is for.
    @test compressed isa
          ClimateRaster{EcoSISTEM.DerivedData{EarthEnv{LandCover}}}
    # ...and a derived layer carries no layer code, for the same reason.
    @test isnothing(compressed.code)
    @test compressed.array isa DimensionalData.AbstractDimArray
    @test Array(compressed.array) == [2 1; 3 2]
end

end
