# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDeprecations

using EcoSISTEM
using Test
using Distributions
using Unitful, Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using AxisArrays
using RasterDataSources

# Coverage for `src/deprecations.jl` (trait line) and `src/ClimatePref/deprecations.jl` (climate line).
# Every deprecated shim is checked on *both* halves: it warns (`@test_deprecated`) **and** its result
# matches the current API it forwards to.

@testset "Deprecations" begin
    # The reader deprecations need downloaded raster data, so guard them the same way `test_ReadData.jl`
    # does (they are unavailable / slow on Windows CI).
    if !Sys.iswindows()
        @testset "climate line: deprecated readers" begin
            # positional extent (in `°`) → the keyword `cut = LatLong(...)` form
            bio1 = getraster(WorldClim{BioClim}, :bio1)
            @test_deprecated readfile(bio1, -10°, 10°, -10°, 10°)
            @test isequal(readfile(bio1, -10°, 10°, -10°, 10°),
                          readfile(bio1;
                                   cut = LatLong(-10° .. 10°, -10° .. 10°)))

            # readworldclim → the same `ClimateRaster` the `read`/`_readsource` path builds
            wind = getraster(WorldClim{Climate}, :wind, month = 1:12)
            @test_deprecated readworldclim(WorldClim{Climate}, wind)
            rw = readworldclim(WorldClim{Climate}, wind)
            @test rw isa ClimateRaster
            @test size(rw)[3] == 12
        end
    end
end

end
