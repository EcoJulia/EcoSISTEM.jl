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
    @testset "trait line: GaussTrait → Bin" begin
        opts = fill(5.0K, 4)
        vars = fill(2.0K, 4)
        # axis form → the same `Normal` `Bin`
        @test_deprecated GaussTrait(MeanTemperature, opts, vars)
        gt = GaussTrait(MeanTemperature, opts, vars)
        @test gt isa Bin{MeanTemperature}
        @test params(getdist(gt, 1)) ==
              params(getdist(Bin(MeanTemperature, Normal, opts, vars), 1))

        # axis-less *bare* form → an `Unclassified` Bin (eltype Float64)
        @test_deprecated GaussTrait([1.0, 2.0], [0.1, 0.2])
        gb = GaussTrait([1.0, 2.0], [0.1, 0.2])
        @test eltype(gb) == Float64
        @test params(getdist(gb, 1)) ==
              params(getdist(Bin(Unclassified, Normal, [1.0, 2.0], [0.1, 0.2]),
                             1))

        # axis-less *unitful* form (doubly deprecated): infers the axis from the unit, and warns about it
        @test_deprecated GaussTrait(opts, vars)                 # K → MeanTemperature
        gu = GaussTrait(opts, vars)
        @test gu isa Bin{MeanTemperature}
        @test params(getdist(gu, 1)) ==
              params(getdist(Bin(MeanTemperature, Normal, opts, vars), 1))
        rain = fill(3.0mm, 4)
        gr = GaussTrait(rain, fill(1.0mm, 4))                   # mm → Precipitation
        @test gr isa Bin{Precipitation}
        # a unit with no canonical axis cannot be inferred — a clear error, not a MethodError
        @test_throws ErrorException GaussTrait(fill(1.0u"kg", 2),
                                               fill(1.0u"kg", 2))
    end

    @testset "trait line: Gauss / Trapeze / Unif → DistRel" begin
        TR = typeof(1.0K)
        @test_deprecated Gauss{TR}()
        @test_deprecated Trapeze{Int64}()
        @test_deprecated Unif{typeof(1.0mm)}()
        # the shims share `DistRel`'s 2-argument density functor
        @test Gauss{TR}()(Normal(1.0, 0.01), 1.0K) ==
              DistRel{TR}()(Normal(1.0, 0.01), 1.0K)
        @test Trapeze{Int64}()(Trapezoid(1, 2, 3, 4), 1) ==
              DistRel{Int64}()(Trapezoid(1, 2, 3, 4), 1)
        @test Unif{typeof(1.0mm)}()(Uniform(1, 2), 1.0mm) ==
              DistRel{typeof(1.0mm)}()(Uniform(1, 2), 1.0mm)
        @test eltype(Gauss{TR}()) == TR
        @test EcoSISTEM.iscontinuous(Unif{typeof(1.0mm)}()) == true

        # restored legacy 3-argument Gaussian functor `Gauss{TR}()(current, opt, sd)`
        g = Gauss{TR}()
        @test ustrip(g(275.0K, 274.0K, 2.0K)) ≈ pdf(Normal(274.0, 2.0), 275.0)
    end

    @testset "climate line: per-source constructors → ClimateRaster" begin
        aa = AxisArray(rand(5, 5),
                       Axis{:latitude}((1:5) .* °),
                       Axis{:longitude}((1:5) .* °))
        cases = [(Worldclim_bioclim, WorldClim{BioClim}),
            (CHELSA_bioclim, CHELSA{BioClim}),
            (Landcover, EarthEnv{LandCover}),
            (Worldclim_monthly, WorldClim{Climate}),
            (CHELSA_monthly, CHELSA{Climate})]
        for (shim, src) in cases
            @test_deprecated shim(aa)
            cr = shim(aa)
            @test cr isa ClimateRaster
            @test cr.array === aa
            @test typeof(cr) == typeof(ClimateRaster(src, aa))
        end
    end

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

    @testset "resource line: Budget → Supply" begin
        # the v0.4.0 `*Budget` layer types are deprecated aliases of the renamed `*Supply` types.
        # (`@deprecate_binding` warns via a channel `@test_deprecated` can't capture, so we assert the
        # binding resolves; the warning still fires under `--depwarn=yes`.)
        @test SolarBudget === SolarSupply
        @test WaterBudget === WaterSupply
        @test VolWaterBudget === VolWaterSupply
        @test SimpleBudget === SimpleSupply
        @test SolarTimeBudget === SolarTimeSupply
        @test WaterTimeBudget === WaterTimeSupply
        @test VolWaterTimeBudget === VolWaterTimeSupply
        @test BudgetCollection2 === SupplyCollection2
        @test EcoSISTEM.AbstractBudget === EcoSISTEM.AbstractSupply
    end
end

end
