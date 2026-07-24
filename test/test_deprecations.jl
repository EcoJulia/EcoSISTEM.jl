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
    @testset "trait line: GaussTrait → NicheTolerance" begin
        opts = fill(5.0K, 4)
        vars = fill(2.0K, 4)
        # axis form → the same `Normal` `NicheTolerance`
        @test_deprecated GaussTrait(Temperature, opts, vars)
        gt = GaussTrait(Temperature, opts, vars)
        @test gt isa NicheTolerance{Temperature}
        @test params(getdist(gt, 1)) ==
              params(getdist(NicheTolerance(Temperature, Normal, opts,
                                            vars), 1))

        # axis-less *bare* form → an `Unclassified` NicheTolerance (eltype Float64)
        @test_deprecated GaussTrait([1.0, 2.0], [0.1, 0.2])
        gb = GaussTrait([1.0, 2.0], [0.1, 0.2])
        @test eltype(gb) == Float64
        @test params(getdist(gb, 1)) ==
              params(getdist(NicheTolerance(Unclassified, Normal, [1.0, 2.0],
                                            [0.1, 0.2]),
                             1))

        # axis-less *unitful* form (doubly deprecated): infers the axis from the unit, and warns about it
        @test_deprecated GaussTrait(opts, vars)                 # K → Temperature
        gu = GaussTrait(opts, vars)
        @test gu isa NicheTolerance{Temperature}
        @test params(getdist(gu, 1)) ==
              params(getdist(NicheTolerance(Temperature, Normal, opts,
                                            vars), 1))
        rain = fill(3.0mm, 4)
        gr = GaussTrait(rain, fill(1.0mm, 4))                   # mm → Precipitation
        @test gr isa NicheTolerance{Precipitation}
        # a unit with no canonical axis cannot be inferred — a clear error, not a MethodError
        @test_throws ErrorException GaussTrait(fill(1.0u"kg", 2),
                                               fill(1.0u"kg", 2))
    end

    @testset "trait line: Gauss / Trapeze / Unif → NicheSuitability" begin
        TR = typeof(1.0K)
        @test_deprecated Gauss{TR}()
        @test_deprecated Trapeze{Int64}()
        @test_deprecated Unif{typeof(1.0mm)}()
        # the shims share `NicheSuitability`'s 2-argument density functor
        @test Gauss{TR}()(Normal(1.0, 0.01), 1.0K) ==
              NicheSuitability{TR}()(Normal(1.0, 0.01), 1.0K)
        @test Trapeze{Int64}()(Trapezoid(1, 2, 3, 4), 1) ==
              NicheSuitability{Int64}()(Trapezoid(1, 2, 3, 4), 1)
        @test Unif{typeof(1.0mm)}()(Uniform(1, 2), 1.0mm) ==
              NicheSuitability{typeof(1.0mm)}()(Uniform(1, 2), 1.0mm)
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

    @testset "resource line: Resource → Supply" begin
        # the v0.4.0 `*Resource` layer types are deprecated aliases of the renamed `*Supply` types.
        # (`@deprecate_binding` warns via a channel `@test_deprecated` can't capture, so we assert the
        # binding resolves; the warning still fires under `--depwarn=yes`.)
        @test SolarBudget === SolarSupply
        @test WaterBudget === WaterSupply
        @test SimpleBudget === SimpleSupply
        @test SolarTimeBudget === SolarTimeSupply
        @test WaterTimeBudget === WaterTimeSupply
        @test BudgetCollection2 === SupplyCollection2
        @test EcoSISTEM.AbstractBudget === EcoSISTEM.AbstractSupply
    end

    @testset "resource line: Requirement → Demand" begin
        # the v0.4.0 `*Requirement` types are deprecated aliases of the renamed `*Demand` types
        @test SimpleRequirement === SimpleDemand
        @test SizeRequirement === SizeDemand
        @test SolarRequirement === SolarDemand
        @test WaterRequirement === WaterDemand
        @test ReqCollection2 === DemandCollection2
        @test EcoSISTEM.AbstractRequirement === EcoSISTEM.AbstractDemand
        # and they still construct the same object
        @test SolarRequirement(fill(2.0kJ / day, 3)) isa SolarDemand
    end

    @testset "condition line: Condition → Regime" begin
        # the v0.4.0 `*Hab`/`HabitatCollection*` condition-layer types → the renamed `*Regime` types
        @test ContinuousHab === ContinuousRegime
        @test ContinuousTimeHab === ContinuousTimeRegime
        @test DiscreteHab === DiscreteRegime
        @test HabitatCollection2 === RegimeCollection2
        @test HabitatCollection3 === RegimeCollection3
    end

    @testset "environment container: AbioticEnv → Condition" begin
        # the v0.4.0 environment container → the renamed `*Condition` (which now means the whole environment;
        # the condition layer that used to be `AbstractHabitat` is now `AbstractRegime`)
        @test GridAbioticEnv === GridHabitat
        @test EcoSISTEM.AbstractAbiotic === EcoSISTEM.AbstractHabitat
    end

    @testset "condition line: Trait → Tolerance" begin
        # the v0.4.0 `*Trait`/`TraitCollection`/`TempBin`/`RainBin` types → the renamed `*Tolerance` types
        @test DiscreteTrait === DiscreteTolerance
        @test LCtrait === LandCoverTolerance
        @test TraitCollection2 === ToleranceCollection2
        @test TraitCollection3 === ToleranceCollection3
        @test TempBin === TempTolerance
        @test RainBin === RainTolerance
        @test EcoSISTEM.AbstractTraits === EcoSISTEM.AbstractTolerance
        @test EcoSISTEM.ContinuousTrait === EcoSISTEM.ContinuousTolerance
    end

    @testset "condition line: matcher → NicheFit / Suitability" begin
        # the v0.4.0 matcher types → the renamed `*Fit`/`*Suitability` types (`DistRel` was new this PR,
        # renamed to `NicheSuitability` with no shim)
        @test Match === MatchSuitability
        @test LCmatch === LandCoverSuitability
        @test NoRelContinuous === NoFitContinuous
        @test NoRelDiscrete === NoFitDiscrete
        @test multiplicativeTR2 === multiplicativeFit2
        @test multiplicativeTR3 === multiplicativeFit3
        @test additiveTR2 === additiveFit2
        @test additiveTR3 === additiveFit3
        @test EcoSISTEM.AbstractTraitRelationship === EcoSISTEM.AbstractNicheFit
    end

    @testset "layer dynamics: HabitatUpdate → LayerUpdate" begin
        # the v0.4.0 (unexported) per-layer update rule → the renamed `LayerUpdate`
        @test EcoSISTEM.HabitatUpdate === EcoSISTEM.LayerUpdate
    end

    @testset "environment constructors: *AE → *habitat" begin
        # each `*AE` constructor is a deprecated forwarder to its `*habitat` rename; the symbol-form
        # `@deprecate` warns (captured by `@test_deprecated`) and forwards to the same method, so the
        # returned habitat matches. All fixtures are in-memory (no downloads).
        grid = (10, 10)
        totalK = 1000.0kJ / m^2 / day
        area = 100.0km^2
        latkm = Axis{:latitude}(collect(1:10) .* km)
        longkm = Axis{:longitude}(collect(1:10) .* km)

        # synthetic constructors
        h = @test_deprecated simplehabitatAE(5.0K, grid, totalK, area)
        @test h isa GridHabitat
        @test typeof(h) == typeof(simplehabitat(5.0K, grid, totalK, area))
        @test (@test_deprecated tempgradAE(-10.0K, 10.0K, grid, totalK, area,
                                           0.01K / month)) isa GridHabitat
        @test (@test_deprecated peakedgradAE(-10.0K, 10.0K, grid, totalK, area,
                                             0.01K / month)) isa GridHabitat
        @test (@test_deprecated raingradAE(0.0mm, 100.0mm, grid, totalK, area,
                                           0.01mm / month)) isa GridHabitat
        @test (@test_deprecated simplenicheAE(4, grid, totalK, area)) isa
              GridHabitat

        # data constructors (in-memory ERA / ClimateRaster fixtures)
        eratemp = ERA(AxisArray(fill(1.0K, 10, 10, 3),
                                Axis{:latitude}((1:10) .* °),
                                Axis{:longitude}((1:10) .* °),
                                Axis{:time}(collect(1:3) .* s)))
        @test (@test_deprecated eraAE(eratemp, totalK, area)) isa GridHabitat
        worldclim = ClimateRaster(WorldClim{Climate},
                                  AxisArray(fill(1.0K, 10, 10, 12), latkm,
                                            longkm,
                                            Axis{:time}(collect(1:12) .* month)))
        @test (@test_deprecated worldclimAE(worldclim, totalK, area)) isa
              GridHabitat
        bio = ClimateRaster(WorldClim{BioClim},
                            AxisArray(fill(1.0K, 10, 10), latkm, longkm))
        @test (@test_deprecated bioclimAE(bio, totalK, area)) isa GridHabitat
        lcr = ClimateRaster(EarthEnv{LandCover},
                            AxisArray(fill(1, 10, 10), latkm, longkm))
        lch = @test_deprecated lcAE(lcr, totalK, area)
        @test lch isa GridHabitat
        @test typeof(lch) == typeof(landcoverhabitat(lcr, totalK, area))
    end

    @testset "land cover: compressLC → compressLandCover" begin
        latkm = Axis{:latitude}(collect(1:10) .* km)
        longkm = Axis{:longitude}(collect(1:10) .* km)
        lcr = ClimateRaster(EarthEnv{LandCover},
                            AxisArray(rand(10, 10, 5), latkm, longkm,
                                      Axis{:time}(collect(1:5))))
        compressed = @test_deprecated compressLC(lcr)
        @test compressed.array == compressLandCover(lcr).array
    end
end

end
