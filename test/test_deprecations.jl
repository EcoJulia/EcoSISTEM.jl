# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDeprecations

using EcoSISTEM
# `[C7-VIS]` C: these are `public` rather than exported — a spec is what a user writes,
# and these are what it materialises into.
using EcoSISTEM: RateChange, SeriesLayerChange, AbsoluteChange
# `[C7-VIS]` B1/B2/B3: these are `public` rather than exported, so they must be named.
using EcoSISTEM: getdist
using Extents: Extent
using Test
using Distributions
using Unitful, Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using RasterDataSources
using DimensionalData: DimensionalData, DimArray, Y, X, Ti, Dim

include("TestCases.jl")

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

        # axis-less *bare* form → a root-axis (`EcoSISTEM.NicheAxis`) NicheTolerance (eltype Float64)
        @test_deprecated GaussTrait([1.0, 2.0], [0.1, 0.2])
        gb = GaussTrait([1.0, 2.0], [0.1, 0.2])
        @test eltype(gb) == Float64
        @test params(getdist(gb, 1)) ==
              params(getdist(NicheTolerance(EcoSISTEM.NicheAxis, Normal,
                                            [1.0, 2.0],
                                            [0.1, 0.2]),
                             1))

        # axis-less *unitful* form (doubly deprecated): infers the axis from the unit, and warns about it
        @test_deprecated GaussTrait(opts, vars)                 # K → Temperature
        gu = GaussTrait(opts, vars)
        @test gu isa NicheTolerance{Temperature}
        @test params(getdist(gu, 1)) ==
              params(getdist(NicheTolerance(Temperature, Normal, opts,
                                            vars), 1))
        rain = fill(3.0mm / day, 4)
        gr = GaussTrait(rain, fill(1.0mm / day, 4))             # mm/d → Precipitation
        @test gr isa NicheTolerance{Precipitation}
        # a unit with no canonical axis cannot be inferred — a clear error, not a MethodError
        @test_throws ErrorException GaussTrait(fill(1.0u"kg", 2),
                                               fill(1.0u"kg", 2))
    end

    @testset "trait line: Gauss / Trapeze / Unif → NicheSuitability" begin
        NF = typeof(1.0K)
        @test_deprecated Gauss{EcoSISTEM.NicheAxis, NF}()
        @test_deprecated Trapeze{EcoSISTEM.NicheAxis, Int64}()
        @test_deprecated Unif{EcoSISTEM.NicheAxis, typeof(1.0mm)}()
        # the shims share `NicheSuitability`'s 2-argument density functor
        @test Gauss{EcoSISTEM.NicheAxis, NF}()(Normal(1.0, 0.01), 1.0K) ==
              NicheSuitability{EcoSISTEM.NicheAxis, NF}()(Normal(1.0, 0.01),
                                                          1.0K)
        @test Trapeze{EcoSISTEM.NicheAxis, Int64}()(Trapezoid(1, 2, 3, 4), 1) ==
              NicheSuitability{EcoSISTEM.NicheAxis, Int64}()(Trapezoid(1, 2, 3,
                                                                       4), 1)
        @test Unif{EcoSISTEM.NicheAxis, typeof(1.0mm)}()(Uniform(1, 2),
                                                         1.0mm) ==
              NicheSuitability{EcoSISTEM.NicheAxis, typeof(1.0mm)}()(Uniform(1,
                                                                             2),
                                                                     1.0mm)
        @test eltype(Gauss{EcoSISTEM.NicheAxis, NF}()) == NF
        @test EcoSISTEM.iscontinuous(Unif{EcoSISTEM.NicheAxis, typeof(1.0mm)}()) ==
              true

        # restored legacy 3-argument Gaussian functor `Gauss{EcoSISTEM.NicheAxis, NF}()(current, opt, sd)`
        g = Gauss{EcoSISTEM.NicheAxis, NF}()
        @test ustrip(g(275.0K, 274.0K, 2.0K)) ≈ pdf(Normal(274.0, 2.0), 275.0)
    end

    @testset "climate line: per-source constructors → ClimateRaster" begin
        # A `DimArray`, not an `AxisArray`. What these shims deprecate is the per-source *wrapper
        # type*, not the array type it wraps — and typing them on `AxisArray` had made all five dead,
        # since `ClimateRaster` takes an `AbstractDimArray` after the migration: they failed with a
        # bare `MethodError` instead of warning and working, which is worse than having no shim.
        da = DimArray(rand(5, 5),
                      (Y((1:5) .* °), X((1:5) .* °)))
        cases = [(Worldclim_bioclim, WorldClim{BioClim}),
            (CHELSA_bioclim, CHELSA{BioClim}),
            (Landcover, EarthEnv{LandCover}),
            (Worldclim_monthly, WorldClim{Climate}),
            (CHELSA_monthly, CHELSA{Climate})]
        for (shim, src) in cases
            @test_deprecated shim(da)
            cr = shim(da)
            @test cr isa ClimateRaster
            @test cr.array === da
            @test typeof(cr) == typeof(ClimateRaster(src, da))
        end
    end

    # The reader deprecations need downloaded raster data, so guard them the same way `test_datasetread.jl`
    # does (they are unavailable / slow on Windows CI).
    if !Sys.iswindows()
        @testset "climate line: deprecated readers" begin
            # positional extent (in `°`) → the keyword `cut = Extent(...)` form
            bio1 = getraster(WorldClim{BioClim}, :bio1)
            @test_deprecated readfile(bio1, -10°, 10°, -10°, 10°)
            @test isequal(readfile(bio1, -10°, 10°, -10°, 10°),
                          readfile(bio1,
                                   cut = Extent(Y = (-10°, 10°),
                                                X = (-10°, 10°))))

            # readworldclim → the same `ClimateRaster` the `read`/`_readsource` path builds
            wind = getraster(WorldClim{Climate}, :wind, month = 1:12)
            @test_deprecated readworldclim(WorldClim{Climate}, wind)
            rw = readworldclim(WorldClim{Climate}, wind)
            @test rw isa ClimateRaster
            @test size(rw)[3] == 12

            # readCRUTS/readCHELSA_monthly/readERA/readCERA → `read(T, ...)` dispatched on the
            # result type (CRUTS/ERA/CERA) or source type (CHELSA{Climate}). `CRUTS`/`ClimateRaster`
            # define no `isequal`/`==` of their own (defaults to identity), so compare `.array`,
            # which does compare structurally.
            winddir = dirname(first(wind))
            @test_deprecated readCRUTS(winddir, "tavg")
            @test isequal(readCRUTS(winddir, "tavg").array,
                          read(CRUTS, winddir, "tavg").array)
            @test_deprecated readCHELSA_monthly(winddir, "wind")
            @test isequal(readCHELSA_monthly(winddir, "wind").array,
                          read(CHELSA{Climate}, winddir, "wind").array)
        end
    end

    @testset "resource line: Resource → Supply" begin
        # the v0.4.0 `*Resource` layer types are deprecated aliases of the renamed `*Supply` types.
        # (`@deprecate_binding` warns via a channel `@test_deprecated` can't capture, so we assert the
        # binding resolves; the warning still fires under `--depwarn=yes`.)
        @test SolarBudget === Supply{SolarRadiation}
        @test WaterBudget === Supply{Precipitation}
        # `SimpleBudget` is **removed**, not deprecated: there is no free/dimensionless supply for
        # it to name, so there is nothing to re-point it at. Asserted as absent so the removal
        # cannot be silently undone by re-adding a shim to something else.
        @test !isdefined(EcoSISTEM, :SimpleBudget)
        @test !isdefined(EcoSISTEM, :SimpleSupply)
        # A supply that varies in time is not a separate type: it is a supply carrying a
        # `SeriesLayerChange`, so both released names resolve to the one static type.
        @test SolarTimeBudget === Supply{SolarRadiation}
        @test WaterTimeBudget === Supply{Precipitation}
        @test BudgetCollection2 === EcoSISTEM.SupplyCollection2
    end

    @testset "resource line: Requirement → Demand" begin
        # the v0.4.0 `*Requirement` types are deprecated aliases of the renamed `*Demand` types
        @test SolarRequirement === Demand{SolarRadiation}
        @test WaterRequirement === Demand{Precipitation}
        @test ReqCollection2 === EcoSISTEM.DemandCollection2
        # and they still construct the same object
        @test SolarRequirement(fill(2.0kJ / day, 3)) isa Demand{SolarRadiation}
    end

    @testset "condition line: Condition → Regime" begin
        # the v0.4.0 `*Hab`/`HabitatCollection*` condition-layer types → the renamed `*Regime` types
        @test ContinuousHab === ContinuousRegime
        @test ContinuousTimeHab === ContinuousRegime
        @test DiscreteHab === CategoricalRegime
        @test HabitatCollection2 === EcoSISTEM.RegimeCollection2
        @test HabitatCollection3 === EcoSISTEM.RegimeCollection3
    end

    @testset "environment container: AbioticEnv → Condition" begin
        # the v0.4.0 environment container → the renamed `*Condition` (which now means the whole environment;
        # the condition layer that used to be `AbstractHabitat` is now `AbstractRegime`)
        @test GridAbioticEnv === GridHabitat
    end

    @testset "condition line: Trait → Tolerance" begin
        # the v0.4.0 `*Trait`/`TraitCollection`/`TempBin`/`RainBin` types → the renamed `*Tolerance` types
        # **`DiscreteTrait` and `LCtrait` are no longer BINDINGS**: the two types they named merged
        # into `SimpleCategoricalTolerance`, whose `penalty` carries the distinction their types used
        # to. So each is now a function shim that pins its own released value, and what has to be
        # asserted is the **behaviour** — an identity test could not tell the two apart, since both
        # would name the same type.
        soft = @test_deprecated DiscreteTrait([1, 2, 3])
        hard = @test_deprecated LCtrait([[1, 2], [3], [1, 3]])
        @test soft isa SimpleCategoricalTolerance
        @test hard isa SimpleCategoricalTolerance
        # The values that matter: `DiscreteTrait` scored 0.5 outside a species' class and `LCtrait`
        # scored 0.0. A shim that inherited the new `0.0` default would silently turn `DiscreteTrait`'s
        # soft exclusion into a hard one — the species could no longer live outside its class at all.
        @test soft.penalty == 0.5
        @test hard.penalty == 0.0
        # …and `DiscreteTrait` took one class per species, which must become a one-element set.
        @test EcoSISTEM.getpref(soft, 2) == [2]
        @test EcoSISTEM.getpref(hard, 1) == [1, 2]
        @test TraitCollection2 === EcoSISTEM.ToleranceCollection2
        @test TraitCollection3 === EcoSISTEM.ToleranceCollection3
        @test TempBin === TempTolerance
        @test RainBin === RainTolerance
    end

    @testset "condition line: matcher → NicheFit / Suitability" begin
        # the v0.4.0 matcher types → the renamed `*Fit`/`*Suitability` types (`DistRel` was new this PR,
        # renamed to `NicheSuitability` with no shim)
        # Both now name the **same** type, and that is the point: they differed only in the weight
        # given outside a species' classes, which is no longer a property of the fit.
        @test Match === CategoricalSuitability{EcoSISTEM.NicheAxis}
        @test LCmatch === CategoricalSuitability{EcoSISTEM.NicheAxis}
        @test NoRelContinuous === NoFitContinuous
        @test NoRelDiscrete === NoFitCategorical
        @test multiplicativeTR2 === EcoSISTEM.MultiplicativeFit2
        @test multiplicativeTR3 === EcoSISTEM.MultiplicativeFit3
        @test additiveTR2 === EcoSISTEM.AdditiveFit2
        @test additiveTR3 === EcoSISTEM.AdditiveFit3
    end

    # **The round trip, and it is the assertion that was missing.** Comparing the *bindings*
    # above passes whatever the alias expands to — so when the member tuple silently moved into the
    # axis slot (step C of the type-families work), all nine gates stayed green while every
    # arity-numbered alias stopped matching what its own constructor built. A test that passes
    # either way proves nothing: these use the aliases **as types**.
    @testset "the arity-numbered aliases still match what they construct" begin
        tK = NicheSuitability{Temperature, typeof(1.0K)}()
        tP = NicheSuitability{Precipitation, typeof(1.0mm / day)}()
        tS = NicheSuitability{SolarRadiation, typeof(1.0kJ / day)}()
        f2 = EcoSISTEM.MultiplicativeFit2(tK, tP)
        f3 = EcoSISTEM.MultiplicativeFit3(tK, tP, tS)
        @test f2 isa EcoSISTEM.MultiplicativeFit2{typeof(tK), typeof(tP)}
        @test f3 isa EcoSISTEM.MultiplicativeFit3{typeof(tK), typeof(tP),
                                           typeof(tS)}
        # …and they discriminate: wrong order, and wrong arity, must NOT match
        @test !(f2 isa EcoSISTEM.MultiplicativeFit2{typeof(tP), typeof(tK)})
        @test !(f3 isa EcoSISTEM.MultiplicativeFit2)
        @test !(f2 isa EcoSISTEM.MultiplicativeFit3)
        a2 = EcoSISTEM.AdditiveFit2(tK, tP)
        @test a2 isa EcoSISTEM.AdditiveFit2{typeof(tK), typeof(tP)}
        @test !(a2 isa EcoSISTEM.AdditiveFit2{typeof(tP), typeof(tK)})
        # the collection aliases carry the same shape, so they get the same round trip
        temp = EcoSISTEM.simpleregime(298.0K, 1.0km, (5, 7), Temperature)
        rain = EcoSISTEM.simpleregime(50.0mm / day, 1.0km, (5, 7),
                                      Precipitation)
        r2 = EcoSISTEM.RegimeCollection2(temp, rain)
        @test r2 isa EcoSISTEM.RegimeCollection2{typeof(temp), typeof(rain)}
        @test !(r2 isa EcoSISTEM.RegimeCollection2{typeof(rain), typeof(temp)})
    end

    @testset "layer dynamics: LayerUpdate → AbstractLayerChange" begin
        # v0.4.0's `HabitatUpdate`, `habitatupdate!` and `budgetupdate!` were **unexported**, so
        # they were internal and owed no shim; theirs were removed 2026-08-07 along with fourteen
        # others. See the rule at the head of `src/deprecations.jl`.

        # `LayerUpdate` became the `AbstractLayerChange` hierarchy. The name survives as a
        # constructor mapping each old change function onto its successor; the ignored third
        # argument (a dimension to check the rate against) is still accepted.
        @test (@test_deprecated EcoSISTEM.LayerUpdate(EcoSISTEM.NoChange,
                                                      0.0 / month_mean_duration,
                                                      Unitful.Dimensions{()})) ===
              NoLayerChange()
        @test (@test_deprecated EcoSISTEM.LayerUpdate(TempChange,
                                                      1.0K /
                                                      month_mean_duration)) ==
              EcoSISTEM.SteadyLayerChange{typeof(1.0K / month_mean_duration)}(1.0K /
                                                                              month_mean_duration)
        @test (@test_deprecated EcoSISTEM.LayerUpdate(RainfallChange,
                                                      1.0mm / day /
                                                      month_mean_duration)) isa
              EcoSISTEM.SteadyLayerChange
        # `cyclic_change` walked a layer's own stored stack; a stored series now lives in the
        # layer's change and is installed from that stack when the layer is built, so there is
        # nothing for `LayerUpdate` — which never sees a layer — to make of it.
        @test_throws ErrorException EcoSISTEM.LayerUpdate(cyclic_change,
                                                          0.0 /
                                                          month_mean_duration)
        @test (@test_deprecated EcoSISTEM.LayerUpdate(EcoSISTEM.HabitatLoss,
                                                      1.0 / month_mean_duration)) isa
              EcoSISTEM.LegacyLoss
        @test_throws ErrorException EcoSISTEM.LayerUpdate(sin,
                                                          1.0 /
                                                          month_mean_duration)

        # `TempFluct` is the one whose *semantics* changed: it used to feed the layer's own
        # values back into `sin`, a path-dependent walk, and now oscillates as a function of
        # elapsed simulation time. That is not a silent substitution, so it warns loudly in its
        # own right on top of the deprecation.
        fluct = @test_logs (:warn, r"semantics have changed") match_mode=:any begin
            EcoSISTEM.LayerUpdate(TempFluct, 1.0K / month_mean_duration)
        end
        @test fluct isa EcoSISTEM.PatternedLayerChange{RateChange}

        # The change functions themselves stay exported and callable, applying whatever change the
        # layer now holds.
        eco = Test1Ecosystem()
        @test_deprecated TempChange(eco, eco.habitat.regime,
                                    1month_mean_duration)
        @test_deprecated RainfallChange(eco, eco.habitat.regime,
                                        1month_mean_duration)
        @test_deprecated TempFluct(eco, eco.habitat.regime,
                                   1month_mean_duration)
        @test_deprecated EcoSISTEM.NoChange(eco, eco.habitat.regime,
                                            1month_mean_duration)
        @test eraChange === cyclic_change
        @test worldclimChange === cyclic_change

        # `resetrate!` → `setchange!`. The old name could only install a constant rate on the whole
        # regime; the new one takes any change and addresses one layer, so it also reaches a
        # sub-layer of a collection, which `resetrate!` never could.
        # (`Test1Ecosystem`'s regime is categorical and dimensionless, so its rate is `𝐓⁻¹`.)
        rate = 0.5 / year
        expected = EcoSISTEM._attachchange(IncrementBy(rate),
                                           eco.habitat.regime)
        eco.habitat.regime.change = NoLayerChange()
        @test_deprecated resetrate!(eco, rate)
        @test eco.habitat.regime.change == expected
    end

    @testset "rainfall-gradient builders → GridHabitat" begin
        # `raingrad`/`raingradhabitat`/`raingradAE` are deprecated as a unit and their bodies live
        # in `src/deprecations.jl` (deferred item 11). Behaviour is preserved rather than forwarded
        # to `GridHabitat`, because `rate` has no equivalent there yet — so this testset is
        # the old `test_GridHabitat.jl` "rainfall gradient" coverage, moved and warning-checked.
        grid = (5, 5)
        area = 25.0km^2
        totalK = 10000.0kJ / km^2 / day
        active = fill(true, grid)

        habitat = @test_deprecated raingradhabitat(0.0mm / day, 100.0mm / day,
                                                   grid, totalK, area,
                                                   0.01mm / day /
                                                   month_mean_duration)
        @test habitat isa GridHabitat
        @test minimum(habitat.regime.matrix) == 0.0mm / day
        @test maximum(habitat.regime.matrix) == 100.0mm / day
        @test size(habitat.regime.matrix) == grid
        @test sum(habitat.supply.matrix) == totalK * area
        @test habitat.active == active

        # the underlying generator is deprecated too, and warns in its own right
        @test (@test_deprecated raingrad(0.0mm / day, 100.0mm / day, 1.0km,
                                         grid,
                                         0.01mm / day / month_mean_duration)) isa
              ContinuousRegime

        # `Precipitation`'s canonical unit is the rate `mm/d`, so an accumulated depth is rejected
        @test_throws Unitful.DimensionError raingradhabitat(0.0mm, 100.0mm,
                                                            grid, totalK, area,
                                                            0.01mm /
                                                            month_mean_duration)

        # `raingradAE` resolves onto the same body directly, so it warns once rather than chaining
        # through `raingradhabitat` and warning twice.
        aeh = @test_deprecated raingradAE(0.0mm / day, 100.0mm / day, grid,
                                          totalK, area,
                                          0.01mm / day / month_mean_duration)
        @test aeh isa GridHabitat
        @test aeh.regime.matrix == habitat.regime.matrix

        # **The `maxsupply`-less form, which nothing reached until now** — the gradient's own
        # rainfall *is* the supply. It was broken: it handed the regime's `mm/day` (an areal rate,
        # `Precipitation`'s condition unit) straight to a constructor wanting per-cell `L/day`, so
        # it had no method at all. v0.4.0 could, because both sides were bare `mm`; the v0.5.0 unit
        # split moved them apart and this line was not moved with them.
        own = @test_deprecated raingradhabitat(0.0mm / day, 100.0mm / day, grid,
                                               area,
                                               0.01mm / day /
                                               month_mean_duration)
        @test own isa GridHabitat
        @test own.supply isa Supply{Precipitation}
        # The rainfall is converted against this grid's own cell area (25 km² over 25 cells, so
        # 1 km² each): 100 mm/day over 1 km² is 10⁸ L/day, and the gradient's foot is dry.
        @test maximum(own.supply.matrix) ≈ uconvert(Unitful.L / day,
                       100.0mm / day * 1.0km^2)
        @test minimum(own.supply.matrix) == 0.0Unitful.L / day
        # The regime is untouched by that conversion — it is still the areal rate it was.
        @test own.regime.matrix == habitat.regime.matrix
    end

    @testset "temperature-gradient builders → GridHabitat" begin
        # `tempgrad`/`tempgradhabitat`/`peakedgradhabitat` (+ the `*AE` names) get the same
        # treatment as the rainfall family above: bodies moved into `src/deprecations.jl`, one
        # warning each, behaviour preserved because `rate` has no `GridHabitat` equivalent.
        # This is the old `test_GridHabitat.jl` coverage, moved and warning-checked.
        grid = (5, 5)
        area = 25.0km^2
        totalK = 10000.0kJ / km^2 / day
        active = fill(true, grid)

        for (f, ae) in ((tempgradhabitat, tempgradAE),
            (peakedgradhabitat, peakedgradAE))
            habitat = @test_deprecated f(-10.0K, 10.0K, grid, totalK, area,
                                         0.01K / month_mean_duration)
            @test habitat isa GridHabitat
            @test minimum(habitat.regime.matrix) == -10.0K
            @test maximum(habitat.regime.matrix) == 10.0K
            @test size(habitat.regime.matrix) == grid
            @test sum(habitat.supply.matrix) == totalK * area
            @test habitat.active == active

            # the `*AE` name resolves onto the same body, so it warns once rather than chaining
            aeh = @test_deprecated ae(-10.0K, 10.0K, grid, totalK, area,
                                      0.01K / month_mean_duration)
            @test aeh.regime.matrix == habitat.regime.matrix
        end

        # peaked really does peak in the middle rather than rising monotonically
        peaked = @test_deprecated peakedgradhabitat(-10.0K, 10.0K, grid, totalK,
                                                    area,
                                                    0.01K / month_mean_duration)
        mid = peaked.regime.matrix[ceil(Int, grid[1] / 2), 1]
        @test mid == maximum(peaked.regime.matrix)

        # the underlying generator is deprecated too
        @test (@test_deprecated tempgrad(-10.0K, 10.0K, 1.0km, grid,
                                         0.01K / month_mean_duration)) isa
              ContinuousRegime
    end

    @testset "environment constructors: *AE → *habitat" begin
        # each `*AE` constructor is a deprecated forwarder to its `*habitat` rename; the symbol-form
        # `@deprecate` warns (captured by `@test_deprecated`) and forwards to the same method, so the
        # returned habitat matches. All fixtures are in-memory (no downloads).
        # **Only the five SYNTHETIC families are here now.** The four data-backed ones
        # (`eraAE`/`worldclimAE`/`bioclimAE`/`lcAE`) are removed, and are covered by the testset
        # below, which asserts their message.
        grid = (10, 10)
        totalK = 1000.0kJ / m^2 / day
        area = 100.0km^2

        # synthetic constructors
        h = @test_deprecated simplehabitatAE(5.0K, grid, totalK, area)
        @test h isa GridHabitat
        @test typeof(h) == typeof(simplehabitat(5.0K, grid, totalK, area))
        @test (@test_deprecated tempgradAE(-10.0K, 10.0K, grid, totalK, area,
                                           0.01K / month_mean_duration)) isa
              GridHabitat
        @test (@test_deprecated peakedgradAE(-10.0K, 10.0K, grid, totalK, area,
                                             0.01K / month_mean_duration)) isa
              GridHabitat
        @test (@test_deprecated raingradAE(0.0mm / day, 100.0mm / day, grid,
                                           totalK, area,
                                           0.01mm / day / month_mean_duration)) isa
              GridHabitat
        @test (@test_deprecated simplenicheAE(4, grid, totalK, area)) isa
              GridHabitat
    end

    # **The four data-backed builders are REMOVED, and these assert the message rather than the
    # type.** They built a regime straight from a raster's own cells with no resampling, which no
    # construction route does now — a `StudyArea` decides the grid and the data is sampled onto it —
    # so there is nothing to redirect to and no shim was written.
    # **Asserting on the message, not just `ErrorException`**, because the message *is* the
    # feature: a released exported name that vanished would give `UndefVarError` and tell the reader
    # nothing. A type-only assertion would let the text rot silently.
    @testset "removed environment builders explain themselves" begin
        for f in (eraAE, worldclimAE, bioclimAE, lcAE, erahabitat,
            worldclimhabitat, bioclimhabitat, landcoverhabitat)
            err = try
                f(1, 2, 3)
                nothing
            catch e
                e
            end
            @test err isa ErrorException
            @test occursin("has been removed", err.msg)
            # Both replacements must be named: whoever hits this is *holding a raster*, so
            # `SourceSpec` alone would not answer them.
            @test occursin("in_memory_raster", err.msg)
            @test occursin("StudyArea", err.msg)
            # …and it must say why, not merely that.
            @test occursin("sampling **no** grid", err.msg)
        end
        # The released `*AE` spelling names its own deprecation release, not v0.5.0.
        @test occursin("v0.4.0",
                       (try
                            eraAE(1, 2, 3)
                        catch e
                            e
                        end).msg)
        @test occursin("v0.5.0",
                       (try
                            erahabitat(1, 2, 3)
                        catch e
                            e
                        end).msg)
    end

    @testset "land cover: compressLC → compress_landcover" begin
        latkm = Y(collect(1:10) .* km)
        longkm = X(collect(1:10) .* km)
        lcr = ClimateRaster(EarthEnv{LandCover},
                            DimArray(rand(10, 10, 5),
                                     (latkm, longkm,
                                      DimensionalData.Dim{:layer}(collect(1:5)))))
        compressed = @test_deprecated compressLC(lcr)
        @test compressed.array == compress_landcover(lcr).array
    end
end

@testset "layer time series: a stack and a cursor become a SeriesLayerChange" begin
    # A layer used to hold every slice and an index into them. It now holds one slice and carries
    # the stack as a `SeriesLayerChange` indexed by elapsed time, so the `(stack, time)` constructors
    # build that pair — honouring `time` as the series' origin rather than dropping it.
    stack = cat((fill(i * 10.0kJ / day, 4, 4) for i in 1:3)..., dims = 3)
    supply = @test_deprecated SolarTimeBudget(stack, 1)
    @test supply isa Supply{SolarRadiation}
    @test ndims(supply.matrix) == 2
    @test all(==(10.0kJ / day), supply.matrix)
    @test supply.change isa SeriesLayerChange{AbsoluteChange}
    # elapsed time zero selects the slice the cursor pointed at…
    @test EcoSISTEM._seriesindex(supply.change, 0.0month_mean_duration) == 1
    @test EcoSISTEM._seriesindex(supply.change, 1.0month_mean_duration) == 2

    # …and a cursor that started part-way through anchors the series there instead.
    later = @test_deprecated SolarTimeBudget(stack, 2)
    @test all(==(20.0kJ / day), later.matrix)
    @test EcoSISTEM._seriesindex(later.change, 0.0month_mean_duration) == 2
    @test EcoSISTEM._seriesindex(later.change, 1.0month_mean_duration) == 3

    water = @test_deprecated WaterTimeBudget(cat((fill(i * 1.0Unitful.L / day,
                                                       4, 4) for i in 1:3)...,
                                                 dims = 3), 1)
    @test water isa Supply{Precipitation}
    @test water.change isa SeriesLayerChange

    regime = @test_deprecated ContinuousTimeHab(cat((fill(i * 1.0K, 4, 4)
                                                     for i in 1:3)...,
                                                    dims = 3), 1, 1.0km,
                                                NoLayerChange())
    @test regime isa ContinuousRegime
    @test ndims(regime.matrix) == 2
    @test regime.change isa SeriesLayerChange{AbsoluteChange}
end

# `up`/`downresolution` are deprecated and live in `src/ClimatePref/deprecations.jl`, so their
# tests belong here rather than with the live raster code.

@testset "up- and down-resolution testing" begin
    ar2 = DimArray(collect(reshape(1.0:81.0, 9, 9)), (Y(1:9), X(1:9)))
    @test all(upresolution(downresolution(ar2, 2), 2) .≈ ar2)
    @test all(downresolution(upresolution(ar2, 3), 3) .≈ ar2)

    ar3 = DimArray(collect(reshape(1.0:25.0, 5, 5, 1)),
                   (Y(1:5), X(1:5), Ti(1:1)))
    @test all(upresolution(downresolution(ar3, 2), 2) .≈ ar3)
    @test all(downresolution(upresolution(ar3, 3), 3) .≈ ar3)

    ar2b = DimArray(collect(reshape(1.0:81.0, 9, 9)), (Y(1:9), X(1:9)))
    art = upresolution(ar2b, 2)
    downresolution!(parent(ar2b), parent(art), 2)
    @test all(ar2b .≈ ar2)
end

@testset "wrapper up/down-resolution (ClimateRaster, ERA)" begin
    arr3 = DimArray(collect(reshape(1.0:75.0, 5, 5, 3)),
                    (Y(collect(1:5) .* °), X(collect(1:5) .* °),
                     Ti(collect(1:3) .* month_mean_duration)))
    cr3 = ClimateRaster(WorldClim{Climate}, arr3)
    @test downresolution(cr3, 2) isa ClimateRaster            # default fn = mean
    @test downresolution(cr3, 2, fn = maximum) isa ClimateRaster
    @test upresolution(cr3, 2) isa ClimateRaster

    # a 2-D ClimateRaster (bioclim) routes to the keyword 2-D method
    cr2 = ClimateRaster(WorldClim{BioClim},
                        DimArray(collect(reshape(1.0:81.0, 9, 9)),
                                 (Y(collect(1:9) .* °),
                                  X(collect(1:9) .* °))))
    @test downresolution(cr2, 2, fn = maximum) isa ClimateRaster

    # The released `ERA(array)` call still works, with a warning, and now yields a raster: the type
    # is a source rather than a container, so `isa ERA` cannot be what it once was.
    @test (@test_deprecated ERA(arr3)) isa ClimateRaster{EcoSISTEM.ERA}
    @test downresolution(ClimateRaster(EcoSISTEM.ERA, arr3), 2) isa
          ClimateRaster{EcoSISTEM.ERA}
end

@testset "up/down-resolution preserves dimension names" begin
    # 2-D Y×X: the result must keep (:Y, :X) in that order, not swap them. The dims are
    # `rebuild`t from the input rather than reconstructed, so a third axis of any name survives too.
    dnames(a) = DimensionalData.name.(DimensionalData.dims(a))
    aa2 = DimArray(collect(reshape(1.0:81.0, 9, 9)),
                   (Y(collect(1:9) .* °), X(collect(1:9) .* °)))
    @test dnames(downresolution(aa2, 3)) == (:Y, :X)
    @test dnames(upresolution(aa2, 3)) == (:Y, :X)

    # 3-D with a `Ti` third axis.
    aa3 = DimArray(collect(reshape(1.0:75.0, 5, 5, 3)),
                   (Y(collect(1:5) .* °), X(collect(1:5) .* °),
                    Ti(collect(1:3) .* month_mean_duration)))
    @test dnames(downresolution(aa3, 2)) == (:Y, :X, :Ti)
    @test dnames(upresolution(aa3, 2)) == (:Y, :X, :Ti)

    # 3-D with a `Dim{:var}` third axis (bioclim): its name is preserved, not forced to `Ti`.
    aav = DimArray(collect(reshape(1.0:75.0, 5, 5, 3)),
                   (Y(collect(1:5) .* °), X(collect(1:5) .* °),
                    Dim{:var}(1:3)))
    @test dnames(downresolution(aav, 2)) == (:Y, :X, :var)
end

# Moved from `test_extractclimate.jl` (v0.5.0): it tests a deprecated name, so it belongs with the
# other deprecations -- and this is now the only test file that reaches into `ClimatePref` at all.
# The shim errors on any arguments, so no raster fixture is needed.
@testset "extractvalues is deprecated, and says how to translate" begin
    old_err = try
        extractvalues(25.0°, 15.0°, nothing)
    catch e
        e
    end
    @test old_err isa ErrorException
    @test occursin("extract_values", old_err.msg)
    @test occursin("year", old_err.msg)      # names what replaced each old argument
end

end
