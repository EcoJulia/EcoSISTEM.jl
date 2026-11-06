# SPDX-License-Identifier: LGPL-3.0-or-later

module TestActions

using EcoSISTEM
using Test
using Distributions
import Dates
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Diversity

include("TestCases.jl")
using EcoSISTEM: hasdata, landcoverclass
using EcoSISTEM: materialise
using Unitful, Unitful.DefaultSymbols
using DimensionalData: DimensionalData, DimArray, X, Y, dims
using Rasters
using RasterDataSources
using ArchGDAL
using Distributions: Normal
using Extents: Extent
using DataFrames: DataFrame
include("rasterfixtures.jl")
include("buildfixtures.jl")

# `TOL` and `DEM` as specs, for the species whose niche and demand no test here is about.
const TOLSPEC = ToleranceSpec(Temperature, μ = TOL[1], σ = TOL[2])
const DEMSPEC = DemandSpec(SolarRadiation, DEM)

# Species built quietly from a pool, or from a pool of `n`, taking `TOLSPEC` and `DEMSPEC` unless
# given others.
function _species(pool::SpeciesPool; kw...)
    return build_species(pool; tolerance = TOLSPEC, demand = DEMSPEC,
                         verbosity = :silent, kw...)
end

function _species(n::Integer; kw...)
    return _species(SpeciesPool(n, verbosity = :silent); kw...)
end

@testset "Simulate functions" begin
    eco = Test1Ecosystem()

    times = 3.0month_mean_duration
    burnin = 1.0month_mean_duration
    interval = 1.0month_mean_duration
    timestep = 1month_mean_duration
    # Run simulation grid
    lensim = length((0.0month_mean_duration):interval:times)
    @testset "simulate" begin
        # Run simulations 10 times
        @test_nowarn generate_storage(eco, lensim, 1)
        abun = generate_storage(eco, lensim, 1)
        @test size(abun) ==
              (size(eco.abundances.matrix, 1), size(eco.abundances.matrix, 2),
               lensim, 1)
        @test_nowarn simulate!(eco, burnin, timestep)
        @test_nowarn simulate!(RecordAbundance(abun), eco, times, timestep,
                               every = EveryInterval(interval))
        # Write the run's output into a temp dir the OS cleans up, so it never
        # touches the repo (and no manual `rm` is needed for the hygiene tests).
        outputdir = mktempdir()
        @test_nowarn simulate!(SaveAbundance(outputdir, "testrun"), eco, times,
                               timestep, every = EveryInterval(interval))
    end
    @testset "simulate! with a callback" begin
        # The starting state first, then the state at exactly each multiple of the interval, over
        # the steps `simulate!` itself takes - three here, so four occurrences at every step.
        step = 1.0month_mean_duration
        eco2 = Test1Ecosystem()
        start = sum(eco2.abundances.matrix)
        seen = []
        @test isnothing(simulate!(eco2, 3step, step) do occurrence
                            return push!(seen,
                                         (occurrence,
                                          sum(eco2.abundances.matrix)))
                        end)
        @test [o.count for (o, _) in seen] == 1:4
        @test [o.elapsed for (o, _) in seen] ≈
              [uconvert(s, k * step) for k in 0:3]
        @test last(first(seen)) == start
        @test all(isnothing(o.date) for (o, _) in seen)

        # A continued run does not observe its starting state again: that was the last occurrence of
        # the run before, and the count starts afresh.
        more = []
        simulate!(eco2, 2step, step) do occurrence
            return push!(more, occurrence)
        end
        @test [o.count for o in more] == 1:2
        @test [o.elapsed for o in more] ≈ [uconvert(s, k * step) for k in 4:5]

        # `every` takes a schedule, and a bare duration means `EveryInterval` of it; a time between
        # steps is acted on at the step that reaches it.
        elapsedof(every) = begin
            e = Test1Ecosystem()
            firings = typeof(1.0s)[]
            simulate!(e, 3step, step, every = every) do occurrence
                return push!(firings, occurrence.elapsed)
            end
            firings
        end
        @test elapsedof(EveryInterval(2step)) ≈
              [uconvert(s, k * step) for k in (0, 2)]
        @test elapsedof(2step) == elapsedof(EveryInterval(2step))
        @test elapsedof(AtTimes([1.5step])) ≈ [uconvert(s, 2step)]

        # What is due at the start acts before the first occurrence sees the state.
        eco3 = Test1Ecosystem()
        seeded = sum(eco3.abundances.matrix) +
                 5 * size(eco3.abundances.matrix, 2)
        atstart = Ref(0)
        simulate!(eco3, 0step, step,
                  intervention = Intervention(AtTime(0.0s), AllCells(),
                                              AddAbundance(1, 5))) do occurrence
            occurrence.count == 1 && (atstart[] = sum(eco3.abundances.matrix))
            return nothing
        end
        @test atstart[] == seeded

        # With an epoch each occurrence carries its date; without one a date schedule is refused.
        eco4 = Test1Ecosystem()
        eco4.epoch = Dates.DateTime(2000, 1, 1)
        dates = []
        simulate!(eco4, step, step) do occurrence
            return push!(dates, occurrence.date)
        end
        @test first(dates) == Dates.DateTime(2000, 1, 1)
        @test_throws "epoch" simulate!(_ -> nothing, Test1Ecosystem(), step,
                                       step,
                                       every = EveryYear())
    end
    @testset "diversity simulate" begin
        # Run diversity simulations 10 times
        qs = collect(1.0:3)
        @test_nowarn generate_storage(eco, length(qs), lensim, 1)
        abun = generate_storage(eco, length(qs), lensim, 1)
        @test size(abun) ==
              (size(eco.abundances.matrix, 2), length(qs), lensim, 1)
        @test_nowarn simulate!(eco, burnin, timestep)
        @test_nowarn simulate!(RecordDiversity(abun, norm_sub_alpha, qs), eco,
                               times, timestep, every = EveryInterval(interval))
    end
end

@testset "build_species" begin
    @testset "scalar fill and defaults" begin
        spp = _species(10, seed = 1)
        @test spp isa SpeciesList
        @test length(spp.names) == 10
        @test spp.demand isa Demand{SolarRadiation}       # kJ/day demand
        @test length(spp.tolerance.dists) == 10
    end

    @testset "per-species values from the pool's columns" begin
        means = [290.0K, 295.0K, 300.0K]
        pool = SpeciesPool(DataFrame(species = ["1", "2", "3"],
                                     mean = ustrip.(K, means)),
                           verbosity = :silent)
        spp = _species(pool,
                       tolerance = ToleranceSpec(Temperature,
                                                 μ = TraitColumn(:mean,
                                                                 unit = K),
                                                 σ = 1.0K),
                       demand = DemandSpec(SolarRadiation, 8.0kJ / day))
        # the built NicheTolerance stores bare means in the canonical (K) frame
        @test EcoSISTEM._nichemeans(spp.tolerance) == ustrip.(K, means)
    end

    # **A multi-layer regime can be inspected, mixing real and synthetic layers.**
    # Inspection must not lag the builder - that is the wrong way round for a function whose whole
    # job is to show what the builder will get. Without this, a tuple falls through to the
    # data-backed path and tries to coerce its first member into a `SourceSpec`, while
    # `GridHabitat` builds exactly such a regime quite happily.
    @testset "materialise takes a multi-layer, mixed regime" begin
        data = _bngraster(WorldClim{BioClim}, fill(285.0K, 9, 9))
        real_ = _reg(data, axis = Temperature)
        syn = UniformSpec(50.0mm / day, axis = Precipitation)
        area = _area(regime = real_)

        # Each kind alone, on the same grid - data sampled onto it, synthetic generated at its shape.
        @test size(EcoSISTEM.materialise(real_, area).matrix) == (9, 9)
        @test size(EcoSISTEM.materialise(syn, area).matrix) == (9, 9)

        # ...and the two together, which is what `GridHabitat` accepts. Several specs give a
        # `LayerCollection`, so the members are reached through the collection interface - which
        # `regimes` gives identically for one layer or many.
        both = EcoSISTEM.materialise((real_, syn), area)
        @test length(values(both)) == 2
        @test all(l -> size(l.matrix) == (9, 9), values(both))
        # A named tuple keeps its names, so members stay identifiable rather than positional.
        named = EcoSISTEM.materialise((temp = real_, rain = syn), area)
        @test keys(NamedTuple(named)) == (:temp, :rain)
        # ...and a role is applied to every member, not just the first.
        @test length(values(EcoSISTEM.materialise((real_, syn), area,
                                                  role = EcoSISTEM.Condition))) ==
              2
    end

    # **A pre-built tolerance is accepted, and that is what makes an axis-less united layer
    # matchable at all.**
    @testset "a pre-built tolerance is used as given" begin
        tol = NicheTolerance(Temperature, Normal, [298.0 2.0; 299.0 2.0])
        sp = _species(2, tolerance = tol, seed = 1)
        @test sp.tolerance === tol
        # A built tolerance carries its own axis.
        @test EcoSISTEM.axisof(sp.tolerance) === Temperature

        # It is used **as given**, so it must already cover every species - a silent mismatch
        # would pair species with the wrong niches.
        @test_throws "must hold one entry per species" _species(5,
                                                                tolerance = tol,
                                                                seed = 1)

        # A pre-built tolerance is used as given, in a frame a `ToleranceSpec` cannot express - here
        # a real unit on `EcoSISTEM.NicheAxis`, where the axis declares none.
        united = NicheTolerance(EcoSISTEM.NicheAxis, Normal,
                                [2.0 0.5; 2.1 0.5], support = u"kJ^2")
        @test _species(2, tolerance = united, seed = 1).tolerance === united
    end

    # **A demand's type comes from its axis, never from its unit.** There is no unit-keyed table to
    # pick it, and there must not be: two quantities can share a unit. The unit is *checked* against
    # the axis, which is a different job from deciding it.
    @testset "demand type from the named axis" begin
        @test _species(4,
                       demand = DemandSpec(Precipitation,
                                           2.0Unitful.L / day)).demand isa
              Demand{Precipitation}
        @test _species(4, demand = DemandSpec(CarbonFlux, 2.0g / day)).demand isa
              Demand{CarbonFlux}

        # **The axis decides, and the unit must agree with it.** A solar-dimensioned demand
        # declared on the water axis is refused, rather than built as a solar demand.
        @test_throws "needs a value in L d⁻¹" _species(4,
                                                       demand = DemandSpec(Precipitation,
                                                                           2.0kJ /
                                                                           day))

        # A **bare number cannot be a demand**: the axis says what unit it expects, and a bare
        # number carries none to compare.
        @test_throws "needs a value in kJ d⁻¹" _species(4,
                                                        demand = DemandSpec(SolarRadiation,
                                                                            3.0))

        # ...and the axis cannot be left out. That is the point: there must be no route by which a
        # demand's meaning can be guessed.
        @test_throws MethodError DemandSpec(2.0kJ / day)
    end

    @testset "names" begin
        anon = _species(3)
        @test anon.names == ["1", "2", "3"]
        binomials = ["Quercus robur", "Fagus sylvatica", "Betula pendula"]
        named = _species(SpeciesPool(binomials, verbosity = :silent))
        @test named.names == binomials
        @test gettypenames(named.types, true) == binomials
    end

    @testset "tolerance and demand are required" begin
        @test_throws ErrorException build_species(SpeciesPool(2,
                                                              verbosity = :silent),
                                                  demand = DEMSPEC,
                                                  verbosity = :silent)
        @test_throws ErrorException build_species(SpeciesPool(2,
                                                              verbosity = :silent),
                                                  tolerance = TOLSPEC,
                                                  verbosity = :silent)
    end

    @testset "one rate fills, a column is per species" begin
        # SpeciesList always stores per-species PopGrowth internally
        s_scalar = _species(2, birth = 0.6 / year)
        @test s_scalar.params.birth == fill(0.6 / year, 2)
        pool = SpeciesPool(DataFrame(species = ["1", "2"], birth = [0.6, 0.5]),
                           verbosity = :silent)
        s_vec = _species(pool, birth = TraitColumn(:birth, unit = year^-1),
                         death = 0.6 / year)
        @test s_vec.params.birth == [0.6, 0.5] ./ year
    end

    @testset "movement mode selection" begin
        species(; kw...) = _species(2; kw...)
        safely(sl) = [EcoSISTEM.dispersesafely(sl.movement, sp) for sp in 1:2]
        @test all(d -> d.mode ==
                       InsectDispersal(),
                  species().movement)
        @test all(d -> d.mode isa AlwaysDispersal,
                  species(movement = MovementSpec(AlwaysDispersal(),
                                                  meandist = 10.0km)).movement)
        # The topology is a property of the **habitat**, not of a species, so what `build_species`
        # decides here is only each species' mode and its `disperse_safely` flag, which the mode
        # carries. Whether the world wraps belongs to the map; whether a disperser survives aiming
        # off it belongs to the disperser.
        fates = SpeciesPool(DataFrame(species = ["1", "2"],
                                      fate = ["safe", "lost"]),
                            verbosity = :silent)
        chosen = ChosenBy(:fate,
                          "safe" =>
                              MovementSpec(AlwaysDispersal(disperse_safely = true),
                                           meandist = 10.0km),
                          "lost" =>
                              MovementSpec(AlwaysDispersal(disperse_safely = false),
                                           meandist = 10.0km))
        @test safely(_species(fates, movement = chosen)) == [true, false]
        # Left unset, the flag is the mode's own: a wind-blown seed aimed at a dead cell is lost.
        @test safely(species(movement = MovementSpec(WindDispersal(GaussianKernel),
                                                     meandist = 10.0km))) ==
              [false, false]
        # Nothing disperses, so no distance is read and nothing can be lost.
        @test species(movement = MovementSpec(NoDispersal())).movement ==
              fill(NoDispersal()(), 2)
        # A family needing more than a mean distance is refused, naming what it lacks.
        @test_throws "needs `scale` in its `MovementSpec`" species(movement = MovementSpec(WindDispersal(),
                                                                                           meandist = 5.0km))
        # Each species' own distance is a column of the pool.
        distances = SpeciesPool(DataFrame(species = ["1", "2"],
                                          meandist = [1.0, 3.0]),
                                verbosity = :silent)
        given = [WindDispersal()(meandist = 1.0km, scale = 2.0km),
            WindDispersal()(meandist = 3.0km, scale = 2.0km)]
        @test _species(distances,
                       movement = MovementSpec(WindDispersal(),
                                               meandist = TraitColumn(:meandist,
                                                                      unit = km),
                                               scale = 2.0km)).movement == given
        @test_throws "`movement` is a `MovementSpec`" species(movement = :birth)
    end
end
@testset "build_ecosystem" begin
    env = GridHabitat(regime = UniformSpec(298.0K, axis = Temperature),
                      supply = SUP,
                      area = _area(extent = (5km, 5km),
                                   cellsize = 1km))
    spp = _species(10, seed = 1)

    @testset "assembly and nichefit inference" begin
        eco = build_ecosystem(BySupply(total = 5000), spp, env, seed = 1)
        @test eco isa Ecosystem
        @test eco.nichefit isa NicheSuitability
        @test sum(eco.abundances.matrix) == 5000
    end

    @testset "an explicit abundance per species" begin
        counts = [10, 20, 30]
        eco = build_ecosystem(BySupply(each = row -> counts[row.key]),
                              _species(3), env, seed = 1)
        @test vec(sum(eco.abundances.matrix, dims = 2)) == counts
    end

    @testset "demand/supply mismatch errors" begin
        water_spp = _species(10,
                             demand = DemandSpec(Precipitation,
                                                 2.0Unitful.L / day),
                             seed = 1)
        @test_throws ErrorException build_ecosystem(BySupply(total = 5000),
                                                    water_spp, env)
    end

    @testset "a built model simulates" begin
        env2 = GridHabitat(regime = UniformSpec(298.0K,
                                                axis = Temperature),
                           supply = UniformSpec(10.0kJ / (km^2 * day),
                                                axis = SolarRadiation),
                           area = _area(extent = (10km, 10km),
                                        cellsize = 1km))
        spp2 = _species(5, demand = DemandSpec(SolarRadiation, 1.0kJ / day),
                        seed = 2)
        eco = build_ecosystem(BySupply(total = 500), spp2, env2, seed = 2)
        @test_nowarn simulate!(eco, 2month_mean_duration, 1month_mean_duration)
    end

    @testset "distributed kwarg (serial process)" begin
        # With MPI not loaded, :auto and false both build a serial Ecosystem; true errors
        # because the MPI extension is absent; a bad value is rejected.
        placement = BySupply(total = 5000)
        @test build_ecosystem(placement, spp, env, seed = 1,
                              distributed = :auto) isa
              Ecosystem
        @test build_ecosystem(placement, spp, env, seed = 1,
                              distributed = false) isa
              Ecosystem
        @test_throws ErrorException build_ecosystem(placement, spp, env,
                                                    seed = 1,
                                                    distributed = true)
        @test_throws ErrorException build_ecosystem(placement, spp, env,
                                                    seed = 1,
                                                    distributed = :nonsense)
    end
end
@testset "multi-trait species and end-to-end multi-regime simulation" begin
    # Projected (BNG) fixtures: simulation needs a metric grid of uniform cells, so a geographic
    # environment is rejected by `build_ecosystem` (see "geographic grids cannot be simulated").
    temp = _bngraster(WorldClim{BioClim}, fill(290.0K, 9, 9))
    rain = _bngraster(WorldClim{BioClim}, fill(50.0mm / day, 9, 9))

    spp = _species(5,
                   tolerance = (ToleranceSpec(Temperature, μ = 290.0K,
                                              σ = 2.0K),
                                ToleranceSpec(Precipitation, μ = 50.0mm / day,
                                              σ = 10.0mm / day)),
                   demand = (DemandSpec(SolarRadiation, 10.0kJ / day),
                             DemandSpec(Precipitation, 5.0Unitful.L / day)),
                   seed = 1)
    @test length(values(spp.tolerance)) == 2
    @test length(values(spp.demand)) == 2

    env = _env((_reg(temp, axis = Temperature),
                _reg(rain, axis = Precipitation)),
               (UniformSpec(10.0kJ / (km^2 * day), axis = SolarRadiation),
                UniformSpec(10.0mm / day, axis = Precipitation)))
    eco = build_ecosystem(BySupply(total = 500), spp, env, seed = 1)
    @test eco.nichefit isa MultiplicativeFit
    @test length(values(eco.nichefit)) == 2
    @test_nowarn simulate!(eco, 2month_mean_duration, 1month_mean_duration)
end
@testset "niche-axis threading + by-axis matching" begin
    # a layer carries the axis of its regime spec
    toy = _area(extent = (5km, 5km), cellsize = 1km)
    env = GridHabitat(regime = GradientSpec(274.0K, 303.0K,
                                            axis = Temperature),
                      supply = SUP, area = toy)
    @test EcoSISTEM.axisof(env.regime) === Temperature
    # a trait carries the axis its spec names
    sp = _species(2,
                  tolerance = ToleranceSpec(Temperature, μ = 288.0K, σ = 5.0K),
                  seed = 1)
    @test sp.tolerance isa NicheTolerance{Temperature}
    @test EcoSISTEM.axisof(sp.tolerance) === Temperature
    # matched axes assemble an ecosystem
    @test build_ecosystem(sp, env, seed = 1) isa Ecosystem
    # a mismatched axis is a clear build-time error (not a silent by-unit match)
    spP = _species(2,
                   tolerance = ToleranceSpec(Precipitation, μ = 50.0mm / day,
                                             σ = 10.0mm / day),
                   seed = 1)
    @test_throws ErrorException build_ecosystem(spP, env, seed = 1)
    # **The root is an axis like any other, and matching is IDENTITY** - so a `Temperature`
    # species against a root-axis environment is refused: `NicheAxis` means *"I claim nothing"*,
    # and a layer that declines to say what it measures must not be silently read as saying
    # whatever the species happens to need.
    rootenv = GridHabitat(regime = GradientSpec(274.0K, 303.0K,
                                                axis = EcoSISTEM.NicheAxis),
                          supply = SUP, area = toy)
    @test_throws ErrorException build_ecosystem(_species(2,
                                                         tolerance = ToleranceSpec(Temperature,
                                                                                   μ = 288.0K,
                                                                                   σ = 5.0K),
                                                         seed = 1),
                                                rootenv, seed = 1)
    # Both sides on the root is identity too, and pairs - but the root declares no density width,
    # so a *continuous* tolerance there is refused where its fit is built, naming the remedy: an
    # axis of one's own beneath it. The tolerance carries its own `support`, since the root
    # declares no canonical unit for a `ToleranceSpec` to build in.
    roottol = NicheTolerance(EcoSISTEM.NicheAxis, Normal, fill(288.0, 2),
                             fill(5.0, 2), support = K)
    rootspp = _species(2, tolerance = roottol, seed = 1)
    @test_throws "declares no `densitywidth`" build_ecosystem(rootspp, rootenv,
                                                              seed = 1)
    @nicheaxis(TestSurveyed<:EcoSISTEM.NicheAxis, condition=K,
               densitywidth=1.0K)
    ownenv = GridHabitat(regime = GradientSpec(274.0K, 303.0K,
                                               axis = TestSurveyed),
                         supply = SUP, area = toy)
    ownspp = _species(2,
                      tolerance = ToleranceSpec(TestSurveyed, μ = 288.0K,
                                                σ = 5.0K),
                      seed = 1)
    @test build_ecosystem(ownspp, ownenv, seed = 1) isa Ecosystem
end
# **An explicit `nichefit` must not bypass the axis check, so this must use the EXPLICIT-`nichefit`
# route.** A caller-supplied `nichefit` walks straight past `_defaultsuitability`, so the inferred
# route - which refuses this anyway - proves nothing about it. And the two axes must **share a
# unit**: `Temperature` and `TemperatureRange` are both `K`, so the `eltype` comparison passes them
# and the axis comparison in `_checkmembers` is the only thing that can catch it. Without it, this
# would assemble an ecosystem pairing a temperature *range* tolerance with an absolute temperature
# regime, silently.
@testset "an explicit nichefit cannot bypass the axis check" begin
    toy = _area(extent = (5km, 5km), cellsize = 1km)
    env = GridHabitat(regime = UniformSpec(290.0K, axis = Temperature),
                      supply = SUP, area = toy)
    sp = _species(2,
                  tolerance = ToleranceSpec(TemperatureRange, μ = 5.0K,
                                            σ = 1.0K),
                  seed = 1)
    # The units really do agree, so nothing but the axis distinguishes these.
    @test eltype(sp.tolerance) === eltype(env.regime)
    @test EcoSISTEM.axisof(sp.tolerance) !== EcoSISTEM.axisof(env.regime)

    fit = NicheSuitability{Temperature, typeof(1.0K)}()
    msg = try
        build_ecosystem(sp, env, nichefit = fit, seed = 1)
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("niche axis", msg)
    @test occursin("TemperatureRange", msg) && occursin("Temperature", msg)
    # ...and the same pairing on one axis assembles, so the refusal is about the axes and not the
    # explicit `nichefit`.
    @test build_ecosystem(_species(2,
                                   tolerance = ToleranceSpec(Temperature,
                                                             μ = 290.0K,
                                                             σ = 2.0K),
                                   seed = 1),
                          env, nichefit = fit, seed = 1) isa Ecosystem
end
@testset "niche-axis threading through the data path" begin
    # the shipped layer table resolves a source+code to its niche axis (no download)
    @test EcoSISTEM._specaxis(SourceSpec(WorldClim{BioClim}, 1)) === Temperature
    @test EcoSISTEM._specaxis(SourceSpec(WorldClim{BioClim}, 12)) ===
          Precipitation
    # ...and a bare `(source, code)` pair is refused, here as everywhere: naming the dataset is what
    # `SourceSpec` is for.
    @test_throws ErrorException EcoSISTEM._specaxis((WorldClim{BioClim}, 1))
    # a bare raster has dropped its code, so no axis
    temp = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5))
    @test EcoSISTEM._specaxis(temp) === EcoSISTEM.NicheAxis
    # an explicit axis on the regime spec tags the data regime; the default leaves it unclassified
    @test EcoSISTEM.axisof(_env(_reg(temp, axis = Temperature),
                                SUP).regime) === Temperature
    @test EcoSISTEM.axisof(_env(_reg(temp), SUP).regime) ===
          EcoSISTEM.NicheAxis
    # a matched-axis species assembles; a mismatched one is a clear build-time error. Assembly needs a
    # projected grid (see "geographic grids cannot be simulated"), so this half uses a BNG fixture -
    # the axis threading being tested is CRS-independent.
    env = _env(_reg(_bngraster(WorldClim{BioClim}, fill(290.0K, 9, 9)),
                    axis = Temperature), SUP)
    sp = _species(3,
                  tolerance = ToleranceSpec(Temperature, μ = 290.0K, σ = 2.0K),
                  demand = DemandSpec(SolarRadiation, 8.0kJ / day), seed = 1)
    @test build_ecosystem(sp, env, seed = 1) isa Ecosystem
    spP = _species(3,
                   tolerance = ToleranceSpec(Precipitation, μ = 50.0mm / day,
                                             σ = 10.0mm / day),
                   demand = DemandSpec(SolarRadiation, 8.0kJ / day), seed = 1)
    @test_throws ErrorException build_ecosystem(spP, env, seed = 1)
end
@testset "DefaultEcosystem fills omitted required inputs" begin
    # Each builder announces every defaulted input with an @info, then delegates to the strict method.
    env = @test_logs (:info,) (:info,) (:info,) (:info,) match_mode=:any build_habitat()
    @test env isa GridHabitat
    @test env.supply isa Supply{SolarRadiation}

    sp = @test_logs (:info,) (:info,) (:info,) match_mode=:any build_species(DefaultEcosystem())
    @test sp isa SpeciesList
    @test length(sp.names) == 10
    @test sp.demand isa Demand{SolarRadiation}

    # An input passed explicitly is not defaulted (and is honoured); only the rest are filled.
    sp3 = build_species(DefaultEcosystem(), numspecies = 3,
                        demand = DemandSpec(Precipitation, 2.0Unitful.L / day))
    @test length(sp3.names) == 3
    @test sp3.demand isa Demand{Precipitation}

    # The full trio assembles into a runnable ecosystem end-to-end.
    eco = build_ecosystem(DefaultEcosystem(), seed = 1)
    @test eco isa Ecosystem
    @test_nowarn simulate!(eco, 2month_mean_duration, 1month_mean_duration)
end
# **`build_habitat`'s first argument says where an unnamed input comes from**, and the mechanism
# is that a keyword default may call a function dispatching on an earlier argument - evaluated at
# call time, so a keyword that *is* named never evaluates its default and never announces.
# That is what makes `build_habitat(h, supply = ...)` a rebuild-with-one-change: the other three
# come from `h`, the grid included.
@testset "build_habitat fills from its source" begin
    toy = build_habitat(verbosity = :silent)
    @test toy isa GridHabitat
    @test size(toy.active) == (10, 10)
    @test toy.topology isa Island

    # An input passed explicitly is honoured, and only the *other* three announce.
    rain = UniformSpec(2.0mm / day, axis = Precipitation)
    wet = @test_logs (:info,) (:info,) (:info,) match_mode = :any build_habitat(supply = rain)
    @test wet.supply isa Supply{Precipitation}

    # --- from an existing habitat -------------------------------------------------
    again = build_habitat(toy, verbosity = :silent)
    @test size(again.active) == size(toy.active)
    @test count(again.active) == count(toy.active)
    @test again.regime.matrix == toy.regime.matrix
    @test again.supply.matrix == toy.supply.matrix
    @test again.topology == toy.topology
    # ...and still **copied**, not aliased: the source habitat must not become live simulation
    # state for the new one. That is `materialise`'s `_ownlayer`, reached through the built-layer
    # path this route always takes.
    @test parent(again.regime.matrix) !== parent(toy.regime.matrix)
    @test parent(again.supply.matrix) !== parent(toy.supply.matrix)

    # ...with one thing changed, and everything else carried across.
    other = build_habitat(toy, supply = rain, verbosity = :silent)
    @test other.regime.matrix == toy.regime.matrix
    @test other.supply isa Supply{Precipitation}
    @test count(other.active) == count(toy.active)
    @test build_habitat(toy, topology = Torus(),
                        verbosity = :silent).topology isa Torus

    # A source that cannot answer says so, rather than raising a `MethodError` on a private
    # helper - the cost of leaving the first argument untyped so a new source is one added method.
    err = try
        build_habitat(:nonsense)
    catch e
        sprint(showerror, e)
    end
    @test occursin("cannot take defaults from a Symbol", err)
end
@testset "a demand on a non-resource axis is refused, not a MethodError" begin
    # `demandtype` returns `nothing` for an axis that declares no resource, so `DemandSpec` must not
    # invoke it directly: that gives `MethodError: objects of type Nothing are not callable` rather
    # than the one wording `_notaresource` owns. `_demandtype` mirrors `_supplytype`, which does the
    # same job on the supply side.
    err = try
        DemandSpec(Temperature, 5.0kJ / day)
        nothing
    catch e
        # The assertion is that it is NOT a `MethodError`, which would name `Nothing`, not the axis.
        @test !(e isa MethodError)
        sprint(showerror, e)
    end
    @test !isnothing(err)
    @test occursin("not a consumable resource", err)
    @test occursin("demand", err)
    # And a genuine resource axis still builds, so the guard has not become a blanket refusal.
    @test _species(4,
                   tolerance = ToleranceSpec(Temperature, μ = 293.0K, σ = 3.0K),
                   movement = MovementSpec(InsectDispersal(), meandist = 1.0km),
                   demand = DemandSpec(SolarRadiation, 5.0kJ / day),
                   seed = 1) isa SpeciesList
end

# Each removed way of building species or an ecosystem raises an error saying how to rewrite the
# call, and the rewrite of a removed `build_species` call is the call itself in the pool's form.
@testset "a removed form says how to rewrite the call" begin
    rewrite = try
        build_species(3, tolerance = TOL, toleranceaxis = Temperature,
                      demand = DEM, demandaxis = SolarRadiation,
                      abundance = 300)
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("`build_species(3; ...)` was removed in v0.11.0", rewrite)
    @test occursin("build_species(SpeciesPool(3); tolerance = " *
                   "ToleranceSpec(Temperature, μ = 298.0 K, σ = 2.0 K)",
                   rewrite)
    @test occursin("build_ecosystem(BySupply(total = 300), species, habitat)",
                   rewrite)
    # A keyword only the removed form took is refused by the `DefaultEcosystem()` shorthand too.
    @test_throws "builds species from a `SpeciesPool` since v0.11.0" build_species(DefaultEcosystem(),
                                                                                   demand = DEM,
                                                                                   demandaxis = SolarRadiation,
                                                                                   verbosity = :silent)
    tol = NicheTolerance(Temperature, Normal, fill(298.0K, 2), fill(2.0K, 2))
    dem = Demand{SolarRadiation}(fill(DEM, 2))
    movement = fill(InsectDispersal()(meandist = 1.0km), 2)
    params = EqualPop(0.6 / year, 0.6 / year, 1.0, 0.2)
    @test_throws "`abun` and `native` arguments were removed in v0.11.0" SpeciesList(2,
                                                                                     tol,
                                                                                     [
                                                                                         10,
                                                                                         10
                                                                                     ],
                                                                                     dem,
                                                                                     movement,
                                                                                     params,
                                                                                     [
                                                                                         true,
                                                                                         true
                                                                                     ])
    spp = SpeciesList(2, tol, dem, movement, params)
    habitat = build_habitat(DefaultEcosystem(), verbosity = :silent)
    fit = NicheSuitability{Temperature, typeof(1.0K)}()
    @test_throws "was removed in v0.11.0" Ecosystem(spp, habitat, fit)
    @test_throws "was removed in v0.11.0" MPIEcosystem(spp, habitat, fit)
    for populate in (populate!, populate_by_tolerance!, emptypopulate!)
        @test_throws "was removed in v0.11.0" populate(nothing, spp, habitat,
                                                       fit, nothing)
    end
end

end
