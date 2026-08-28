# SPDX-License-Identifier: LGPL-3.0-or-later

module TestNicheFit

using EcoSISTEM
# `[C7-VIS]` B1/B2/B3: these are `public` rather than exported, so they must be named.
using EcoSISTEM: getdist, getpref
# `nichefitcombine` is `public`, not exported.
using EcoSISTEM: nichefitcombine
using EcoSISTEM: getregime
using EcoSISTEM.Units
using Test
using Distributions
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM: hasdata, landcoverclass
using EcoSISTEM: materialise
using Unitful, Unitful.DefaultSymbols
using DimensionalData: DimensionalData, DimArray, X, Y, dims
using Rasters
using RasterDataSources
using ArchGDAL
using Distributions: Normal
using Extents: Extent
include("rasterfixtures.jl")
include("buildfixtures.jl")

@testset "Trait relationships" begin
    @test_nowarn NicheSuitability{EcoSISTEM.NicheAxis, Unitful.Temperature}()
    @test_nowarn CategoricalSuitability{EcoSISTEM.NicheAxis, Int64}()
    @test_nowarn NoFitContinuous{EcoSISTEM.NicheAxis, Int64}()
    @test_nowarn NoFitCategorical{EcoSISTEM.NicheAxis, Int64}()

    # `NicheSuitability` evaluates a distribution's density at the (unit-stripped) regime value.
    @test NicheSuitability{EcoSISTEM.NicheAxis, Unitful.Temperature}()(Normal(1.0,
                                                                              0.01),
                                                                       1.0K) >
          0.0
    @test NicheSuitability{EcoSISTEM.NicheAxis, typeof(1.0mm)}()(Uniform(1, 2),
                                                                 1.0mm) == 1.0
    @test NicheSuitability{EcoSISTEM.NicheAxis, Int64}()(Trapezoid(1, 2, 3, 4),
                                                         1) == 0.0

    # A concrete `NF` forces a real `uconvert`, not a bare strip: a dimensionally-compatible but
    # differently-scaled `current` is corrected rather than silently misread (1000.0μm ≡ 1.0mm).
    @test NicheSuitability{EcoSISTEM.NicheAxis, typeof(1.0mm)}()(Uniform(1, 2),
                                                                 1000.0Unitful.μm) ==
          1.0
    @test_throws Unitful.DimensionError NicheSuitability{EcoSISTEM.NicheAxis,
                                                         typeof(1.0mm)}()(Uniform(1,
                                                                                  2),
                                                                          1.0K)
    @test EcoSISTEM.iscontinuous(NicheSuitability{EcoSISTEM.NicheAxis,
                                                  Unitful.Temperature}()) ==
          true
    @test eltype(NicheSuitability{EcoSISTEM.NicheAxis, Unitful.Temperature}()) ==
          Unitful.Temperature

    @test CategoricalSuitability{EcoSISTEM.NicheAxis, Int64}()(1, 1) == 1.0
    @test EcoSISTEM.iscontinuous(CategoricalSuitability{EcoSISTEM.NicheAxis,
                                                        Int64}()) == false
    @test eltype(CategoricalSuitability{EcoSISTEM.NicheAxis, Int64}()) == Int64

    @test NoFitContinuous{EcoSISTEM.NicheAxis, Int64}()(1, 1, 1) == 1.0
    @test EcoSISTEM.iscontinuous(NoFitContinuous{EcoSISTEM.NicheAxis, Int64}()) ==
          true
    @test eltype(NoFitContinuous{EcoSISTEM.NicheAxis, Int64}()) == Int64

    @test NoFitCategorical{EcoSISTEM.NicheAxis, Int64}()(1, 1) == 1.0
    @test EcoSISTEM.iscontinuous(NoFitCategorical{EcoSISTEM.NicheAxis, Int64}()) ==
          false
    @test eltype(NoFitCategorical{EcoSISTEM.NicheAxis, Int64}()) == Int64

    # **Named**: every fit here is on the root axis, so derived names cannot tell them apart —
    # what is under test is the *combining*, not the naming.
    tr2 = MultiplicativeFit((a = NoFitContinuous{EcoSISTEM.NicheAxis, Int64}(),
                             b = NoFitCategorical{EcoSISTEM.NicheAxis, Int64}()))
    @test map(EcoSISTEM.iscontinuous, values(tr2)) == (true, false)
    @test map(eltype, values(tr2)) == (Int64, Int64)
    tr3 = MultiplicativeFit((a = NoFitContinuous{EcoSISTEM.NicheAxis, Int64}(),
                             b = NoFitCategorical{EcoSISTEM.NicheAxis, Int64}(),
                             c = NicheSuitability{EcoSISTEM.NicheAxis,
                                                  Unitful.Temperature}()))
    @test map(EcoSISTEM.iscontinuous, values(tr3)) == (true, false, true)
    @test map(eltype, values(tr3)) == (Int64, Int64, Unitful.Temperature)

    # `nichefitcombine` is now a whole-tuple function of the per-layer results, not a binary operator
    @test EcoSISTEM.nichefitcombine(tr2) == prod
    @test EcoSISTEM.nichefitcombine(tr3) == prod
    @test EcoSISTEM.nichefitcombine(tr2)((a = 2.0, b = 3.0)) == 6.0

    tr2 = AdditiveFit((a = NoFitContinuous{EcoSISTEM.NicheAxis, Int64}(),
                       b = NoFitCategorical{EcoSISTEM.NicheAxis, Int64}()))
    @test map(EcoSISTEM.iscontinuous, values(tr2)) == (true, false)
    @test map(eltype, values(tr2)) == (Int64, Int64)
    tr3 = AdditiveFit((a = NoFitContinuous{EcoSISTEM.NicheAxis, Int64}(),
                       b = NoFitCategorical{EcoSISTEM.NicheAxis, Int64}(),
                       c = NicheSuitability{EcoSISTEM.NicheAxis,
                                            Unitful.Temperature}()))
    @test map(EcoSISTEM.iscontinuous, values(tr3)) == (true, false, true)
    @test map(eltype, values(tr3)) == (Int64, Int64, Unitful.Temperature)

    @test EcoSISTEM.nichefitcombine(tr2) == sum
    @test EcoSISTEM.nichefitcombine(tr3) == sum
    @test EcoSISTEM.nichefitcombine(tr2)((a = 2.0, b = 3.0)) == 5.0
end

# --- the suitability functions, which live in `NicheFit.jl` ---

grid = (10, 10)
area = 25.0km^2
totalK = 10000.0kJ / km^2 / day

# Grid decided first: 25 km² over 10 × 10 is 0.5 km cells.
studyarea = StudyArea(extent = (sqrt(area), sqrt(area)),
                      cellsize = sqrt(area) / grid[1], verbosity = :silent)
supply = UniformSpec(totalK, axis = SolarRadiation)

@testset "Trait functions" begin
    # Every regime in this file is scaffolding. What is under test is `_suitability` — the step
    # that reads a cell out of a regime and pairs it with a species' tolerance row — and all it
    # needs is known values on a known axis, so the simplest spec that gives them is the right one.
    habitat1 = GridHabitat(regime = UniformSpec(1.0,
                                                axis = EcoSISTEM.NicheAxis),
                           supply = supply, area = studyarea)
    habitat2 = GridHabitat(regime = GradientSpec(-10.0K, 10.0K,
                                                 axis = Temperature),
                           supply = supply, area = studyarea)

    regime = LayerCollection((habitat1.regime, habitat2.regime))
    tolerance = SpeciesRequirementCollection((NicheTolerance(EcoSISTEM.NicheAxis,
                                                             Normal,
                                                             fill(1.0, 10),
                                                             fill(0.1, 10)),
                                              NicheTolerance(Temperature,
                                                             Normal,
                                                             fill(1.0K, 10),
                                                             fill(0.1K, 10))))
    nichefit = MultiplicativeFit((NicheSuitability{EcoSISTEM.NicheAxis,
                                                   Float64}(),
                                  NicheSuitability{Temperature,
                                                   Unitful.Temperature}()))
    @test_nowarn EcoSISTEM._suitability(regime, tolerance, nichefit, 1, 1)
    # Members by name — **the axis names**, since neither side was named by hand and the two
    # axes here (the root and `Temperature`) are distinguishable. It read `.one`/`.two` until
    # 2026-08-18; those were the positional fallback, which now fires for nothing at all.
    @test tolerance.NicheAxis === values(tolerance)[1]
    @test tolerance.Temperature === values(tolerance)[2]
    @test nichefit.NicheAxis === values(nichefit)[1]
    @test NamedTuple(nichefit) ==
          (NicheAxis = nichefit.NicheAxis, Temperature = nichefit.Temperature)
    # the deprecated symbol-keyed accessors still reach the same members
    @test (@test_deprecated getpref(tolerance, :NicheAxis)) ===
          tolerance.NicheAxis
    @test (@test_deprecated getregime(regime, :Temperature)) ===
          regime.Temperature

    regime = GridHabitat(regime = UniformSpec(1.0K, axis = Temperature),
                         supply = supply, area = studyarea).regime
    @test_nowarn EcoSISTEM._suitability(regime, tolerance.Temperature,
                                        nichefit.Temperature, 1, 1)
    tolerance = TempTolerance(Array(hcat(fill(collect(1:4), 10)...)'))
    nichefit = NicheSuitability{Temperature, Unitful.Temperature}()
    @test_nowarn EcoSISTEM._suitability(regime, tolerance, nichefit, 1, 1)
    @test getpref(tolerance, 1) == params(getdist(tolerance, 1)) ==
          (1.0, 2.0, 3.0, 4.0)

    # The fixture this replaces held `1.0mm` — a **length**, where `Precipitation`'s canonical
    # unit is a **rate**. Only the deprecated builder allowed that, by taking its raster verbatim;
    # measured, `GridHabitat` refuses the length outright with a `DimensionError`. So the
    # regime is a rate now and the nichefit is typed to match. The tolerance's own parameters are
    # unitless — `_suitability` strips before evaluating — so both assertions are unchanged.
    regime = GridHabitat(regime = UniformSpec(1.0mm / day,
                                              axis = Precipitation),
                         supply = supply, area = studyarea).regime
    # a Uniform response takes 2 parameters, so each species' row has 2 entries
    tolerance = RainTolerance(Array(hcat(fill(collect(1:2), 10)...)'))
    nichefit = NicheSuitability{Precipitation, typeof(1.0mm / day)}()
    @test_nowarn EcoSISTEM._suitability(regime, tolerance, nichefit, 1, 1)
    @test getpref(tolerance, 1) == params(getdist(tolerance, 1)) == (1.0, 2.0)
end

@testset "NicheTolerance trait function is distribution-generic + type-stable" begin
    regime = GridHabitat(regime = UniformSpec(1.0K, axis = Temperature),
                         supply = supply, area = studyarea).regime

    # A Normal response (neither Trapezoid nor Uniform) proves any continuous distribution works.
    bin = NicheTolerance(Temperature, Normal,
                         Array(hcat(fill([1.0, 2.0], 10)...)'))
    nichefit = NicheSuitability{Temperature, typeof(1.0K)}()
    @test EcoSISTEM._suitability(regime, bin, nichefit, 1, 1) ==
          pdf(Normal(1.0, 2.0), 1.0)
    # hot-loop type stability (hard constraint): the `D(row...)` splat must not leak Any/Union.
    @test @inferred(EcoSISTEM._suitability(regime, bin, nichefit, 1, 1)) isa
          Float64

    # the nichefit's unit `NF` is imputed from the trait's axis / the regime — not typed by hand
    @test typeof(NicheSuitability(bin)) ==
          NicheSuitability{Temperature, typeof(1.0K)}
    @test eltype(NicheSuitability(bin)) == eltype(bin)
    @test EcoSISTEM._defaultsuitability(bin, regime) isa
          NicheSuitability{Temperature, typeof(1.0K)}
end

@testset "_defaultsuitability derives NF from the tolerance for every tolerance kind" begin
    # Generalises the `NicheTolerance` fix: a `SimpleCategoricalTolerance` must also take its
    # nichefit's `NF` from the tolerance, not `eltype(regime)` — otherwise nichefit trivially mirrors
    # whatever the regime happens to be, and a genuine tolerance/regime disagreement goes uncaught.
    disc = SimpleCategoricalTolerance(fill(1, 5), axis = EcoSISTEM.NicheAxis)
    # Any categorical regime will do here — what is under test is which nichefit a tolerance
    # picks, not the regime's contents.
    regime = GridHabitat(regime = NicheSpec(3, axis = EcoSISTEM.NicheAxis),
                         supply = UniformSpec(10000.0kJ / km^2 / day,
                                              axis = SolarRadiation),
                         area = StudyArea(extent = (5.0km, 5.0km),
                                          cellsize = 1.0km,
                                          verbosity = :silent)).regime
    @test EcoSISTEM._defaultsuitability(disc, regime) isa
          CategoricalSuitability{EcoSISTEM.axisof(disc), eltype(disc)}

    lc = SimpleCategoricalTolerance([[0.2, 0.8] for _ in 1:5],
                                    axis = EcoSISTEM.NicheAxis)
    @test EcoSISTEM._defaultsuitability(lc, regime) isa
          CategoricalSuitability{EcoSISTEM.axisof(lc), eltype(lc)}

    # **The merge's real payoff, and it could not be written before.** Neither
    # `CategoricalTolerance` nor `LandCoverTolerance` carried a niche axis, so a categorical
    # tolerance could never be checked against the regime it was paired with — only the continuous
    # branch was. The merged type carries one, so the same mismatch is now refused on both.
    # `NicheSpec(3)` above is declared on the ROOT axis, and so are the tolerances above — which
    # is why they are accepted: matching is **identity**, and root against root is identity. The
    # root is not a wildcard; a root regime against a `LandCoverTypology` tolerance is refused just
    # as firmly as any other mismatch.
    typed = GridHabitat(regime = NicheSpec(3, axis = LandCoverTypology),
                        supply = UniformSpec(10000.0kJ / km^2 / day,
                                             axis = SolarRadiation),
                        area = StudyArea(extent = (5.0km, 5.0km),
                                         cellsize = 1.0km,
                                         verbosity = :silent)).regime
    @test EcoSISTEM._defaultsuitability(SimpleCategoricalTolerance(fill(1, 5),
                                                                   axis = LandCoverTypology),
                                        typed) isa CategoricalSuitability
    # …and a tolerance on a *different* declared axis is refused, with the same tailored message
    # the continuous branch gives.
    @test_throws ErrorException EcoSISTEM._defaultsuitability(SimpleCategoricalTolerance(fill(1,
                                                                                              5),
                                                                                         axis = Temperature),
                                                              typed)
end

end
