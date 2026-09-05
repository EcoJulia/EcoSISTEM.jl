# SPDX-License-Identifier: LGPL-3.0-or-later

module TestClimate

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using DimensionalData
using RasterDataSources
using Test
using Plots
# Named, not a blanket `using Dates`: it would collide with `EcoSISTEM.Units` on `day`/`week`/`year`.
using Dates: Dates, Date
using SimpleTraits: SimpleTraits, istrait

# The three netCDF archives are **sources**, not containers: a reader returns a `ClimateRaster`
# carrying the source in its parameter, exactly as it does for a RasterDataSources dataset. That is
# what lets one name a source at all - a container type cannot.
@testset "ECMWF sources build rasters, and satisfy IsRasterData" begin
    temp = DimArray(fill(1.0K, 10, 10, 3),
                    (Y(1:10), X(1:10), Ti(collect(1:3) .* s)))
    for S in (EcoSISTEM.ERA, EcoSISTEM.CERA, EcoSISTEM.CRUTS)
        @test istrait(EcoSISTEM.IsRasterData{S})
        r = EcoSISTEM._timeseriesraster(S, temp)
        @test r isa ClimateRaster
        @test typeof(r).parameters[1] === S
    end
end

# A **`Dates.Date`** axis. All three constructors share `_istimeaxis`, and it must be typed
# `Dates.TimeType` rather than `Dates.AbstractDateTime` - which excludes `Date`, since Julia branches
# `TimeType -> {AbstractDateTime -> DateTime, Date, Time}`. So a daily series could not be built at all,
# and neither could any `CFTime` calendar (`AbstractCFDateTime <: Dates.TimeType`), which is what
# NCDatasets produces from a netCDF whose `calendar` attribute is not Gregorian.
# `Date` is tested rather than a CF type deliberately: it exercises exactly the widened predicate
# and needs no new dependency, where `CFTime` is only an indirect one (via `NCDatasets`).
@testset "a Date axis is accepted - TimeType, not AbstractDateTime" begin
    @test Date <: Dates.TimeType
    @test !(Date <: Dates.AbstractDateTime)      # the reason the old predicate refused it
    dated = DimArray(fill(1.0K, 10, 10, 3),
                     (Y(1:10), X(1:10),
                      Ti([Date(2000, 1, 1), Date(2000, 1, 2), Date(2000, 1, 3)])))
    for S in (EcoSISTEM.ERA, EcoSISTEM.CERA, EcoSISTEM.CRUTS)
        @test_nowarn EcoSISTEM._timeseriesraster(S, dated)
    end
    # ...and the predicate is a widening rather than a removal: a third axis that is not time still
    # fails. `ClimateRaster` itself has no such guard, so `_timeseriesraster` is where it lives and
    # the readers are the only door to it.
    untimed = DimArray(fill(1.0K, 10, 10, 3),
                       (Y(1:10), X(1:10), Ti(collect(1:3))))
    @test_throws ErrorException EcoSISTEM._timeseriesraster(EcoSISTEM.ERA,
                                                            untimed)
end

@testset "Worldclim types" begin
    temp = DimArray(fill(1.0K, 10, 10, 12),
                    (Y(1:10), X(1:10), Ti(collect(1:12) .* s)))
    @test_nowarn ClimateRaster(WorldClim{Climate}, temp)
    @test_nowarn ClimateRaster(CHELSA{Climate}, temp)
end

@testset "Bioclim types" begin
    temp = DimArray(fill(1.0K, 10, 10, 19),
                    (Y(1:10), X(1:10), Dim{:vars}(collect(1:19))))
    @test_nowarn ClimateRaster(WorldClim{BioClim}, temp)
end

# `SourceSpec.code` is never `nothing`. `SourceSpec(dataset)` resolves the dataset's own code list
# at construction, so there is one code shape instead of three ("one", "several", "all of them") and
# every layer's identity is known before anything is read - which is what lets each keep its own unit.
@testset "SourceSpec resolves its code list at construction" begin
    single = SourceSpec(WorldClim{BioClim}, :bio1)
    @test single.code === :bio1
    @test single.unit == u"°C"

    # A whole-dataset spec expands, rather than carrying a late-bound sentinel.
    whole = SourceSpec(EarthEnv{LandCover})
    @test whole.code isa Vector
    @test length(whole.code) == 12
    # ...and it may honestly claim the unit its layers agree on.
    @test whole.unit == EcoSISTEM.layerunit(EarthEnv{LandCover},
                              first(whole.code))

    # An explicit codes vector is the middle ground between one layer and all of them.
    pair = SourceSpec(WorldClim{Climate}, [:tmin, :tmax])
    @test pair.code == [:tmin, :tmax]
    @test pair.unit == u"°C"          # they agree, so the spec can say so

    # Layers that disagree get a neutral *placeholder*, not a claim - four of the seven shipped
    # datasets are heterogeneous, so refusing them here would rule out the flagship sources.
    # `WorldClim{BioClim}` alone mixes six units.
    mixed = SourceSpec(WorldClim{BioClim})
    @test length(mixed.code) == 19
    @test mixed.unit == NoUnits
    @test EcoSISTEM._sharedunit(WorldClim{BioClim}, mixed.code) ==
          NoUnits
end

@testset "ConstructedRasterSpec declares when its combine runs" begin
    # Nothing is read here - the stage is a property of the spec, resolved at construction.
    @test ConstructedRasterSpec(() -> nothing, axis = EcoSISTEM.NicheAxis).combinestage isa
          CombineOnTargetGrid
    early = ConstructedRasterSpec(() -> nothing, axis = EcoSISTEM.NicheAxis,
                                  combinestage = CombineOnSourceGrid())
    @test early.combinestage isa CombineOnSourceGrid
    # A type, not a symbol, so a wrong stage is refused by the signature where it was written
    # rather than by a check that runs later and somewhere else.
    @test_throws TypeError ConstructedRasterSpec(() -> nothing,
                                                 axis = EcoSISTEM.NicheAxis,
                                                 combinestage = :early)
    # It is independent of the **axis**: an early combine that produces ordinary numbers (a ratio,
    # say) is on no particular axis, and `combinestage` says nothing about what the values are.
    # There is deliberately no `valuetype` field to assert on: whether a layer holds class codes is
    # the axis's to say.
    @test early.axis === EcoSISTEM.NicheAxis
    @test !EcoSISTEM.iscategorical(early.axis)
end

@testset "ShapeSpec: URL vs local path" begin
    # A leading URL scheme defers to a `CachedAsset` (downloaded only at materialisation); an
    # already-local path is kept as-is - neither variant touches the filesystem or network here.
    @test ShapeSpec("https://example.com/shape.zip").path isa
          EcoSISTEM.CachedAsset
    @test ShapeSpec("http://example.com/shape.zip").path isa
          EcoSISTEM.CachedAsset
    @test ShapeSpec("/local/path/shape.zip").path isa String
    @test ShapeSpec("relative/shape.zip").path isa String
end

# **A raster broadcasts and stays a raster**, which is what lets a `ConstructedRasterSpec` combine name
# no array type - neither `.array` going in nor a constructor coming out. Without it the combine
# contract could only be met by hand-wrapping, which puts our array type in *user* code.
@testset "a raster broadcasts and stays a raster" begin
    yx = (Y(1:2), X(1:2))
    a = ClimateRaster(EarthEnv{LandCover}, DimArray(Float64[1 2; 3 4], yx), 1)
    b = ClimateRaster(EarthEnv{LandCover}, DimArray(Float64[10 20; 30 40], yx),
                      2)
    same = ClimateRaster(EarthEnv{LandCover},
                         DimArray(Float64[10 20; 30 40], yx), 1)

    @test a .+ b isa ClimateRaster
    @test parent((a .+ b).array) == Float64[11 22; 33 44]
    # `sum` over the varargs of a multi-band combine is the shape that matters most, and it
    # reduces with `+` - so this is the whole of `sum(bands)` working.
    @test parent(sum((a, b)).array) == Float64[11 22; 33 44]
    # A predicate gives a `Bool`-valued raster: that *is* a mask, not a different container.
    @test eltype(a .!= 1.0) == Bool
    @test (a .!= 1.0) isa ClimateRaster

    # The dims survive, which is what the result is later sampled and cropped by.
    ydim(r) = DimensionalData.dims(r.array, Y)
    @test ydim(a .+ b) == ydim(a)

    # **The metadata rule: a derived raster inherits NOTHING.** Not the code, not the value type,
    # not the source. A combine is free to change all three - `argmax` over twelve *continuous*
    # land-cover bands gives a *class code* - so even unanimous inputs say nothing about the output.
    # An earlier version kept whatever every input agreed on, which is the same reasoning that put
    # a false source on every derived raster. What a derived layer *is* comes from its spec's `axis`.
    @test isnothing((a .* 2).code)        # one input, and still no inherited code
    @test isnothing((a .+ same).code)     # two that agree - agreement is not the test
    @test isnothing((a .+ b).code)        # two that do not
end

# A source need not be a `RasterDataSource`: what may name one is a **trait**, so a type this
# package has never seen becomes usable with one `@traitimpl` line and no change here.
struct AForeignSource end
SimpleTraits.@traitimpl EcoSISTEM.IsRasterData{AForeignSource}
struct AnUnmarkedType end

@testset "which types may name a source" begin
    @test istrait(EcoSISTEM.IsRasterData{WorldClim{BioClim}})      # marked via the RDS supertype
    @test istrait(EcoSISTEM.IsRasterData{EcoSISTEM.SyntheticData})        # ...and via ours
    @test istrait(EcoSISTEM.IsRasterData{AForeignSource})          # ...and one marked from outside
    @test !istrait(EcoSISTEM.IsRasterData{AnUnmarkedType})

    yx = (Y(1:2), X(1:2))
    arr = DimArray(Float64[1 2; 3 4], yx)
    @test_nowarn ClimateRaster(AForeignSource, arr)
    # A source with no `codetype` cannot be given a code - it has no layers to name.
    @test_throws ErrorException ClimateRaster(EcoSISTEM.SyntheticData, arr,
                                              :bio4)
end

@testset "a layer code is stored in one canonical spelling" begin
    yx = (Y(1:2), X(1:2))
    arr = DimArray(Float64[1 2; 3 4], yx)
    # All three spellings resolve to one layer in the catalogue, so they must normalise to one
    # code: otherwise two rasters of the same layer hold different codes and compare unequal.
    bysym = ClimateRaster(WorldClim{BioClim}, arr, :bio4)
    byint = ClimateRaster(WorldClim{BioClim}, arr, 4)
    bystr = ClimateRaster(WorldClim{BioClim}, arr, "bio4")
    @test bysym.code == byint.code == bystr.code
    @test typeof(bysym) === typeof(byint) === typeof(bystr)

    # The canonical spelling is the dataset's OWN, from `RasterDataSources.layers` - an `Int` here,
    # a `Symbol` for a source that names its layers that way. It is deliberately not uniform.
    @test bysym.code isa Int
    @test ClimateRaster(WorldClim{Climate}, arr, "tmin").code === :tmin
    # ...including where the two spellings differ only by case, which `layers` settles.
    @test ClimateRaster(EarthEnv{HabitatHeterogeneity}, arr, :contrast).code ===
          :Contrast

    # A code naming no layer is refused where it is written, not left to fail later.
    @test_throws ErrorException ClimateRaster(WorldClim{BioClim}, arr, :nope)
    # ...as is something that could not be a code at all.
    @test_throws ErrorException ClimateRaster(WorldClim{BioClim}, arr, 1.5)
    # Mixed spellings in one stack normalise element by element.
    @test ClimateRaster(WorldClim{BioClim}, arr, [:bio4, 3]).code == [4, 3]
end

@testset "derived and synthetic rasters name themselves honestly" begin
    yx = (Y(1:2), X(1:2))
    a = ClimateRaster(EarthEnv{LandCover}, DimArray(Float64[1 2; 3 4], yx), 1)
    src(r) = typeof(r).parameters[1]

    # A combine derives, whether or not its inputs agree - the result is a new quantity whose
    # meaning comes from the spec's `axis`, not from what went in.
    @test src(a .* 2) === EcoSISTEM.DerivedData{EarthEnv{LandCover}}
    # ...and carries no code: inheriting the parent's would attach that layer's whole catalogue row
    # - its unit, axis and accumulation period - to a quantity none of them describe.
    @test isnothing((a .* 2).code)

    # **Nesting collapses at every depth.** The collapse is `derivedfrom`'s job and happens at the
    # *type* level: a normalising constructor would never run, since `DerivedData` is only ever a
    # type parameter and is never instantiated. Measured before that was so - three chained
    # broadcasts gave three levels of wrapper, silently, each one a perfectly good type.
    r = a
    for _ in 1:4
        r = r .* 2
        @test src(r) === EcoSISTEM.DerivedData{EarthEnv{LandCover}}
    end
    @test EcoSISTEM._derivedfrom(EcoSISTEM.DerivedData{EarthEnv{LandCover}}) ===
          EcoSISTEM.DerivedData{EarthEnv{LandCover}}
    # The *type* does not self-normalise - only `derivedfrom` does, which is why nothing should
    # name `DerivedData{...}` directly.
    @test EcoSISTEM.DerivedData{EcoSISTEM.DerivedData{EarthEnv{LandCover}}} !==
          EcoSISTEM.DerivedData{EarthEnv{LandCover}}
    # Deriving from something that could never have been a source is a mistake worth naming.
    @test_throws ErrorException EcoSISTEM._derivedfrom(AnUnmarkedType)
end

@testset "a combine over several sources records all of them" begin
    yx = (Y(1:2), X(1:2))
    lc = ClimateRaster(EarthEnv{LandCover}, DimArray(Float64[1 2; 3 4], yx), 1)
    cl = ClimateRaster(WorldClim{Climate}, DimArray(Float64[5 6; 7 8], yx),
                       :tmin)
    src(r) = typeof(r).parameters[1]

    # **Order-independent, and that is the whole point.** Taking the *first* input's source made
    # `lc .+ cl` and `cl .+ lc` different types for the same computation, and silently discarded the
    # other lineage - the same defect, one level up, as a derived raster claiming its parent's
    # identity. The lineage is sorted, so the two agree.
    @test typeof(lc .+ cl) === typeof(cl .+ lc)
    @test src(lc .+ cl) ===
          EcoSISTEM.DerivedData{Tuple{EarthEnv{LandCover}, WorldClim{Climate}}}

    # Deriving again neither nests nor grows: the lineage flattens, and a source already in it is
    # not added twice.
    @test src((lc .+ cl) .* 2) === src(lc .+ cl)
    @test src((lc .+ cl) .+ lc) === src(lc .+ cl)
    # ...while a genuinely new source does extend it.
    @test length(EcoSISTEM._origins(src((lc .+ cl) .+
                                        ClimateRaster(WorldClim{BioClim},
                                                      DimArray(Float64[1 1;
                                                                       1 1],
                                                               yx),
                                                      4)))) == 3
end

# **The ERA/CERA plot recipes live in `src/Climate.jl`**, beside the types they plot, not in an extension:
# a `@recipe` needs only `RecipesBase`, a hard dependency, so they are defined unconditionally and are
# inert until a backend loads. These assertions still need `Plots` to render with, which is why
# this file loads it; `Plots` is in the test target, and is not a weak dependency.
@testset "ERA/CERA plot recipes" begin
    temp = DimArray(fill(1.0K, 11, 11, 21),
                    (Y((0°):(1°):(10°)), X((0°):(1°):(10°)),
                     Ti((2000year):(1year):(2020year))))
    era = EcoSISTEM._timeseriesraster(EcoSISTEM.ERA, temp)
    @test plot(era, 2000year).n == 1
    @test plot(era, 2002year, 1° .. 4°, 5° .. 10°).n == 1
    cera = EcoSISTEM._timeseriesraster(EcoSISTEM.CERA, temp)
    @test plot(cera, 2001year).n == 1
    @test plot(cera, 2020year, 1° .. 4°, 5° .. 10°).n == 1
end

end
