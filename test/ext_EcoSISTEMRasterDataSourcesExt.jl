# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests for `ext/EcoSISTEMRasterDataSourcesExt`, run as part of `core_ext.jl`:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["core_ext.jl"])'
#
# This file must go through `Pkg.test`, not `julia --project=.. ext_...jl`: `RasterDataSources` is a
# weak dependency and so is not in the package environment.
#
# **What is asserted is the *seam*, not the reading.** Whether a WorldClim layer downloads and
# parses correctly is tested wherever that layer is used; what can only be tested here is that each
# hook the parent declares gains exactly one method when the extension loads, and - the half a
# same-process test cannot see - that the parent still works without it.

using Test
using EcoSISTEM
using EcoSISTEM: layerinfo, SourceSpec
# Only for the deprecated per-source wrappers below, which stay in the submodule.
using EcoSISTEM: ClimatePref
# `public`, not exported, so `using` does not bring them in - they have to be named.
using RasterDataSources
import Rasters
using SimpleTraits
using Unitful
using Unitful.DefaultSymbols
using DimensionalData
using DimensionalData: Ti
using EcoSISTEM.Units
using Plots

include("rasterfixtures.jl")

@testset "EcoSISTEMRasterDataSourcesExt" begin
    @testset "the extension is loaded" begin
        @test !isnothing(Base.get_extension(EcoSISTEM,
                                            :EcoSISTEMRasterDataSourcesExt))
        # Its `__init__` is what defaults the download directory into EcoSISTEM's scratch space.
        @test haskey(ENV, "RASTERDATASOURCES_PATH")
    end

    @testset "hooks gain their sole method" begin
        # The trait is what every door - `SourceSpec`, `ClimateRaster` - actually gates on, and
        # marking the abstract type is what covers the whole hierarchy in one line.
        @test istrait(EcoSISTEM.IsRasterData{WorldClim{BioClim}})
        @test istrait(EcoSISTEM.IsRasterData{EarthEnv{LandCover}})
        @test !istrait(EcoSISTEM.IsRasterData{Int})

        # Derived from `RasterDataSources.layers`, so it cannot drift from it: `Int` for BioClim and
        # land cover, `Symbol` for the rest.
        @test EcoSISTEM._codetype(WorldClim{BioClim}) === Int
        @test EcoSISTEM._codetype(WorldClim{Climate}) === Symbol
        @test EcoSISTEM._codetype(Int) === Nothing

        @test length(EcoSISTEM._alllayercodes(EarthEnv{LandCover})) == 12
        # An unmarked source cannot answer, and says so rather than returning an empty list.
        @test_throws ErrorException EcoSISTEM._alllayercodes(Int)

        # Any spelling a caller may write resolves to the one the dataset uses.
        @test EcoSISTEM._preferredcode(WorldClim{BioClim}, "bio4") == 4
        @test isnothing(EcoSISTEM._preferredcode(WorldClim{BioClim}, nothing))
    end

    @testset "the catalogue stays in the parent, keyed on any Type" begin
        # These are `EcoSISTEM.ClimatePref`'s own, not the extension's - they use their argument
        # only to find a shipped CSV. The extension supplies just the one per-dataset correction.
        @test layerunit(WorldClim{BioClim}, 1) == u"°C"
        @test layerinfo(WorldClim{Climate}, :srad).unit == u"kJ*m^-2"
        @test EcoSISTEM._documentedceiling(WorldClim{BioClim}, 1) == 100.0
        @test EcoSISTEM._documentedceiling(EarthEnv{LandCover}, 1) == 1.0
        # An unknown type is refused by the table lookup, naming the file it looked for.
        @test_throws ErrorException layerunit(Int, 1)
    end

    @testset "SourceSpec is gated by the trait, not by a package" begin
        one = SourceSpec(WorldClim{BioClim}, 1)
        @test one.source === WorldClim{BioClim}
        @test one.code == 1
        # No code means the whole dataset, expanded through `_alllayercodes` at construction.
        @test length(SourceSpec(EarthEnv{LandCover}).code) == 12
        # The tailored error, not a bare `MethodError` leaking `SimpleTraits.Not`.
        @test_throws ErrorException SourceSpec(Int, 1)
    end

    @testset "dataset-typed methods are supplied here" begin
        @test hasmethod(read, Tuple{Type{WorldClim{BioClim}}})
        @test !isempty(methods(EcoSISTEM.sourcecrs))
        @test !isempty(methods(EcoSISTEM.compress_landcover))
        @test !isempty(methods(EcoSISTEM.landcoverclass))
        # Looked up by name in the shipped table, never hardcoded.
        @test EcoSISTEM.landcoverclass(:open_water) isa Int
        # The deprecated wrapper constructors keep their names in the parent and gain methods here.
        @test !isempty(methods(ClimatePref.Worldclim_bioclim))
        # **There is no `_bioclimhabitat`, and no dataset-typed method here for any of the four
        # data-backed `*habitat` builders** - nor a method-less stub in the parent for them, since
        # those builders sampled no grid, which is not how a layer is built. The released *names*
        # live in `src/deprecations.jl` as errors that explain themselves, and are asserted in
        # `test_deprecations.jl` rather than here, because they need nothing from this extension.
        @test !isdefined(EcoSISTEM, :_bioclimhabitat)
    end

    # **`read(::Type{CRUTS}, ...)` is this extension's own reader, and had no test at all.** It is
    # here rather than in the parent because it takes its variable codes and units from
    # `WorldClim{Climate}`'s shipped layer table - CRU TS has none of its own - so it needs the
    # trigger package that defines that dataset type.
    #
    # Synthetic files, so this runs on every platform and downloads nothing: what is asserted is
    # the seam, not CRU TS's own data.
    @testset "read(CRUTS, dir, var) gives a ClimateRaster carrying the source" begin
        dir = mktempdir()
        # Twelve monthly `.tif`s, which is the shape `_readmonthlydir` expects. Built through the
        # shared `_testraster` fixture rather than by hand, so they carry what a real GeoTIFF does -
        # a CRS, and `Intervals(Start)` sampling. A hand-rolled `Rasters.Raster` gets neither, and
        # the package refuses both by name: no CRS gives `LOCAL_CS[...UNIT["unknown"]]`, and `Points`
        # sampling has no cell extent to label.
        # North-up, as a real GeoTIFF is: latitude descends from the top-left corner. Written
        # south-up, GDAL warns on every file - and the files would then exercise an orientation no
        # CRU TS release actually has, which is the opposite of what a fixture is for.
        for m in 1:12
            src = _testraster(WorldClim{Climate}, fill(Float32(m), 5, 5),
                              lat = (4.0:-1.0:0.0) .* °)
            r = Rasters.Raster(parent(src.array), dims(src.array))
            Rasters.write(joinpath(dir, "cruts_$(lpad(m, 2, '0')).tif"), r,
                          force = true)
        end
        cr = read(CRUTS, dir, "tavg")

        # **A raster carrying the source in its parameter, not a container type of its own** - the
        # shape every other data source already had, and what makes `CRUTS` nameable as a source.
        @test cr isa ClimateRaster
        @test typeof(cr).parameters[1] === EcoSISTEM.CRUTS
        # Three dimensions, the third being time: the guard `_timeseriesraster` applies, and every
        # reader and recipe for this source slices by it.
        @test ndims(cr.array) == 3
        @test size(cr.array, 3) == 12
        # **The unit comes from `WorldClim{Climate}`'s table for the same short code**, which is the
        # whole reason this reader lives in the extension: CRU TS has no layer table of its own.
        # Compared on `dimension`, because `tavg`'s table entry is an *affine* temperature (a
        # position on a scale) while the values arrive on the absolute one - `unit` has no method
        # for an affine `FreeUnits` at all, which is what makes the dimension the honest comparison.
        @test dimension(eltype(cr.array)) ==
              dimension(layerunit(WorldClim{Climate}, "tavg"))
    end

    # **The monthly-climate plot recipes belong to this extension**, because they dispatch on
    # `WorldClim{Climate}`/`CHELSA{Climate}` - types only this extension's trigger provides. They
    # cannot live in a `Plots` extension: one may see only its own trigger and the parent's
    # dependencies, so it could not name those types at all.
    @testset "monthly climate rasters plot, and name the right month" begin
        # Built the way `_mkstackaxis` really builds it - `Ti((1:12) .* month_mean_duration)`,
        # 1-based. A hand-built `Ti(January:1month_mean_duration:December)` is **0-based**, which no
        # reader produces, and it hides a month bug in the recipe: on a real 1-based axis
        # `At(0month)` matches no slice at all, while a value that does match is labelled one month
        # late. Each slice holds its own month number, so a test can tell *which* slice was plotted
        # - the thing a constant-filled array cannot show.
        temp = DimArray([1.0K * k for i in 1:11, j in 1:11, k in 1:12],
                        (Y((0°):(1°):(10°)), X((0°):(1°):(10°)),
                         Ti((1:12) .* month_mean_duration)))
        wm = ClimateRaster(WorldClim{Climate}, temp)
        # A month is an ordinal, so it is named by an ordinal. `January`...`December` come from
        # `EcoSISTEM.Units` (re-exported from `Dates`), so naming a month needs no `using Dates` - which
        # matters because `Dates` and `Units` both export `day`/`month`/`week`/`year`.
        @test plot(wm, January).n == 1
        @test plot(wm, September, 0° .. 2°, 5° .. 10°).n == 1
        chelsa = ClimateRaster(CHELSA{Climate}, temp)
        @test plot(chelsa, February).n == 1
        @test plot(chelsa, December, 3° .. 10°, 0° .. 3°).n == 1

        # A duration is still accepted, for a caller who has one in hand - including one in an
        # equal-but-differently-named unit, which needs the `uconvert` in `_monthindex` to cancel.
        @test plot(wm, 1month_mean_duration).n == 1
        @test plot(wm, 1month_mean_duration).n == 1

        # The labels are what actually pin the direction: `January` is January, not February. Under the
        # old `(time + 1month) / month` these came out one month late, and no test noticed.
        # `_monthindex` itself is tested in `test_Units.jl`, where it now lives; what matters here
        # is that the *recipe* uses it, which the titles below show.
        @test plot(wm, January).subplots[1][:title] == "Jan"
        @test plot(wm, December).subplots[1][:title] == "Dec"
        # And the data really is that month's: slice `n` is filled with `n` K, so the plotted series
        # for `March` must be all 3s. This is the assertion the old 0-based axis made impossible.
        @test all(==(3.0), plot(wm, March).series_list[1][:z].surf)

        # Out-of-range and fractional months are refused, with a message naming the months rather than
        # the old `error("NO")`.
        @test_throws ErrorException plot(wm, 0)
        @test_throws ErrorException plot(wm, 13)
        @test_throws ErrorException plot(wm, 0month_mean_duration)
        @test_throws ErrorException plot(wm, 1.5month_mean_duration)
    end
end
