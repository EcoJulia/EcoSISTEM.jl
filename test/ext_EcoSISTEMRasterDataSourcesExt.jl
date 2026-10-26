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

    @testset "every file a dataset can give has a title of its own" begin
        RDS = RasterDataSources
        # The keywords that tell one dataset's files apart, as `getraster` takes them.
        kwsets(T) = T <: WorldClim{Climate} ?
                    [(res = r, month = m) for r in RDS.resolutions(T)
                     for m in 1:12] :
                    T <: Union{WorldClim, EarthEnv{HabitatHeterogeneity}} ?
                    [(res = r,) for r in RDS.resolutions(T)] :
                    T <: CHELSA{Climate} ? [(month = m,) for m in 1:12] :
                    [NamedTuple()]
        paths = Dict{Type, Set{String}}()
        for T in (WorldClim{BioClim}, WorldClim{Climate}, WorldClim{Elevation},
            CHELSA{BioClim}, CHELSA{BioClimPlus}, CHELSA{Climate},
            EarthEnv{LandCover}, EarthEnv{LandCover{:DISCover}},
            EarthEnv{HabitatHeterogeneity})
            rec = EcoSISTEM.datasetinfo(T)
            titles = Dict{String, String}()
            # The layers both the source package and the catalogue know, which are those a spec
            # can name; a table may serve several datasets, each with only some of its rows.
            tabled = Set(a
                         for r in EcoSISTEM._catalogue()
                         if r.dataset === EcoSISTEM._tablename(T)
                         for a in r.aliases)
            codes = [c for c in RDS.layers(T) if string(c) in tabled]
            @test !isempty(codes)
            for code in codes, kw in kwsets(T)
                path = RDS.rasterpath(T, code; kw...)
                titles[path] = EcoSISTEM._filetitle(rec, layerinfo(T, code),
                                                    path)
            end
            @test allunique(values(titles))
            paths[T] = Set(keys(titles))
        end
        # No file is a layer of two rows, so none can carry two datasets' facts.
        @test all(isempty(intersect(paths[A], paths[B]))
                  for A in keys(paths), B in keys(paths) if A !== B)
    end

    # The four `rsds` files RasterDataSources cannot fetch, named as CHELSA's server does not name
    # them, which `_refusemisnamed` refuses before trying. When the names fail here,
    # RasterDataSources names them as the server does: delete that method, this testset and the
    # rows' notes in `BioClimPlus.csv`.
    @testset "CHELSA-BIOCLIM+'s rsds files are refused while their names are wrong" begin
        RDS = RasterDataSources
        T = CHELSA{BioClimPlus}
        for stat in ("max", "min", "mean", "range")
            code = Symbol("rsds_", stat)
            @test RDS.rastername(T, code) ==
                  "CHELSA_rsds_$(stat)_1981-2010_V.2.1.tif"
            isfile(RDS.rasterpath(T, code)) && continue
            @test_throws "CHELSA_rsds_1981-2010_$(stat)_V.2.1.tif" EcoSISTEM._fetchfiles(T,
                                                                                         code)
        end
        # A dataset with no such file is passed on to `getraster` as before.
        ext = Base.get_extension(EcoSISTEM, :EcoSISTEMRasterDataSourcesExt)
        @test isnothing(ext._refusemisnamed(T, (:gdd0, :swb)))
        @test isnothing(ext._refusemisnamed(CHELSA{BioClim}, (1, 12)))
    end

    @testset "moving the asset cache moves only a download path the extension chose" begin
        ext = Base.get_extension(EcoSISTEM, :EcoSISTEMRasterDataSourcesExt)
        saved = (ext._SETPATH[], get(ENV, "RASTERDATASOURCES_PATH", nothing))
        # The preference goes into the first environment on the load path that has EcoSISTEM,
        # so one made here keeps it out of every real environment.
        env = mktempdir()
        write(joinpath(env, "Project.toml"),
              "[deps]\nEcoSISTEM = \"$(Base.PkgId(EcoSISTEM).uuid)\"\n")
        pushfirst!(LOAD_PATH, env)
        try
            # As `__init__` leaves things where it chose the path.
            ENV["RASTERDATASOURCES_PATH"] = ext._SETPATH[] = EcoSISTEM.assetdir(RasterDataSources)
            cache = mktempdir()
            EcoSISTEM.set_assetdir!(cache)
            @test ENV["RASTERDATASOURCES_PATH"] ==
                  joinpath(cache, "RasterDataSources")
            @test isdir(ENV["RASTERDATASOURCES_PATH"])
            EcoSISTEM.set_assetdir!(nothing)
            @test ENV["RASTERDATASOURCES_PATH"] ==
                  EcoSISTEM.assetdir(RasterDataSources)
            # A path set since, or set before the extension loaded, stays where it was set.
            mine = mktempdir()
            ENV["RASTERDATASOURCES_PATH"] = mine
            EcoSISTEM.set_assetdir!(cache)
            @test ENV["RASTERDATASOURCES_PATH"] == mine
            ext._SETPATH[] = nothing
            EcoSISTEM.set_assetdir!(nothing)
            @test ENV["RASTERDATASOURCES_PATH"] == mine
        finally
            EcoSISTEM.set_assetdir!(nothing)
            filter!(!=(env), LOAD_PATH)
            ext._SETPATH[] = saved[1]
            if isnothing(saved[2])
                delete!(ENV, "RASTERDATASOURCES_PATH")
            else
                ENV["RASTERDATASOURCES_PATH"] = saved[2]
            end
        end
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
        # A layer the table leaves out is not part of the whole dataset: BIOCLIM+ less V2.1's 19.
        plus = EcoSISTEM._alllayercodes(CHELSA{BioClimPlus})
        @test length(plus) ==
              length(RasterDataSources.layers(CHELSA{BioClimPlus})) - 19
        @test !any(c -> occursin(r"^bio\d+$", string(c)), plus)
        @test "gdd0" in string.(plus)
        # An unmarked source cannot answer, and says so rather than returning an empty list.
        @test_throws ErrorException EcoSISTEM._alllayercodes(Int)

        # Any spelling a caller may write resolves to the one the dataset uses.
        @test EcoSISTEM._preferredcode(WorldClim{BioClim}, "bio4") == 4
        @test isnothing(EcoSISTEM._preferredcode(WorldClim{BioClim}, nothing))
    end

    @testset "the catalogue stays in the parent, keyed on any Type" begin
        # These are the parent's own, not the extension's - they use their argument
        # only to find a shipped CSV, the per-dataset corrections included.
        @test layerunit(WorldClim{BioClim}, 1) == u"°C"
        @test layerinfo(WorldClim{Climate}, :srad).unit == u"kJ*m^-2"
        @test EcoSISTEM._documentedceiling(WorldClim{BioClim}, 4) == 100.0
        @test EcoSISTEM._documentedceiling(WorldClim{BioClim}, 1) == 1.0
        @test EcoSISTEM._documentedceiling(EarthEnv{LandCover}, 1) == 1.0
        @test EcoSISTEM.datasetinfo(WorldClim{BioClim}).fetch === :getraster
        @test EcoSISTEM._stackaxis(WorldClim{Climate}) == Ti
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

    # Empty files under a temporary download directory, so nothing is fetched: what is asserted is
    # that a file already present is found whichever spelling names its layer, and recorded.
    @testset "files on disk are found by a layer's name or code, and recorded" begin
        E = EcoSISTEM
        mktempdir() do dir
            withenv("RASTERDATASOURCES_PATH" => dir) do
                for (T, name, code, kw) in ((EarthEnv{LandCover}, :cultivated_and_managed, 7, (;)),
                    (CHELSA{BioClim}, :bio1, 1, (;)),
                    (WorldClim{BioClim}, :bio1, 1, (; res = "10m")))
                    path = RasterDataSources.rasterpath(T, code; kw...)
                    mkpath(dirname(path))
                    touch(path)
                    @test E._localfiles(T, name) == [path]
                    @test E._localfiles(T, code) == [path]
                    record = only(E._readinputs(SourceSpec(T, name)))
                    @test record.path == basename(path) &&
                          record.role === :habitat
                end
                # Several codes give the files present among them, and a layer with no file yet
                # gives nothing rather than an error.
                @test E._localfiles(CHELSA{BioClim}, [:bio1, :bio3]) ==
                      [RasterDataSources.rasterpath(CHELSA{BioClim}, 1)]
                @test isempty(E._localfiles(CHELSA{BioClim}, :bio3))
                # A range names one file a value, as `getraster` maps over it, where `rasterpath`
                # given the range would name one file that cannot exist.
                for m in 1:12
                    path = RasterDataSources.rasterpath(WorldClim{Climate},
                                                        :wind,
                                                        res = "10m", month = m)
                    mkpath(dirname(path))
                    touch(path)
                end
                @test length(E._localfiles(WorldClim{Climate}, :wind,
                                           res = "10m",
                                           month = 1:12)) == 12
                @test length(E._sourcepaths(WorldClim{Climate}, :wind,
                                            res = "10m", month = 2:4)) == 3
            end
        end
    end

    @testset "dataset-typed methods are supplied here" begin
        @test !isempty(methods(EcoSISTEM.sourcecrs))
        @test !isempty(methods(EcoSISTEM.compress_landcover))
        @test !isempty(methods(EcoSISTEM.landcoverclass))
        # Looked up by name in the shipped table, never hardcoded.
        @test EcoSISTEM.landcoverclass(:open_water) isa Int
    end

    # **`read(::Type{CRUTS}, ...)` is this extension's own reader, and had no test at all.** It is
    # here rather than in the parent because it takes its variable codes and units from
    # `WorldClim{Climate}`'s shipped layer table - CRU TS has none of its own - so it needs the
    # trigger package that defines that dataset type.
    #
    # Synthetic files, so this runs on every platform and downloads nothing: what is asserted is
    # the seam, not CRU TS's own data.
    @testset "a CRU TS directory reads as a ClimateRaster carrying the source" begin
        dir = mktempdir()
        # Twelve monthly `.tif`s, which is the shape a directory read expects. Built through the
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
        cr = read(SourceSpec(CRUTS, "tavg", directory = dir))

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
