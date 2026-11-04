# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `ext/EcoSISTEMDataPipelineExt/`, a tracked spec's files through a FAIR data pipeline run.
# Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["ext_EcoSISTEMDataPipelineExt.jl"])'
# or as part of `core_ext.jl`.
#
# A run's handle exists only inside a `fair run` with a registry behind it, so the calls themselves
# - `link_register!`, `link_write!`, `link_read!`, `identify`, `isregistered` - are exercised by hand,
# by `examples/pipeline/`. What each file is registered as is decided by the plain functions of
# `src/registration.jl`, which `test_registration.jl` tests; this file checks the extension supplies
# the routes and the parts of them that need no registry.

module TestEcoSISTEMDataPipelineExt

using Test
using EcoSISTEM
using EcoSISTEM: Track, ERA, Provenancer, SourceFacts, DepositIsFile
using DataPipeline
using DataPipeline: DataRegistryHandle, ConfigDataProduct
using RasterDataSources

const EXT = Base.get_extension(EcoSISTEM, :EcoSISTEMDataPipelineExt)

@testset "the extension loads on DataPipeline and supplies each route for a run's handle" begin
    @test !isnothing(EXT)
    tracked = Track{<:DataRegistryHandle}
    for (hook, args) in ((EcoSISTEM._trackedfetch, (tracked, Any)),
        (EcoSISTEM._trackedprovide!, (tracked, Any, Any)),
        (EcoSISTEM._trackedmissing, (tracked, Any)))
        # The parent's method refuses any other handle; the extension adds the one for a run's.
        @test length(methods(hook)) == 2
        @test which(hook, Tuple{args...}).module === EXT
    end
    # The archive helper it once supplied is gone with its dependency.
    @test !isdefined(EcoSISTEM, :unziptemp)
end

@testset "every file is an output of its source's deposit, one source a run" begin
    track = Track(nothing, "era5")
    # Each file is named as an output, which the run's `write:` entry resolves.
    name = EXT._registrationname(track, "/x/era5_t2m_1990s.nc")
    @test name isa ConfigDataProduct &&
          name.data_product == "era5/data/era5_t2m_1990s"
    # The deposit is named in its source's namespace, and a run takes one source.
    run = Track(Ref(0), "era5")
    era5 = SourceFacts(source = "ERA5", doi = "10.24381/cds.f17050d7",
                       namespace = "ECMWF")
    @test EXT._claimsource!(run, era5, "/x/a.nc") == "ECMWF:era5/deposit"
    @test EXT._claimsource!(run, era5, "/x/b.nc") == "ECMWF:era5/deposit"
    other = Track(run.handle, "worldclim")
    worldclim = SourceFacts(source = "WorldClim", doi = "10.1002/joc.5086",
                            namespace = "WorldClim")
    @test_throws "derives every output of a run from all of its inputs" EXT._claimsource!(other,
                                                                                          worldclim,
                                                                                          "/x/c.tif")
    @test EXT._claimsource!(Track(Ref(0), "worldclim"), worldclim,
                            "/x/c.tif") ==
          "WorldClim:worldclim/deposit"
    # A file registered as its source itself is its run's one source, named as the file, which no
    # `write:` entry resolves; a second file, or a deposit after it, is refused.
    given = SourceFacts(source = "Given", key = "given", identifies = :data,
                        namespace = "Mine", deposit = DepositIsFile())
    one = Track(Ref(0), "given")
    @test EXT._claimsource!(one, given, "/x/table.csv") ==
          "Mine:given/data/table"
    @test EXT._claimsource!(one, given, "/x/table.csv") ==
          "Mine:given/data/table"
    @test_throws "derives every output" EXT._claimsource!(one, given,
                                                          "/x/other.csv")
    @test_throws "derives every output" EXT._claimsource!(Track(one.handle,
                                                                "era5"),
                                                          era5, "/x/a.nc")
    @test EXT._registrationname(one, "/x/table.csv", given) ==
          "Mine:given/data/table"
    @test EXT._registrationname(run, "/x/a.nc", era5) ==
          EXT._registrationname(run, "/x/a.nc")
    # A source with no namespace of its own is registered in the run's default.
    @test EXT._depositof(Track(Ref(0), "box"), nothing) == "box/deposit"
end

@testset "a tracked spec's source is its catalogue row's" begin
    @test EcoSISTEM._trackedsource(SourceSpec(WorldClim{BioClim}, :bio1)).facts.namespace ==
          "WorldClim"
    ne = EcoSISTEM._trackedsource(NaturalEarthSpec("Scotland"))
    @test ne.facts.key == "naturalearth" && ne.fileversioned
    file = RasterFileSpec(tempname() * ".tif", axis = EcoSISTEM.NicheAxis)
    @test_throws "no catalogued source" EcoSISTEM._trackedsource(file)
    # The routes for a file of the user's own, and for a download, are the extension's too.
    tracked = Track{<:DataRegistryHandle}
    for (hook, args) in ((Provenancer._trackedrecord, (tracked, Any, Any)),
        (Provenancer._trackedinputmissing,
        (tracked, AbstractString)),
        (Provenancer._trackedfetchfile,
        (tracked, Provenancer.URLFile)))
        @test length(methods(hook)) == 2
        @test which(hook, Tuple{args...}).module === EXT
    end
end

@testset "what is done once a run is remembered by its handle" begin
    run = Ref(1)
    @test !EXT._isdone(run, "citations")
    EXT._markdone!(run, "citations")
    @test EXT._isdone(run, "citations") && !EXT._isdone(Ref(1), "citations")
end

end
