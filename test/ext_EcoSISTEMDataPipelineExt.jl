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
using EcoSISTEM: Track, ERA
using DataPipeline
using DataPipeline: DataRegistryHandle, ConfigDataProduct
using RasterDataSources

const EXT = Base.get_extension(EcoSISTEM, :EcoSISTEMDataPipelineExt)

@testset "the extension loads on DataPipeline and supplies each route for a run's handle" begin
    @test !isnothing(EXT)
    tracked = Track{<:DataRegistryHandle}
    for (hook, args) in ((EcoSISTEM._trackedfetch, (tracked, Any)),
        (EcoSISTEM._trackedprovide!, (tracked, Any, Any)),
        (EcoSISTEM._trackedmissing, (tracked, Any)),
        (EcoSISTEM._trackedcovers, (tracked, Any)))
        # The parent's method refuses any other handle; the extension adds the one for a run's.
        @test length(methods(hook)) == 2
        @test which(hook, Tuple{args...}).module === EXT
    end
    # The archive helper it once supplied is gone with its dependency.
    @test !isdefined(EcoSISTEM, :unziptemp)
end

@testset "every file is an output of its source's deposit, one source a run" begin
    track = Track(nothing, "era5")
    rec = EcoSISTEM.datasetinfo(ERA)
    # Each file is named as an output, which the run's `write:` entry resolves.
    name = EXT._registrationname(track, "/x/era5_t2m_1990s.nc")
    @test name isa ConfigDataProduct &&
          name.data_product == "era5/data/era5_t2m_1990s"
    # The deposit is named in its source's namespace, and a run takes one.
    run = Track(Ref(0), "era5")
    @test EXT._claimdeposit!(run, rec) == "ECMWF:era5/deposit"
    @test EXT._claimdeposit!(run, rec) == "ECMWF:era5/deposit"
    other = Track(run.handle, "worldclim")
    @test_throws "derives every output of a run from all of its inputs" EXT._claimdeposit!(other,
                                                                                           EcoSISTEM.datasetinfo(WorldClim{BioClim}))
    @test EXT._claimdeposit!(Track(Ref(0), "worldclim"),
                             EcoSISTEM.datasetinfo(WorldClim{BioClim})) ==
          "WorldClim:worldclim/deposit"
end

@testset "a tracked spec must read a catalogued source" begin
    @test EXT._sourcerow(SourceSpec(WorldClim{BioClim}, :bio1)).namespace ==
          "WorldClim"
    @test EXT._sourcerow(NaturalEarthSpec("Scotland")).dataset == "NaturalEarth"
    file = RasterFileSpec(tempname() * ".tif", axis = EcoSISTEM.NicheAxis)
    @test_throws "no catalogued source" EXT._sourcerow(file)
end

@testset "what is done once a run is remembered by its handle" begin
    run = Ref(1)
    @test !EXT._isdone(run, "citations")
    EXT._markdone!(run, "citations")
    @test EXT._isdone(run, "citations") && !EXT._isdone(Ref(1), "citations")
end

end
