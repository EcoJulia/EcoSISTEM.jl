# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Download every raster the test suite reads, and nothing else. Not a test: it makes no assertions
# and is named so that neither `runtests.jl`'s `core_*`/`extras_*` loops nor `core_test.jl`'s
# `test_*` filter picks it up.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["primecache.jl"])'
#
# It exists so that CI can fetch the data in a job of its own, and the job that runs the tests can
# restore it instead. That is not only about time. Downloading while the suite is at its peak is what
# made a GitHub runner unsurvivable: `test_datasetread.jl` alone reached 7.0 GB resident of the
# runner's 16 GB, and a concurrent ~1 GB of downloads filling the page cache on top of it took the
# machine far enough into swap that the Actions agent was shut down. The tests then never reached
# their own cache-saving step, so the next run downloaded again -- a loop that could not break itself.
#
# Most entries call `getraster`, which downloads and returns paths without reading pixels, so they
# cost a few hundred MB whatever they fetch. Two of them deliberately DO read, to precompute the
# aggregated forms the tests use: that is the expensive part, at over 10 GB resident, and it is why
# this job has to exist rather than merely being an optimisation. Nothing else runs here, so the peak
# has the machine to itself -- in the test job it did not, and the runner was shut down.
#
# An incomplete list is SAFE, and that is deliberate. `getraster`/`download` each test their own
# target file's existence before fetching, so anything missed here is simply downloaded by the test
# job as before -- the cost of a gap is a few downloads, never a wrong answer. So this list is
# maintained by hand without a mechanism to keep it honest, and should be extended when a test starts
# reading something new. The same reasoning is written up beside the cache key in
# `.github/workflows/testing.yaml`.
#
# CHELSA{BioClimPlus} is deliberately absent: those layers are 7 GB and `heavydata()` already keeps
# them off a runner entirely, so priming them would download what nothing on CI reads.

module PrimeCache

using EcoSISTEM
using RasterDataSources

# Each entry fetches one dataset, or the part of one that the suite actually reads. Kept as thunks so
# a failure in one is reported and the rest still run: a partial cache is useful and a missing one is
# not, so there is no reason for the first network hiccup to discard the others.
const WANTED = [
    "WorldClim BioClim (all layers)" => () -> getraster(WorldClim{BioClim}),
    "WorldClim Climate :wind (12 months)" =>
        () -> getraster(WorldClim{Climate},
                        :wind,
                        month = 1:12),
    "WorldClim Elevation" => () -> getraster(WorldClim{Elevation}),
    "EarthEnv LandCover (all classes)" => () -> getraster(EarthEnv{LandCover}),
    "EarthEnv HabitatHeterogeneity" =>
        () -> getraster(EarthEnv{HabitatHeterogeneity}),
    "CHELSA BioClim layer 1" => () -> getraster(CHELSA{BioClim}, 1),
    # These READ rather than download, so that each file's aggregated form is computed here and
    # cached alongside the rasters.
    #
    # A read naming a `scale` and no `cut` aggregates the WHOLE source file, memoised under
    # `assetdir()/aggregates` keyed on (file, scale, reducer, unit). Measured cold, EarthEnv
    # LandCover costs 10.6 to 12.2 GB resident -- all twelve global bands at full resolution --
    # whatever scale is asked for, since the scale only decides what comes out. That is most of a
    # 16 GB runner, and it is what killed the Linux jobs about ninety seconds into the suite.
    #
    # Paying it here fixes it everywhere: each test job restores the memoised result and loads a
    # JLD2 instead. The scales MUST match the ones the tests ask for, or the key differs and the
    # work is done again where it cannot be afforded -- `test_datasetread.jl` chooses them.
    #
    # A windowed read is deliberately absent: `_cachedlayer` bypasses the cache whenever a `cut` is
    # given and reads only the window, so there is nothing to prime.
    #
    # CHELSA is absent too, and cannot be added. Its `bio1` is a 43200 x 20880 grid whose aggregate
    # peaks at 11.7 GB read whole, 10.8 GB through the lazy fallback that exists to be the cheap
    # path, and 25 GB through the full read pipeline -- none of which fits 16 GB even in a job
    # running nothing else. This job was killed attempting exactly that. `test_datasetread.jl` skips
    # those reads on a runner instead.
    "EarthEnv LandCover, aggregated at scale 40" =>
        () -> read(EarthEnv{LandCover},
                   scale = 40),
    "WorldClim BioClim, aggregated at scale 4" =>
        () -> read(WorldClim{BioClim},
                   scale = 4),
    # Used by `examples/ScottishCultivatedLand.jl`, so it is only reached by
    # `extras_examples`; fetched here for the same reason as the rasters.
    "Scotland land-cover shapefile" =>
        () -> EcoSISTEM.assetpath(ShapeSpec("https://gis-downloads.nature.scot/LSCMAP_SCOTLAND_SHP_27700.zip").path)
]

failed = String[]
for (name, fetch) in WANTED
    print("    * ", name, " ... ")
    try
        fetch()
        println("ok")
    catch e
        println("FAILED")
        @error "Could not prime $name" exception = (e, catch_backtrace())
        push!(failed, name)
    end
end

if isempty(failed)
    @info "Raster cache primed: all $(length(WANTED)) entries present."
else
    # Reported rather than thrown. A partly primed cache is still worth saving, and the test job
    # fetches whatever is missing exactly as it did before this job existed.
    @warn "Raster cache primed with gaps; the test job will fetch these itself." failed
end

end
