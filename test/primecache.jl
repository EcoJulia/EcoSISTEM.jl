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
# `getraster` downloads and returns paths; it does not read pixels, so this runs in a few hundred MB
# whatever it fetches. That is the whole reason a separate job can survive what the test job could not.
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
