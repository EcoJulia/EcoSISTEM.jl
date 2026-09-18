# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The unit tests that read real rasters: the five `test/test_*.jl` files named in `testsets.jl`,
# which `core_test.jl` leaves out.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["core_rasters.jl"])'
#
# A set of their own so that on a CI runner, where they run one at a time to keep within memory,
# they have a job of their own too.

using Random
using Test
using EcoSISTEM
using ParallelTestRunner: find_tests, runtests

# See `checkmem.jl`, and `core_test.jl` for why it is included here too.
include(joinpath(@__DIR__, "checkmem.jl"))
include(joinpath(@__DIR__, "testsets.jl"))

Random.seed!(1234)

@testset "Raster-reading unit tests" begin
    println()
    @info "Running tests for files:"
    foreach(t -> println("    = $t.jl"), RASTERTESTS)
    println()

    suite = filter(kv -> kv.first in RASTERTESTS, find_tests(@__DIR__))
    # A name in `RASTERTESTS` matching no file would silently drop that file from every set, so the
    # names are asserted, on every run rather than only on a runner.
    @test isempty(setdiff(RASTERTESTS, keys(suite)))
    # `ECOSISTEM_SERIAL_RASTERS=true` serialises them on a development machine too, which is the
    # remedy when one does hit the starvation: a laptop with a dozen live workers can report almost
    # no free memory whatever its true capacity, and the resulting error names Rasters rather than
    # the pool.
    runtests(EcoSISTEM, setargs(), testsuite = suite,
             init_worker_code = RELAXRASTERMEMCHECK,
             serial = serialonrunner(RASTERTESTS,
                                     override = "ECOSISTEM_SERIAL_RASTERS"))
end
