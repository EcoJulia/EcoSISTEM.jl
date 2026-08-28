# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The unit tests: every `test/test_*.jl`, which test the matching `src/*.jl`.
#
# Run this set on its own — the longest chunk of the suite, and the one worth isolating:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["core_test.jl"])'
#
# That is the point of this file existing. Any *single* file could already be named
# (`test_args = ["test_Layer.jl"]`), but there was no way to say "all the unit tests and nothing
# else" short of naming all of them — so the only alternative to one file was the whole suite,
# extras included.

using Random
using Test
using EcoSISTEM
using ParallelTestRunner: find_tests, parse_args, runtests

# See `checkmem.jl`: Rasters' pre-read memory guard misreports free memory on a macOS runner. Included
# here as well as in `runtests.jl` because this file is documented as runnable on its own, and handed
# to the workers below because the setting is per-process.
include(joinpath(@__DIR__, "checkmem.jl"))

# The raster-reading files, run one at a time rather than alongside each other. Measured peak RSS on
# a GitHub runner: test_rasters 4140 MB, test_datasetread 4127 MB, test_StudyArea 3933 MB,
# test_deprecations 3753 MB, against a few hundred for most of the rest. ParallelTestRunner sizes its
# worker pool at one per 2 GiB of available memory, which is a fair budget for an ordinary test file
# and roughly half what any of these needs, so four of them in flight at once exhausted a 16 GB Linux
# runner and the job was killed outright.
#
# The set is those that read real rasters, not those that merely allocate: what makes them expensive
# is holding a downloaded layer, so a file joins this list when it starts naming a dataset, not when
# it gets slower. Regenerate the candidates with
# `grep -l "WorldClim{\|EarthEnv{\|CHELSA{" test/*.jl`, which over-reports -- most of those touch
# only a fixture -- and keep the ones that actually materialise a layer.
#
# `:before` rather than the default position, so that test_datasetread performs the downloads while
# nothing else is running. Several workers fetching into one scratchspace at the same time is the
# race that ordering, here and in `runtests.jl`, exists to avoid.
const SERIALTESTS = ["test_datasetread", "test_rasters", "test_StudyArea",
                     "test_deprecations", "test_GridHabitat"]

# Identify files in test/ that are testing matching files in src/
#  - src/Source.jl will be matched by test/test_Source.jl
let filebase = String[]
    for (root, dirs, files) in walkdir(joinpath(@__DIR__, "..", "src"))
        append!(filebase,
                map(file -> replace(file, r"(.*).jl" => s"\1"),
                    filter(file -> occursin(r".*\.jl", file), files)))
    end

    testbase = map(file -> replace(file, r"test_(.*).jl" => s"\1"),
                   filter(str -> occursin(r"^test_.*\.jl$", str),
                          readdir(@__DIR__)))

    # Identify tests with no matching file
    superfluous = filter(f -> f ∉ filebase, testbase)
    if length(superfluous) > 0
        println()
        @info "Potentially superfluous tests:"
        for f in superfluous
            println("    + $f.jl")
        end
        println()
    end

    # Identify files with no matching test
    notest = filter(f -> f ∉ testbase, filebase)
    if length(notest) > 0
        println()
        @info "Potentially missing tests:"
        for f in notest
            println("    - $f.jl")
        end
        println()
    end

    # Seeded here rather than only in `runtests.jl`, or this file run on its own would not be
    # reproducible — which is most of why it exists.
    Random.seed!(1234)

    @testset "Unit tests" begin
        @test isfile(pkgdir(EcoSISTEM, "test", "runtests.jl"))
        println()
        @info "Running tests for files:"
        for t in testbase
            println("    = $t.jl")
        end
        println()

        suite = filter(kv -> startswith(kv.first, "test_"), find_tests(@__DIR__))

        # ParallelTestRunner silently drops a `serial` entry that names nothing in the suite, so a
        # renamed file would stop being serialised with no signal at all -- and the symptom would be
        # an out-of-memory kill on a runner, days later and nowhere near the rename. Assert the names
        # instead, so the failure is here and says which one went.
        @test isempty(setdiff(SERIALTESTS, keys(suite)))

        runtests(EcoSISTEM, parse_args(String[]), testsuite = suite,
                 init_worker_code = RELAXRASTERMEMCHECK,
                 serial = SERIALTESTS, serial_position = :before)
    end
end
