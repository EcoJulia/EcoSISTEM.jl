# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Run the canonical suite — the blessed-result tests in `test/canonical/`.
#
# On its own:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_canonical.jl"])'
#
# To re-bless after an intended change (read `test/canonical/README.md` first):
#
#     ECOSISTEM_BLESS=true julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_canonical.jl"])'
#
# Kept out of the unit suite because these run whole simulations and read real rasters. The
# directory was, until now, run by *nothing at all*: `runtests.jl` discovers `test_*.jl` directly
# inside `test/`, so a `canonical/` subdirectory matches nothing and a file there is never run —
# which is how one can sit dead long enough to accumulate assertions that no longer pass. This file
# is what makes the directory reachable.
#
# Discovery here is by glob too (`canonical/test_*.jl`), which is the same trap: a file dropped
# into `canonical/` under any other name is not run and will rot exactly as that one did.

module ExtrasCanonical

using Test
using EcoSISTEM
using ParallelTestRunner: find_tests, parse_args, runtests

include(joinpath(@__DIR__, "canonical", "canonical.jl"))
using .Canonical

# See `checkmem.jl`. Included here as well as in `runtests.jl` because this file is documented as a
# standalone script, and passed to the workers below because the guard is per-process.
include(joinpath(@__DIR__, "checkmem.jl"))

@testset "Canonical results" begin
    dir = joinpath(@__DIR__, "canonical")
    files = sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"),
                        readdir(dir)))
    println()
    @info "Running canonical tests" * (blessing() ? " — RE-BLESSING" : "")
    if blessing()
        # **Blessing stays SERIAL, deliberately.** Each worker is a separate process with its own
        # `RECORDED`, so the parent's `writereference()` would see an empty one — and worse,
        # `_writethrough` is a read-modify-write of `reference.toml` per value, so several workers
        # blessing at once would race and silently drop results. Blessing is rare and interactive;
        # checking is what runs on every build, and that is what is worth parallelising.
        for f in files
            println("    * ", f, " ...")
            include(joinpath(dir, f))
        end
        # Written once, after every file, and merged rather than replaced — see `writereference`.
        writereference()
    else
        # `find_tests` on the **canonical subdirectory**, so the keys are bare file names; the
        # `test_` filter then excludes `canonical.jl` itself, which is a helper rather than a test.
        # `parse_args(String[])` — see `core_test.jl` for why forwarding `ARGS` runs nothing.
        runtests(EcoSISTEM, parse_args(String[]),
                 testsuite = filter(kv -> startswith(kv.first, "test_"),
                                    find_tests(dir)),
                 init_worker_code = RELAXRASTERMEMCHECK)
    end
end

end
