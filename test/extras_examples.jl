# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Run every top-level `examples/*.jl` script as an integration test, in the `examples/` environment.
#
# On its own, either through the test harness (which provisions the test environment):
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_examples.jl"])'
#
# or as a bare script, which works here because this file activates its own environment:
#
#     cd test && julia --project=.. extras_examples.jl
#
# `runtests.jl` also picks it up automatically, along with every other `test/extras_*.jl` — but only
# after the unit tests have passed, since a failing testset throws before the loop is reached.

module ExtrasExamples

# Headless plots — see the fuller note in `runtests.jl`. Repeated here because this file is
# documented as a standalone script, and running it *is* running tests. `get!` is idempotent, so
# setting it in both places is harmless.
get!(ENV, "GKSwstype", "100")

using Test
using EcoSISTEM
using Pkg

let
    # Restore whatever environment we were called with. This file activates `examples/`, and
    # nothing here activated it back before — which is a real bug, not just untidiness: `Git`,
    # `JuliaFormatter` and `ResearchSoftwareMetadata` are in the **test** environment only, so
    # `extras_clean.jl` would fail to load them if it ran after this with `examples/` still active.
    # It has never fired only because the suite has never got this far with the unit tests green.
    original_project = dirname(Base.active_project())
    try
        @testset "Examples folder" begin
            println()
            @info "Running from examples folder ..."
            Pkg.activate(pkgdir(EcoSISTEM, "examples"))
            Pkg.resolve()

            # Load these HERE, before the loop, purely to move one piece of stderr output
            # outside `@test_nowarn` — which fails on *any* stderr output, not just warnings.
            #
            # An extension is built lazily on the first `using` of its trigger, and EcoSISTEM is
            # already loaded from a different environment, whose build ID is absent from the
            # examples environment's cache. Julia therefore logs "Module EcoSISTEM ... is missing
            # from the cache" while building it, and the first example to trigger one failed on
            # that alone — on every cold cache (a fresh CI runner), while passing locally on a warm
            # one, which is what made it look intermittent.
            #
            # **`Plots` is not an extension trigger**: the plot recipes live in the package proper,
            # since a `@recipe` needs only `RecipesBase`, a hard dependency. The trigger that matters
            # here is `RasterDataSources`, a weak dependency the examples use heavily. Both are
            # imported because `Plots` pulls a large stack in, and doing it here rather
            # than inside the first example keeps that noise out of `@test_nowarn` too.
            #
            # `Pkg.precompile()` does **not** fix this — measured, both ways: with it alone the
            # build still happens at `using` time and the test still fails; with only these imports
            # and no `Pkg.precompile()` it passes. So the imports are the whole fix, and they keep
            # `@test_nowarn`'s real job — catching warnings the examples themselves emit — intact.
            @eval using RasterDataSources
            @eval using Plots
            # **Every `.jl` at the top level of `examples/` is run** — there is no `test_` prefix
            # to opt in with, which is the point: an example that is not exercised rots, and the
            # directory says so in its `README.md`. Anything not meant to run every time lives in a
            # subdirectory (`other/`, `HPC/`, `scripts/`, …), which this deliberately does not walk.
            exampledir = pkgdir(EcoSISTEM, "examples")
            examplefiles = sort(filter(str -> occursin(r"\.jl$", str),
                                       readdir(exampledir)))
            # Timed per file and reported at the end, cheapest to keep than to re-derive: running
            # one example on its own pays full compilation, so a standalone timing says little about
            # its share of a whole run. This is what says where the gate's minutes actually go.
            timings = Pair{String, Float64}[]
            for fn in examplefiles
                # An absolute path, so this works whatever the caller's working directory is —
                # `include` resolves a relative path against *this* file, not the cwd.
                println("    * Testing $fn ...")
                # **Still `@test_nowarn`, and the reason is worth recording.** It rejects *any*
                # stderr, which briefly made the package's own `DefaultEcosystem()` builders
                # unusable here — they deliberately **announce** each default they fill in. The fix
                # went into the *package* rather than this harness: those builders now take
                # `verbosity = :silent`, matching `StudyArea`. So an example can use them and stay
                # silent, and this keeps its stricter and simpler check.
                elapsed = @elapsed @test_nowarn include(joinpath(exampledir,
                                                                 fn))
                push!(timings, fn => elapsed)
            end

            # Slowest first: the only ordering anyone reads this for is "what should I shrink?".
            println("\n    Time by example (total $(round(sum(last, timings), digits = 1)) s):")
            for (fn, elapsed) in sort(timings, by = last, rev = true)
                println("      ", rpad(fn, 30),
                        lpad(round(elapsed, digits = 1), 7),
                        " s")
            end
        end
    finally
        Pkg.activate(original_project)
    end
end

end
