# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Run every `test/clean_*.jl` repo-hygiene check - the JuliaFormatter and ResearchSoftwareMetadata
# crosswalk gates.
#
# On its own:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_clean.jl"])'
#
# **Not** runnable as a bare script, and that is a dependency fact rather than an oversight: the
# checks need `Git`, `JuliaFormatter` and `ResearchSoftwareMetadata`, which are `[extras]` in
# `Project.toml`'s `test` target, plus package deps like `Phylo` as *direct* deps - a combination only
# `Pkg.test` provisions. Going through it means the environment is right by construction rather than
# reconstructed here and drifting the moment a dependency changes.
#
# `runtests.jl` also picks it up automatically, along with every other `test/extras_*.jl` - but only
# after the unit tests have passed, since a failing testset throws before the loop is reached. There
# it *skips* on a CI runner, because the dedicated `metadata.yaml` job asks for this file by name and
# there is no reason for every matrix entry to redo the crosswalk.
#
# These checks read the *working tree*, so they only mean anything on a tree whose changes are
# staged or committed: `is_repo_clean` counts unstaged tracked changes, and the crosswalk regenerates
# `codemeta.json`. A dirty tree fails them for reasons that have nothing to do with the code.

module ExtrasClean

using Test
using EcoSISTEM

# Run when this file was asked for by name (`test_args = ["extras_clean.jl"]`, which is what the
# `metadata.yaml` job passes), and when running locally. Skip on a CI runner that merely reached this
# as part of the whole suite - no environment variable needed to tell those apart, because "was I
# asked for?" is the question that actually distinguishes them.
const ASKED_FOR = any(a -> occursin("extras_clean", a), ARGS)

if !(ASKED_FOR || !haskey(ENV, "RUNNER_OS"))
    @info "Skipping the hygiene checks: this is a CI runner and they were not asked for directly, " *
          "so the dedicated metadata job runs them instead."
else
    # Each `clean_Foo.jl` checks one kind of repo hygiene, named after what it checks.
    cleanbase = map(file -> replace(file, r"clean_(.*).jl$" => s"\1"),
                    filter(str -> occursin(r"^clean_.*\.jl$", str),
                           readdir(@__DIR__)))

    if length(cleanbase) > 0
        println()
        @info "Crosswalk and clean testing:"
        for c in cleanbase
            println("    = $c")
        end
        println()

        @testset for c in cleanbase
            # Relative to *this file*, which `include` resolves against - so the caller's working
            # directory does not matter.
            fn = joinpath(@__DIR__, "clean_$c.jl")
            println("    * Verifying $c.jl ...")
            include(fn)
        end
    end
end

end
