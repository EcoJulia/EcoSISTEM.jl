# SPDX-License-Identifier: LGPL-3.0-or-later

using Random
using Test
using EcoSISTEM
using Pkg
using ParallelTestRunner: find_tests, parse_args, runtests

# Rasters' pre-read memory guard misreports on a CI runner; see the file for the measurement and for
# what relaxing it gives up. Included here so the main process is covered, and re-run inside each
# ParallelTestRunner worker through `init_worker_code` at every `runtests` call in the suite.
include("checkmem.jl")

# A test argument names one test file to run *instead of* the whole suite:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_clean.jl"])'
#
# Going through `Pkg.test` rather than running the file directly is the whole point: it is what
# provisions the test environment. `Git`, `JuliaFormatter` and `ResearchSoftwareMetadata` are
# `[extras]` in `Project.toml`'s `test` target and nothing else supplies them, so a bare
# `julia test/extras_clean.jl` dies on `using Git` - while a package dependency like `Phylo` (which
# `clean_ResearchSoftwareMetadata.jl` uses) needs to be a *direct* dep of the active project too.
# `Pkg.test` gets both right by construction; anything else reconstructs it and drifts.
#
# The `.jl` is optional, so `test_args = ["extras_clean"]` works too. Any test file may be named -
# `test_Layer.jl` as readily as a whole set.
#
# **The suite is eight nameable sets**, which is what lets a full run go in parallel rather than in
# sequence:
#
#     core_test  core_ext
#     extras_canonical  extras_clean  extras_docs  extras_examples  extras_notebooks  extras_pkg
#
# The split is semantic: the **core** sets test this package against itself, the **extras** check it
# against something outside - the examples, the notebooks, the repo's own hygiene, the blessed
# results, another package.
#
# `core_test` is by far the longest (32 files, ~10 min) and most others ~1-2 min each, so running them
# concurrently takes about as long as `core_test` alone rather than the sum. `extras_pkg` is the
# exception waiting to happen: cross-validation against another package tends to mean extensive
# randomised dataset comparison and can be *very* slow, which is exactly why it is separately
# nameable.
#
# **Running them in parallel gives up the ordering guarantee below** - the extras then run even
# when the unit tests are failing, so one broken thing reports as several. That is a fair trade when
# chosen deliberately (the ten minutes are worth more than the duplicate failures), and it is why the
# sequential order stays the default here rather than the other way round.
#
# **One practical tip if you do parallelise**: let the first invocation get through precompilation
# before starting the rest, or every process compiles the same package at once and they contend.

# **The examples run small under the suite and full-scale outside it**, and this line is the whole
# mechanism: reaching this file *is* the signal, so there is nothing to detect. Examples read
# `ECOSISTEM_SCALE`, and a direct `julia --project=examples examples/landscapes.jl` leaves it unset,
# which they take as full scale.
#
# Two examples skip themselves **entirely** when it is `small`, rather than merely shrinking:
# `examples/biodiversity.jl` (~20 whole simulations) and `examples/models.jl` (50-year runs that draw
# plots). Both are far too slow for any routine gate, so a normal `Pkg.test` never executes them.
#
# `get!` rather than assignment, so an explicit setting from the caller always wins:
# `ECOSISTEM_SCALE=large julia --project -e 'using Pkg; Pkg.test()'` runs the big configurations
# under the suite. Subprocesses inherit `ENV`, so this reaches spawned workers too.
get!(ENV, "ECOSISTEM_SCALE", "small")

# **Plots stay headless under the suite, and pop up when you run an example by hand** - the same
# mechanism as `ECOSISTEM_SCALE` above, for the same reason: reaching this file *is* the signal.
# `GKSwstype = "100"` is GR's null workstation, so a figure is built and discarded instead of opening
# a window. **A CI runner has no display**, and the suite genuinely draws: two top-level examples,
# two Pluto notebooks, and (at full scale) `examples/models.jl`.
#
# `get!` again, so a caller who has set `GKSwstype` keeps it - and a direct
# `julia --project=examples examples/AvailableGround.jl` leaves it unset, which is what makes the
# window appear when a person is watching.
get!(ENV, "GKSwstype", "100")

# **Unit rendering is made the same on every platform, because assertions depend on it.** Unitful
# decides between `m s^-1` and `m s^-1` with
# `get(ENV, "UNITFUL_FANCY_EXPONENTS", Sys.isapple() ? "true" : "false")` -- so the fancy form is the
# default on macOS and the plain one on Linux and Windows. Three tests assert an error message or a
# plot title containing a unit, and so passed on macOS and failed everywhere else; measured on a
# Windows runner as `"WindSpeed (m s^-1)" == "WindSpeed (m s⁻¹)"`.
#
# Set rather than defaulted with `get!`, unlike the two above, and the difference is the point: those
# are preferences a caller may legitimately override, while the assertions here are only meaningful
# for one of the two renderings. A caller who set it to `false` would break the suite, so the suite
# decides.
ENV["UNITFUL_FANCY_EXPONENTS"] = "true"

requested = map(a -> endswith(a, ".jl") ? a : a * ".jl", ARGS)
for fn in requested
    isfile(joinpath(@__DIR__, fn)) ||
        error("`$fn` was asked for with `test_args`, but there is no such file in `test/`.")
end

if !isempty(requested)
    @info "Running only the requested test file(s): " * join(requested, ", ")
    @testset for fn in requested
        println("    * Running $fn ...")
        include(fn)
    end
else
    # Two loops, and nothing else. Each `core_*.jl` and `extras_*.jl` is a standalone set that can be
    # run on its own by name (see above); this file only decides the order they go in.
    #
    # The extras run **after** the core sets deliberately: a failing `@testset` throws at its end,
    # so the extras are reached only once the unit and extension tests pass. There is no point
    # running the examples, the notebooks or the hygiene checks against a broken package.
    corebase = sort(filter(str -> occursin(r"^core_.*\.jl$", str), readdir()),
                    order = Base.Order.Reverse)
    extrabase = sort(filter(str -> occursin(r"^extras_.*\.jl$", str),
                            readdir()))

    println()
    @info "Running the core test sets:"
    foreach(f -> println("    = $f"), corebase)
    println()

    @testset "EcoSISTEM.jl" begin
        @testset for fn in corebase
            println("    * Running $fn ...")
            include(fn)
        end
    end

    # **The extras do not run on a Windows CI runner.** They check this package against things
    # outside it - the examples, the notebooks, the repo's hygiene, the blessed results - and none of
    # that is OS-specific in a way the other runners do not already cover, while Windows is the
    # slowest of the three and was pushing the job past its timeout. `extras_notebooks.jl` already
    # skipped itself there; this generalises that to the whole group rather than to one file.
    #
    # Deliberately keyed on `RUNNER_OS`, not on `Sys.iswindows()`: a developer running the suite
    # on Windows locally still gets the extras, because the reason to skip is the CI budget rather
    # than anything being broken.
    skipextras = get(ENV, "RUNNER_OS", "") == "Windows"
    if skipextras
        println()
        @info "Skipping the extra test suites on a Windows runner."
    elseif !isempty(extrabase)
        println()
        @info "Running the extra test suites:"
        foreach(f -> println("    = $f"), extrabase)
        println()

        # **Three of the extras run CONCURRENTLY; the tests inside each stay serial.**
        #
        # **Why these three and not all of them.** `core_test`, `core_ext` and `extras_canonical`
        # are already parallel *internally*, so putting them in here too would nest parallelism and
        # oversubscribe the machine. These three are the ones still serial inside, which is exactly
        # what makes them the ones worth running side by side.
        #
        # **Why set-level rather than per-example.** Each of these activates a *different*
        # environment (`examples/`, `notebooks/`), and one process per set isolates that for free -
        # no `Pkg.activate` restore dance, and no two examples racing to fetch the same raster into
        # one cache. The one thing each worker is told is `checkmem.jl`, which is a property of the
        # runner rather than of the set.
        #
        # **The download race is handled by the ORDER, not by luck**: `core_test` runs first and
        # `test_GridHabitat.jl`/`test_datasetread.jl` already fetch the big EarthEnv layer the examples
        # use, so the cache is warm before any of this starts. Moving the extras ahead of the core
        # sets would reintroduce the race.
        parallelextras = ["extras_docs", "extras_examples", "extras_notebooks"]
        serialextras = filter(fn -> chop(fn, tail = 3) ∉ parallelextras,
                              extrabase)

        # The rest first, in order, exactly as before.
        @testset for fn in serialextras
            println("    * Running $fn ...")
            include(fn)
        end

        suite = filter(kv -> kv.first in parallelextras, find_tests(@__DIR__))
        if !isempty(suite)
            println()
            @info "Running these concurrently: " *
                  join(sort(collect(keys(suite))), ", ")
            println()
            runtests(EcoSISTEM, parse_args(String[]), testsuite = suite,
                     init_worker_code = RELAXRASTERMEMCHECK)
        end
    end
end
