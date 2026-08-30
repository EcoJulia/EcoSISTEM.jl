# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The extension tests: every `test/ext_*.jl`, which test the matching `ext/*.jl`.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["core_ext.jl"])'
#
# MPI is a weak dependency but *is* installed in the test environment, so `EcoSISTEMMPIExt` loads
# and `ext_EcoSISTEMMPIExt.jl` launches `SmallMPItest.jl` under `mpiexec` at 1, 2 and 4 ranks — the
# slowest thing in this set by some way.

using Random
using Test
using EcoSISTEM
# `[C7-VIS]` B1/B2/B3: these are `public` rather than exported, so they must be named.
using EcoSISTEM: empty_mpi_gridlandscape, synchronise_from_cols!,
                 synchronise_from_rows!
using ParallelTestRunner: find_tests, parse_args, runtests

# **Where an extension's docstrings are allowed to live.**
#
# `docs/src/api.md` is an `@autodocs` block over the package's own modules, so a docstring inside an
# extension is invisible to the manual — which is why every public name whose implementation moves to
# one leaves its docstring behind on a method-less stub in the parent.
#
# **The failure this catches is silent.** Copying a docstring to the stub instead of *moving* it
# leaves the name documented twice; Documenter renders two entries with the same HTML anchor and
# still exits 0. That happened to `sourcecrs` in 3e and was found by reading the file, not by a gate.
#
# The exception list is by **file and count**, not a total, so a new docstring has to be looked at
# and a removed one shrinks the list rather than passing silently.
@testset "extension docstrings" begin
    # `Base.read`'s dataset methods are the one thing that cannot take a parent stub: the parent
    # would have to define `read(::Type, …)`, pirating `Base.read` for every type in Julia. They are
    # reached instead by naming their module in `api.md`'s `Modules`.
    allowed = Dict("EcoSISTEMRasterDataSourcesExt/read.jl" => 3,
                   # **Pre-existing and NOT endorsed** — `EcoSISTEMMPIExt` has carried these since
                   # before the rule existed, and moving its files out of `src/` in 3f is what first
                   # put them under this guard. Measured: five names (`MPIEcosystem`,
                   # `MPIGridLandscape`, `empty_mpi_gridlandscape`, `synchronise_from_cols!`,
                   # `synchronise_from_rows!`) are documented in **both** here and `src/EcoSISTEM.jl`,
                   # and at least one (`rows_matrix`) only here, so it never reaches the manual at
                   # all. Harmless *today* only because this module is not named in `api.md`'s
                   # `Modules`, so the duplicates are invisible rather than rendered twice.
                   # Listed with real counts, not exempted: a **new** docstring in these files
                   # still fails this test and has to be looked at.
                   "EcoSISTEMMPIExt/Ecosystem.jl" => 2,
                   "EcoSISTEMMPIExt/dynamics.jl" => 7,
                   "EcoSISTEMMPIExt/Landscape.jl" => 4)
    extdir = joinpath(@__DIR__, "..", "ext")
    for (root, _, files) in walkdir(extdir)
        for file in filter(endswith(".jl"), files)
            path = joinpath(root, file)
            rel = replace(relpath(path, extdir), '\\' => '/')
            src = read(path, String)
            # `#= … =#` blocks are parked code, not documentation
            src = replace(src, r"#=(?:.|\n)*?=#" => "")
            # The `m` flag is load-bearing: without it `^` anchors to the start of the whole
            # string, every count is zero, and the guard silently passes anything.
            n = count(!isnothing, eachmatch(r"^\"\"\""m, src)) ÷ 2
            @test n == get(allowed, rel, 0)
        end
    end
end

# Identify files in test/ that are testing matching files in ext/
#  - ext/SourceExt.jl will be matched by test/ext_SourceExt.jl
let filebase = String[]
    # **Extension *modules*, not every `.jl` under `ext/`.** An extension may be laid out as a
    # single `ext/Name.jl` or as a directory `ext/Name/Name.jl` including further files (see
    # `EcoSISTEMRasterDataSourcesExt`), and only the module is a thing a test can name — counting its
    # parts would report each one as an untested "extension".
    let extdir = joinpath(@__DIR__, "..", "ext")
        for entry in readdir(extdir)
            if endswith(entry, ".jl") && isfile(joinpath(extdir, entry))
                push!(filebase, replace(entry, r"(.*).jl" => s"\1"))
            elseif isdir(joinpath(extdir, entry)) &&
                   isfile(joinpath(extdir, entry, "$entry.jl"))
                push!(filebase, entry)
            end
        end
    end

    extbase = map(file -> replace(file, r"ext_(.*).jl" => s"\1"),
                  filter(str -> occursin(r"^ext_.*\.jl$", str),
                         readdir(@__DIR__)))

    # Identify tests with no matching file
    superfluous = filter(f -> f ∉ filebase, extbase)
    if length(superfluous) > 0
        println()
        @info "Potentially superfluous extension tests:"
        for f in superfluous
            println("    + $f.jl")
        end
        println()
    end

    # Identify files with no matching test
    notest = filter(f -> f ∉ extbase, filebase)
    if length(notest) > 0
        println()
        @info "Potentially missing extension tests:"
        for f in notest
            println("    - $f.jl")
        end
        println()
    end

    Random.seed!(1234)

    @testset "Extension tests" begin
        println()
        @info "Running tests for extensions:"
        for t in extbase
            println("    = $t.jl")
        end
        println()

        # `parse_args(String[])`, **not** `parse_args(ARGS)`: `ARGS` is the `test_args` that
        # selected *this file*, and `filter_tests!` keeps a test only if its name starts with one of
        # them — so forwarding them filters the suite to empty and reports success over zero tests.
        runtests(EcoSISTEM, parse_args(String[]),
                 testsuite = filter(kv -> startswith(kv.first, "ext_"),
                                    find_tests(@__DIR__)))
    end
end
