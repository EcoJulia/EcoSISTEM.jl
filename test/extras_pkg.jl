# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Cross-validation against other packages: every `test/pkg_*.jl` checks EcoSISTEM's results against
# the package it is named for (`test/pkg_Package.jl` validates against `Package`).
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_pkg.jl"])'
#
# **An extra, not a core set, and the distinction is semantic rather than cosmetic**: the core sets
# test *this* package against itself, while this one checks it against someone else's — a different
# question, answered against a moving target, and one that only makes sense once the core sets pass.
# Being an `extras_*` also puts it after them in `runtests.jl`, which is where it belongs.
#
# **There are no `pkg_*.jl` files today**, so this set is currently a no-op. It is kept as a named
# family rather than deleted for two reasons: adding a cross-validation file needs no wiring at all
# (drop in `pkg_Foo.jl` and it runs), and **this is the set most worth being able to name**,
# because cross-validation can be *very* slow — it tends to mean extensive randomised dataset
# comparison — so running it alone, or everything else without it, matters more here than anywhere.

using Random
using Test
using EcoSISTEM
using ParallelTestRunner: find_tests, parse_args, runtests

let pkgbase = map(file -> replace(file, r"pkg_(.*).jl$" => s"\1"),
                  filter(str -> occursin(r"^pkg_.*\.jl$", str),
                         readdir(@__DIR__)))
    if length(pkgbase) > 0
        Random.seed!(1234)
        @info "Cross validation packages:"
        @testset "Cross-validation" begin
            for p in pkgbase
                println("    = $p")
            end
            println()

            runtests(EcoSISTEM, parse_args(String[]),
                     testsuite = filter(kv -> startswith(kv.first, "pkg_"),
                                        find_tests(@__DIR__)))
        end
    end
end
