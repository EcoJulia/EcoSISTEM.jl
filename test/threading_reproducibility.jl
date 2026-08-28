# SPDX-License-Identifier: LGPL-3.0-or-later

## Run a seeded simulation and save the final abundance matrix to the path given
## in ARGS[1]. This script is launched as a subprocess with different thread
## counts by the "Multithreaded reproducibility" testset in test_dynamics.jl, so
## that the parallel `update!` path is genuinely exercised with more than one
## thread. Because the ecosystem is seeded (one RNG stream per species), the
## result must be identical regardless of the number of threads used.

using EcoSISTEM
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using JLD2

include(pkgdir(EcoSISTEM, "test", "TestCases.jl"))

# **Checking that an argument EXISTS is not enough, and this file destroyed itself proving it
# (2026-08-13).** Named in `Pkg.test(test_args = ["threading_reproducibility.jl"])` — which is how
# every *real* test file here is run — `ARGS` becomes that list, so `ARGS[1]` was this script's own
# name. The length check passed and `@save ARGS[1]` wrote a 136 KB JLD2 over the source.
#
# So the guard is on the argument's **shape**, not its presence: an output path must end `.jld2`,
# which the harness in `test_dynamics.jl` supplies and a stray `test_args` never will.
length(ARGS) == 1 && endswith(ARGS[1], ".jld2") ||
    error("threading_reproducibility.jl is NOT a test file — it is launched as a subprocess by " *
          "the \"Multithreaded reproducibility\" testset in test_dynamics.jl, and writes its " *
          "result to the `.jld2` path given as its only argument.\n" *
          "usage: julia threading_reproducibility.jl <output.jld2>\n" *
          "got ARGS = $(ARGS)")

# The seed makes both the initial abundances and the per-species simulation RNGs
# deterministic, so the whole run is reproducible.
eco = Test1Ecosystem(seed = 1234)
for _ in 1:50
    EcoSISTEM.update!(eco, 1month_mean_duration)
end

matrix = copy(eco.abundances.matrix)
@save ARGS[1] matrix

println("Saved reproducibility result with $(Threads.nthreads()) threads to $(ARGS[1])")
