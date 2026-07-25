# SPDX-License-Identifier: LGPL-3.0-or-later

## Run a seeded simulation and save the final abundance matrix to the path given
## in ARGS[1]. This script is launched as a subprocess with different thread
## counts by the "Multithreaded reproducibility" testset in test_Generate.jl, so
## that the parallel `update!` path is genuinely exercised with more than one
## thread. Because the ecosystem is seeded (one RNG stream per species), the
## result must be identical regardless of the number of threads used.

using EcoSISTEM
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using JLD2

include(pkgdir(EcoSISTEM, "test", "TestCases.jl"))

length(ARGS) >= 1 || error("usage: threading_reproducibility.jl <output.jld2>")

# The seed makes both the initial abundances and the per-species simulation RNGs
# deterministic, so the whole run is reproducible.
eco = Test1Ecosystem(seed = 1234)
for _ in 1:50
    EcoSISTEM.update!(eco, 1month)
end

matrix = copy(eco.abundances.matrix)
@save ARGS[1] matrix

println("Saved reproducibility result with $(Threads.nthreads()) threads to $(ARGS[1])")
