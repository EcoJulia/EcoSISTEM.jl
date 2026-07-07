# SPDX-License-Identifier: LGPL-3.0-or-later

## Multithreaded reproducibility check, run as a subprocess so that it is
## guaranteed to exercise the threaded `update!` path with more than one thread
## (see the "Multithreaded reproducibility" testset in test_Generate.jl).

using EcoSISTEM
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Random

include(EcoSISTEM.path("TestCases.jl"))

# Seed before construction: the constructor draws the initial abundances, so the
# seed must be set first for the whole run to be reproducible.
function run_sim(seed)
    Random.seed!(seed)
    eco = Test1Ecosystem()
    for _ in 1:50
        EcoSISTEM.update!(eco, 1month)
    end
    return copy(eco.abundances.matrix)
end

Threads.nthreads() > 1 ||
    error("reproducibility check must run with more than one thread, got $(Threads.nthreads())")

a = run_sim(1234)
b = run_sim(1234)

a == b ||
    error("Non-reproducible results across seeded repeats with $(Threads.nthreads()) threads")

println("Reproducible across seeded repeats with $(Threads.nthreads()) threads")
