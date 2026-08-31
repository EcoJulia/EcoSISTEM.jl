# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Canonical results for the fixture the **distributed** test runs — blessed from a **serial** run.
#
# **Why this file exists.** `test/ext_EcoSISTEMMPIExt.jl` compares MPI at 1, 2 and 4 ranks against
# each other, and every one of those takes the duplicated hot loop in `ext/EcoSISTEMMPIExt/`. A
# change made there and not in `src/` therefore agrees with itself at every rank count and passes.
# Comparing a duplicate against itself cannot detect that it has diverged from the original.
#
# So both paths are pinned to **one** blessed number instead. This file records what the serial run
# produces; `test/SmallMPItest.jl` reads the same keys and asserts the gathered distributed matrix
# against them, read-only, at every rank count.
#
# **What that buys, and it is the whole point — it separates three failure modes:**
#
#   - serial *and* distributed both move off the blessed value: a **model change**. Explain it, then
#     re-bless deliberately.
#   - **only the distributed** result moves: an **MPI bug**, exactly A22, now caught at 2 and 4 ranks
#     rather than only where a serial run happens to sit alongside it.
#   - **only the serial** result moves: a serial bug, or an intended change whose distributed half
#     was forgotten.
#
# No MPI is needed to bless these: it is a serial run, so this set stays MPI-free and cheap.
#
# **Re-blessing `mpi/…` is only ever right when the serial keys moved for a reason already
# understood.** Re-blessing to make a red distributed test green would silently re-admit A22, which
# is the one thing this file is here to prevent.

module CanonicalMPIFixture

using Test
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

include("canonical.jl")
using .Canonical
include(joinpath(@__DIR__, "..", "varyingcase.jl"))

@testset "canonical: the distributed fixture, run serially" begin
    eco = mpifixture_ecosystem()
    simulate!(eco, MPIFIXTURE_BURNIN, MPIFIXTURE_TIMESTEP)
    abun = eco.abundances.grid              # species × Y × X

    canonical("mpi/total_abundance", sum(abun))
    canonical("mpi/abundance_by_species", vec(sum(abun, dims = (2, 3))))

    # The spatial signatures, blessed separately: a change that moves individuals around the grid
    # while preserving both the grand total and the per-species totals still moves one of these.
    # A distributed run partitions by species *and* by cells, so both axes are worth pinning.
    canonical("mpi/abundance_by_row", vec(sum(abun, dims = (1, 3))))
    canonical("mpi/abundance_by_column", vec(sum(abun, dims = (1, 2))))

    # --- properties that must hold whatever the blessed numbers are -----------------------------
    # A blessed number says *something changed*; a property says *the model is still right*.
    @test sum(abun) > 0
    @test all(>=(0), abun)

    # The dimension-order guard: 7 × 11 on purpose, so a transposed index is a shape error rather
    # than a plausible wrong number.
    @test size(abun)[2:3] == (VARYING_NY, VARYING_NX)
end

end
