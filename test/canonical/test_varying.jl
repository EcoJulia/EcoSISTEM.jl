# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Canonical results from a **non-uniform, non-square, time-varying** simulated ecosystem - the
# companion to `test_simulated.jl`, which is uniform, square and static.
#
# **This file exists because that one could not break.** Every blessed number in the suite came
# from an environment where every cell was identical and nothing changed over time, so no spatial
# arithmetic and no layer change was pinned by anything. See `test/varyingcase.jl` for the shape of
# the fixture and why it is orthogonal (regime varies down `Y`, supply across `X`).
#
# Additive by design: no key here replaces one in `test_simulated.jl`, so the historical baseline
# those numbers provide is untouched.

module CanonicalVarying

using Test
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

include("canonical.jl")
using .Canonical
include(joinpath(@__DIR__, "..", "varyingcase.jl"))

@testset "canonical: varying, non-square ecosystem" begin
    eco = varying_ecosystem()
    simulate!(eco, 2year, 1month_mean_duration)
    abun = eco.abundances.grid          # species × Y × X

    canonical("varying/total_abundance", sum(abun))
    canonical("varying/abundance_by_species", vec(sum(abun, dims = (2, 3))))

    # **The spatial signatures, and the reason this file was worth writing.** Row and column
    # profiles are blessed separately: a change that moves individuals around the grid while
    # preserving both the grand total and the per-species totals still moves one of these. On the
    # uniform square grid neither profile carried any information at all.
    canonical("varying/abundance_by_row", vec(sum(abun, dims = (1, 3))))
    canonical("varying/abundance_by_column", vec(sum(abun, dims = (1, 2))))

    # The end state of the *environment*, which is what pins the layer change: the regime warms at
    # a rate, so this is a function of elapsed time and would move if a change were ever applied per
    # step instead of per unit time.
    canonical("varying/final_regime_mean",
              ustrip(K,
                     sum(eco.habitat.regime.matrix) /
                     length(eco.habitat.regime.matrix)))

    # --- properties that must hold whatever the blessed numbers are -----------------------------
    # A blessed number says *something changed*; a property says *the model is still right*.
    # Re-blessing silences the first and must never silence the second.
    @test sum(abun) > 0
    @test all(>=(0), abun)

    # **The dimension-order guard.** The grid is deliberately 7 × 12, so a transposed index is a
    # shape error rather than a plausible wrong number.
    @test size(abun)[2:3] == (VARYING_NY, VARYING_NX)
    @test size(eco.habitat.regime.matrix) == (VARYING_NY, VARYING_NX)

    # The regime varies **down** the grid and the supply **across** it, so a transposition swaps
    # these two counts. Asserted on the built environment, before any simulation.
    @test length(unique(eco.habitat.regime.matrix)) == VARYING_NY
    @test length(unique(EcoSISTEM.values(eco.habitat.supply)[1].matrix)) ==
          VARYING_NX

    # **The warming is analytically predictable, so this blessed number is *checked* rather than
    # merely reproducible** - the strongest form a canonical entry can take. The spatial mean of a
    # 288-302 K gradient is 295 K, and `IncrementBy(0.5K/year)` adds 0.5 K per elapsed year.
    #
    # It is **25** months of warming after 24 steps, and that is the documented behaviour, not an
    # off-by-one: the clock advances *before* the layers change, so the last change sees the time it
    # is changing *to* (see CLAUDE.md). Pinning it here means that ordering cannot be altered silently.
    mean_regime = sum(eco.habitat.regime.matrix) /
                  length(eco.habitat.regime.matrix)
    @test mean_regime ≈ 295.0K + 0.5K / year * (25 * month_mean_duration) rtol=1e-6

    # Species sort along the regime's gradient: with optima spread across it, the warm end must
    # not hold the same assemblage as the cold end. This is the ecological content of the fixture,
    # and it is what a uniform environment cannot express.
    warm = sum(abun[:, 1, :], dims = 2)         # row 1 is the 302 K end
    cold = sum(abun[:, end, :], dims = 2)       # the 288 K end
    @test argmax(vec(warm)) != argmax(vec(cold))

    # Reproducibility is the premise, so it is asserted rather than assumed.
    second = varying_ecosystem()
    simulate!(second, 2year, 1month_mean_duration)
    @test second.abundances.grid == abun
end

end
