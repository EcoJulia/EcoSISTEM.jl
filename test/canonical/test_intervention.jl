# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Canonical results for **ecosystem-level interventions** - the second mechanism of change, and one
# with no blessed coverage at all before this file.
#
# **Why it needs its own canonical run.** An intervention is deliberately *not* a layer change: it
# mutates the ecosystem (the active mask, abundances) rather than being a pure function of elapsed
# time, so it must be applied once and identically everywhere. Two things that can only be pinned by
# blessing a result:
#
#   - **the selection stream** - `RandomCells` draws from a counter-based stream
#     `hash((seed, :intervention, k, step))`, so which cells are chosen is reproducible and a change
#     to that scheme silently redistributes everything;
#   - **the ordering** - the clock advances, *then* interventions run, *then* layers update, so a
#     `SetChange` bites the same step rather than one late. Moving an intervention either side of the
#     layer update changes the answer without changing any single component's arithmetic.
#
# Runs on the non-square, spatially varying fixture (`test/varyingcase.jl`), so *which* cells are
# deactivated matters to the outcome - on a uniform grid every choice of cells is equivalent and the
# selection stream could be replaced wholesale without moving a number.

module CanonicalIntervention

using Test
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

include("canonical.jl")
using .Canonical
include(joinpath(@__DIR__, "..", "varyingcase.jl"))

const NCELLS = VARYING_NX * VARYING_NY

@testset "canonical: interventions" begin
    # A fixed *count* of random cells, not a rate, so the number deactivated is exact and any
    # movement in these results is about *which* cells were chosen, never how many.
    eco = varying_ecosystem()
    kill = Intervention(AtTime(6.0month_mean_duration), RandomCells(20),
                        Deactivate())
    simulate!(eco, 2year, 1month_mean_duration, intervention = kill)
    abun = eco.abundances.grid

    canonical("intervention/total_abundance", sum(abun))
    canonical("intervention/abundance_by_species",
              vec(sum(abun, dims = (2, 3))))
    # The spatial profiles are the point: deactivating 20 of 84 cells is a *spatial* event, so a
    # change to which cells the stream picks moves these even when the total barely shifts.
    canonical("intervention/abundance_by_row", vec(sum(abun, dims = (1, 3))))
    canonical("intervention/abundance_by_column", vec(sum(abun, dims = (1, 2))))
    canonical("intervention/active_cells", count(eco.habitat.active))

    # --- properties ------------------------------------------------------------------------------
    # Exactly 20 cells were deactivated, and they stayed deactivated.
    @test count(eco.habitat.active) == NCELLS - 20

    # **`Deactivate` kills what lives in the cell, and must**: a deactivated cell is skipped by the
    # hot loop, so anything left in it would neither breed nor die. Asserted rather than assumed -
    # this is a model requirement, not an implementation detail.
    # Indexed through the flat `matrix` (species × location) rather than the 3-D `grid`: the two
    # are views of the same memory, and location is column-major over `(Y, X)`, so the flattened mask
    # lines up exactly. Masking the grid directly needs a 2-D index over dims 2:3 and is easy to get
    # subtly wrong.
    deadcells = .!vec(parent(eco.habitat.active))
    @test all(iszero, eco.abundances.matrix[:, deadcells])
    @test sum(abun) > 0                      # ...but the run as a whole survived

    # The selection must be reproducible from the seed alone - the premise of blessing any of the
    # numbers above.
    second = varying_ecosystem()
    simulate!(second, 2year, 1month_mean_duration,
              intervention = Intervention(AtTime(6.0month_mean_duration),
                                          RandomCells(20), Deactivate()))
    @test second.habitat.active == eco.habitat.active
    @test second.abundances.grid == abun

    # ...and it must be a *different* selection from a different seed, or "reproducible" would be
    # satisfied by a constant, which is the failure mode this guards.
    other = varying_ecosystem(seed = 12345)
    simulate!(other, 2year, 1month_mean_duration,
              intervention = Intervention(AtTime(6.0month_mean_duration),
                                          RandomCells(20), Deactivate()))
    @test other.habitat.active != eco.habitat.active
end

@testset "canonical: intervention ordering" begin
    # `AddAbundance` on a schedule, blessed separately from the deactivation run so that a change
    # to *when* interventions are applied relative to the clock and the layer update is visible on
    # its own rather than tangled with a spatial selection.
    eco = varying_ecosystem()
    add = Intervention(AtTimes([3.0month_mean_duration, 9.0month_mean_duration]),
                       AllCells(), AddAbundance(1, 5))
    simulate!(eco, 1year, 1month_mean_duration, intervention = add)

    canonical("intervention/added_species_total",
              sum(eco.abundances.grid[1, :, :]))
    canonical("intervention/added_total_abundance", sum(eco.abundances.grid))

    @test sum(eco.abundances.grid) > 0
    @test all(>=(0), eco.abundances.grid)
end

end
