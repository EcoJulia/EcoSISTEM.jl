# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The published worked examples: whole simulations exploring temperature optima, niche widths,
# resource supply, area, grid resolution, dispersal distance and large species pools, followed by
# diversity through time under the climate, habitat-loss and invasion scenarios.
#
# **What this is for, against its neighbours.** `examples/biodiversity.jl` asserts the *end state*
# of the same ecology as property tests; this follows the same investigations *through time* with a
# full suite of diversity measures and draws the figures. Same science, different output.
#
# **It never runs under the test suite or in CI.** `test/extras_examples.jl` executes every
# top-level `examples/*.jl`, but `runtests.jl` sets `ECOSISTEM_SCALE=small` and this file skips
# itself whenever it is set. These are 10–50-year runs over hundreds of species that write figures —
# far too slow for a routine gate. Run it deliberately:
#
#     julia --project=examples examples/models.jl
#
# Everything it needs lives in `examples/models/`:
#
#   - `models/ecosystems.jl`  — `uniform_environment`/`community`, on `GridHabitat`/
#     `build_species`/`build_ecosystem`
#   - `models/simulations.jl` — `diversity_through_time`, a do-block over `simulate_action!`, with
#     the measures taken from **Diversity.jl** rather than EcoSISTEM's own
#   - `models/experiments.jl` — the nine investigations, and the figures
#
# **The scenarios are gone, and that is the point of the rewrite.** `models/scenarios.jl` defined
# seven mutating callbacks; six of them were referenced by nothing at all, the seventh
# (`TempIncrease!`) only ever at a rate of `0.0K/10year`, and it additionally pirated three
# `runscenario!` methods from the package — character-for-character copies that differed only in
# being more specific, so they silently won dispatch. All of it is deleted. The scenarios are now
# declarative `Intervention`s, reused verbatim from `examples/interventions/`.
#
# **A module, deliberately.** `test/extras_examples.jl` includes every top-level example into one
# module, and `examples/interventions.jl` defines the same names this file reuses. Without the
# wrapper the two would redefine each other whenever both ran — which `@test_nowarn` would catch as
# stderr output.

module Models

if get(ENV, "ECOSISTEM_SCALE", "large") == "small"
    # `println`, not `@info` — see the note in `biodiversity.jl`: logging goes to stderr, which
    # `@test_nowarn` in `test/extras_examples.jl` treats as a failure.
    println("Skipping the published model examples: they are 10–50-year runs that draw figures, " *
            "and never execute under the test suite. Set ECOSISTEM_SCALE=large to run them.")
else
    using EcoSISTEM
    using EcoSISTEM.Units
    using Unitful
    using Unitful.DefaultSymbols
    using Distributions
    using Random
    using Plots

    gr()

    include(joinpath(@__DIR__, "models", "ecosystems.jl"))
    include(joinpath(@__DIR__, "models", "simulations.jl"))
    # Reused as they stand — these *are* the modern replacements for the deleted scenarios.
    # Their `configurations.jl` first: it carries the scale and the shared landscape.
    include(joinpath(@__DIR__, "interventions", "configurations.jl"))
    include(joinpath(@__DIR__, "interventions", "climate.jl"))
    include(joinpath(@__DIR__, "interventions", "habitat_loss.jl"))
    include(joinpath(@__DIR__, "interventions", "invasion.jl"))
    include(joinpath(@__DIR__, "models", "experiments.jl"))
end

end
