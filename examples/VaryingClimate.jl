# SPDX-License-Identifier: LGPL-3.0-or-later
#
# **A layer that changes while the simulation runs.** Three ways of saying so, side by side, on the
# same landscape and the same species — so what differs between the runs is only how the environment
# moved.
#
#   1. **Static** — the control. Nothing changes.
#   2. **A steady drift** — `Varying(spec, IncrementBy(rate))`, a warming climate.
#   3. **A seasonal cycle** — `Varying(spec, OffsetBy(PatternedChange(...)))`, a temperature that
#      rises and falls through the year and returns to where it began.
#
# **A layer that varies in time is not a different type.** It is an ordinary layer holding the
# values current *now*, plus an [`AbstractLayerChange`](@ref) that is a **pure function of elapsed
# time**. That purity is a requirement rather than a tidiness: layers are updated redundantly on
# every MPI rank, so a change that consulted the ecosystem — or drew a random number — would let
# ranks silently disagree. It is also why twelve one-month steps and one twelve-month step land in
# the same place.
#
# **Contrast with an `Intervention`**, which `examples/interventions.jl` covers: that is the
# *other* mechanism, for what happens to the **ecosystem** — abundances, the active mask — rather
# than to one layer's values. A layer change is a function of time; an intervention is an edit.
#
# The clock advances *before* the layers change, so a change sees the time it is changing **to**;
# and the layer update runs *after* the demographics, so step `n`'s dynamics see the environment as
# of step `n − 1`.
#
# Run it directly for a longer, more visible run:
#
#     julia --project=examples examples/VaryingClimate.jl

# **A module, deliberately** — `test/extras_examples.jl` includes every top-level example into one
# module, and more than one of them defines names like `runscale()`.
module VaryingClimate

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

const SMALL = get(ENV, "ECOSISTEM_SCALE", "large") == "small"
const YEARS = SMALL ? 4year : 20year
const TIMESTEP = 1month_mean_duration

# Non-square deliberately: a square grid cannot see a y/x transposition.
const AREA = StudyArea(extent = (20.0km, 30.0km), cellsize = 5.0km,
                       verbosity = :silent)

const SUNLIGHT = UniformSpec(1.0e4kJ / (km^2 * day), axis = SolarRadiation)

# The three regimes. All start at the same 285 K, so any divergence is the change and nothing else.
const BASE = UniformSpec(285.0K, axis = Temperature)

# `IncrementBy` accumulates a **rate** each step — the only shape that takes a bare constant,
# because "the same increment every step" is the one thing a constant can unambiguously mean.
const WARMING = Varying(BASE, IncrementBy(0.1K / year))

# `OffsetBy` displaces the layer from its own values rather than accumulating, which is what makes
# it the right shape for something that must come back: a cycle that accumulated would drift.
# `PatternedChange(amplitude, period)` is a sinusoid, so this returns to 285 K every year.
const SEASONAL = Varying(BASE, OffsetBy(PatternedChange(8.0K, 1year)))

# One species, wide enough to survive all three so the comparison is about where it thrives rather
# than whether it lives.
function _species()
    return build_species(20, tolerance = (285.0K, 6.0K),
                         toleranceaxis = Temperature,
                         demand = 5.0kJ / day, demandaxis = SolarRadiation,
                         dispersal = 5.0km, abundance = 200_000, seed = 1)
end

# Run one regime and report where it ended. Seeded identically each time, so the runs are comparable.
function _run(name, regime)
    environment = GridHabitat(regime = regime, supply = SUNLIGHT,
                              area = AREA)
    eco = build_ecosystem(_species(), environment, seed = 1)
    simulate!(eco, YEARS, TIMESTEP)
    # Read the temperature off the **built habitat** after the run, not off the spec — the layer
    # holds the values current now, which is the whole point.
    #
    # And read the clock rather than assuming it says `YEARS`. `simulate!` takes
    # `length((0s):timestep:duration)` steps — an **inclusive** range from zero — so a 4-year run at
    # monthly steps advances **49** months, not 48. A layer change is a function of elapsed time, so
    # that one extra step is visible in the values: worth knowing before comparing against a period.
    return (name = name, abundance = sum(eco.abundances.matrix),
            temperature = first(environment.regime.matrix),
            elapsed = EcoSISTEM.simulationtime(eco))
end

# `println`, not `@info`: logging goes to stderr and `test/extras_examples.jl` wraps each example
# in `@test_nowarn`, which fails on *any* stderr output.
println("Three climates over $(YEARS), from the same 285.0 K start:")

const RESULTS = [_run("static", BASE), _run("warming", WARMING),
    _run("seasonal", SEASONAL)]

for r in RESULTS
    println("  ", rpad(r.name, 9), "ends at ",
            round(typeof(1.0K), r.temperature, digits = 2), " after ",
            round(typeof(1.0year), r.elapsed, digits = 2), ", ",
            r.abundance, " individuals")
end

# --- what the three must show -------------------------------------------------------------

const STATIC, WARM, SEASON = RESULTS

# The control has not moved at all.
STATIC.temperature == 285.0K ||
    error("the static regime changed, and nothing should have moved it")

# Warming has accumulated its rate over the elapsed time — read from the clock, not assumed.
const EXPECTED = 285.0K + 0.1K / year * WARM.elapsed
isapprox(WARM.temperature, EXPECTED, atol = 0.01K) ||
    error("warming reached $(WARM.temperature), expected about $(EXPECTED)")

# **The property that distinguishes a cycle from a drift, and it is not "ends where it started"**
# — that is only true at whole periods, and this run stops one step past one. The real distinction is
# that a cycle stays **bounded** by its amplitude however long it runs, while `IncrementBy` grows
# without limit. Run it for a century and this still holds; the warming assertion above does not.
abs(SEASON.temperature - 285.0K) <= 8.0K ||
    error("the seasonal regime left its amplitude: $(SEASON.temperature)")

# …and it is exactly the sinusoid evaluated at the clock, which is what "a pure function of elapsed
# time" means. `sinpi(2φ)` is `PatternedChange`'s own default shape.
const PHASE = SEASON.elapsed / 1year
isapprox(SEASON.temperature, 285.0K + 8.0K * sinpi(2 * PHASE), atol = 0.01K) ||
    error("the seasonal regime is not its declared shape at elapsed $(SEASON.elapsed)")

println("\n✓ static held, warming accumulated its rate, and the cycle stayed bounded by its",
        " amplitude —\n  landing exactly where its shape says it should at $(round(typeof(1.0year), SEASON.elapsed, digits = 2)).")

end
