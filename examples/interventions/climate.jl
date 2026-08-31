# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Climate scenarios from `examples/models/scenarios.jl`, recreated as declarative interventions.
#
# Each was a *callback* - `SimpleScenario(TempIncrease!, 1.0K/year)` - whose behaviour was
# invisible until it ran and which could not be dispatched on, validated or reported. Each is now a
# declaration: a schedule, a region and one of six operations.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

# A small warming/drying ecosystem: temperature *and* rainfall as a two-layer regime, so that the
# rainfall scenario has a layer to name.
function climate_ecosystem(; numspecies = configuration().numspecies,
                           individuals = configuration().individuals, seed = 1)
    env = GridHabitat(regime = (temperature = UniformSpec(285.0K,
                                                          axis = Temperature),
                                rainfall = UniformSpec(5.0mm / day,
                                                       axis = Precipitation)),
                      supply = UniformSpec(supply_density(individuals),
                                           axis = SolarRadiation),
                      area = study_area())
    # Both sides named, and `build_ecosystem` insists: a named regime against a positional
    # tolerance is refused outright, because it cannot then check them layer for layer.
    species = build_species(numspecies,
                            tolerance = (temperature = (285.0K, 5.0K),
                                         rainfall = (5.0mm / day,
                                                     2.0mm / day)),
                            toleranceaxis = (temperature = Temperature,
                                             rainfall = Precipitation),
                            demand = intervention_demand(),
                            demandaxis = SolarRadiation,
                            abundance = individuals,
                            seed = seed)
    return build_ecosystem(species, env, seed = seed)
end

# ---------------------------------------------------------------------------
# `TempIncrease!` - steady warming
# ---------------------------------------------------------------------------
# **One line, and two things fall away.** The original re-installed the same change on *every*
# step (a `SimpleScenario` runs every step), which was harmless but pointless; `AtTime(0.0s)`
# installs it once. And its trailing `matrix[matrix .< 0K] .= 0K` clamp is gone: a temperature's
# floor is now a property of the axis (`bounds(::Type{<:TemperatureAxis}) == (0.0K, nothing)`), enforced
# automatically - and `check_bounds` refuses a run that *would* breach it before the first step
# rather than clamping silently part-way through.
function temperature_increase(rate = 1.0K / year)
    return Intervention(AtTime(0.0s), AllCells(),
                        SetChange(:temperature,
                                  IncrementBy(rate)))
end

# ---------------------------------------------------------------------------
# `RainDecrease!` - steady drying
# ---------------------------------------------------------------------------
# **A rewrite, not a port, and the original was broken.** Its line
# `supply.two.matrix .= ustrip(depth) * L/day` *assigned the whole water supply to a constant* - one
# timestep's worth of change - rather than decrementing it, so with a decreasing rate that constant
# was negative and the clamp below zeroed the entire grid's water in one step. Nothing ever called
# it, which is why nobody noticed.
#
# Declared, the intent is unambiguous: rainfall *decrements* at a rate, and naming the layer is
# what `resetrate!` could never do - it reached the regime as a whole, so a two-layer environment had
# no way to say *which* variable was drying.
function rain_decrease(rate = -0.1mm / day / year)
    return Intervention(AtTime(0.0s),
                        AllCells(),
                        SetChange(:rainfall,
                                  IncrementBy(rate)))
end

# ---------------------------------------------------------------------------
# `TempFluct!` - seasonal oscillation
# ---------------------------------------------------------------------------
# The original stepped a 12-point sine table by `mod(currentstep/timestep, 12)`, so its result
# depended on the *number of steps taken* rather than on elapsed time - twelve one-month steps and
# one twelve-month step disagreed. `PatternedChange` is a pure function of elapsed time, so they
# agree. Numbers therefore differ from v0.4.0; that is the point, and `TempFluct` already carries
# a loud semantics-changed deprecation warning saying so.
function temperature_fluctuation(amplitude = 2.0K, period = 1.0year)
    return Intervention(AtTime(0.0s),
                        AllCells(),
                        SetChange(:temperature,
                                  OffsetBy(PatternedChange(amplitude,
                                                           period))))
end

# All three at once - an `InterventionSet` applies them in the order written. The old `MultiScenario`
# held exactly two, and hard-coded their types.
function changing_climate(; warming = 1.0K / year, drying = -0.1mm / day / year,
                          amplitude = 2.0K)
    return InterventionSet(temperature_increase(warming), rain_decrease(drying),
                           temperature_fluctuation(amplitude))
end
