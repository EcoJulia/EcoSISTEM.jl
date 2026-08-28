# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Runs every scenario recreation in `examples/interventions/`, so that they are exercised rather than
# merely written. Discovered by `test/extras_examples.jl`, which includes each `examples/test_*.jl`.
#
# These recreate the scenarios of `examples/models/scenarios.jl` using the modern declarative
# interface. They deliberately do **not** reproduce the old numbers: two of the originals drew from
# the global RNG (so were not reproducible at all), and `TempFluct!` indexed a sine table by step
# count rather than elapsed time. The old files stay where they are and keep running, so nothing is
# lost while these are proven.

# **A module, deliberately.** `test/extras_examples.jl` includes every top-level example into
# **one** module, and more than one of them defines `runscale()` and `CONFIGURATIONS`. Julia 1.12
# lets a later `const` silently overwrite an earlier one — no warning — so without the wrapper these
# examples would quietly reconfigure each other depending on the order they happened to run in.

module Interventions

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

const _INTERVENTIONS = joinpath(@__DIR__, "interventions")

# First: it defines the scale and the shared landscape the three below build on.
include(joinpath(_INTERVENTIONS, "configurations.jl"))
include(joinpath(_INTERVENTIONS, "climate.jl"))
include(joinpath(_INTERVENTIONS, "habitat_loss.jl"))
include(joinpath(_INTERVENTIONS, "invasion.jl"))

# --- climate: warming, drying and a seasonal cycle, together -----------------
let eco = climate_ecosystem()
    before = eco.habitat.regime.temperature.matrix[1, 1]
    simulate!(eco, 5.0year, 1.0month_mean_duration,
              intervention = changing_climate())
    @assert eco.habitat.regime.temperature.matrix[1, 1] > before
    @assert all(>=(0.0mm / day), eco.habitat.regime.rainfall.matrix)
end

# --- habitat loss: destroying a cell kills what lives there -------------------
let eco = loss_ecosystem()
    simulate!(eco, 5.0year, 1.0month_mean_duration,
              intervention = random_loss(0.05 / year))
    lost = findall(.!vec(parent(eco.habitat.active)))
    @assert !isempty(lost)
    @assert iszero(sum(eco.abundances.matrix[:, lost]))
end

# --- …and it replays exactly, which the original never did --------------------
let run = () -> begin
        eco = loss_ecosystem()
        simulate!(eco, 1.0year, 1.0month_mean_duration,
                  intervention = clustered_loss(0.05 / year))
        findall(vec(parent(eco.habitat.active)))
    end
    @assert run() == run()
end

# --- land conversion: clear, then plant, on ONE resolved region ---------------
let eco = loss_ecosystem()
    crop = length(eco.spplist.names)
    simulate!(eco, 2.0year, 1.0month_mean_duration,
              intervention = convert_to_crop(crop, 500))
    converted = findall(.!vec(parent(eco.habitat.active)))
    @assert !isempty(converted)
    @assert all(eco.abundances.matrix[crop, converted] .>= 500)
end

# --- invasion: a species ARRIVES, carrying its own niche ----------------------
let outcomes = map((:generalist, :specialist)) do kind
        eco = invasion_ecosystem()
        # Room for the arrival, since the recording is sized before the run.
        storage = generate_storage(eco, 6, 1,
                                   maxspecies = configuration().numspecies)
        simulate_record!(storage, eco, 5.0year, 1.0year,
                         1.0month_mean_duration,
                         intervention = invasion(kind,
                                                 at = 1.0month_mean_duration))
        @assert length(eco.spplist.names) == configuration().numspecies
        @assert last(eco.spplist.names) == string(kind)
        @assert !last(eco.spplist.native)
        return sum(eco.abundances.matrix[end, :])
    end
    # The generalist, with the widest niche in the assemblage, outlasts the narrow specialist.
    @assert first(outcomes) > last(outcomes)
end

end
