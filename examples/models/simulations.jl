# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Running the experiments and recording diversity through time.
#
# **This is now a single do-block over the package's own `simulate_action!`.** The original had
# its own `runsim!`/`dispersalrun!`/`simulate_record_diversity!` loop that called `update!` and then
# `runscenario!` by hand, plus a JLD2 cache of every abundance matrix. `simulate_action!` already is
# that loop — "run, and call this every `interval`" — and it takes an `intervention`, so nothing
# here needs to know how the environment changes.
#
# It also drops three `runscenario!` methods the original defined. Those were **method piracy**:
# character-for-character copies of the package's own (`src/Scenarios.jl`), differing only in taking
# `Ecosystem` rather than `AbstractEcosystem` — which made them *more* specific, so they silently won
# dispatch while doing exactly what the package would have done anyway.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Diversity
using Diversity.Ecology

# **The measures come from Diversity.jl**, not from EcoSISTEM's own copies. The metacommunity
# partition (`norm_meta_alpha`…`meta_gamma`) takes the ecosystem and a vector of diversity orders;
# the classical indices in `Diversity.Ecology` take a plain species × cells abundance matrix. Only
# **`faith_pd` is deliberately NOT here.** These ecosystems come from `build_species`, which builds
# no phylogeny, so their similarity is `UniqueTypes` and there is no tree to span. A `faith_pd`
# defined as `meta_gamma(eco, 0.0)` accepts that happily and reports **species richness at q = 0,
# not Faith's PD**. Faith's PD is Diversity's, supplied for an ecosystem by `EcoSISTEMPhyloExt` and
# constrained to a species list that actually holds a `PhyloBranches`, so the meaningless call is
# refused rather than answered.
#
# EcoSISTEM's `meta_speciesrichness`/`meta_shannon`/`meta_simpson`/`sorenson`/`mean_abun`/
# `geom_mean_abun` — which the original recorded — are on their way out, so this deliberately reaches
# past them to `Diversity.Ecology`'s `richness`/`shannon`/`simpson`/`pielou`.
# The abundance matrix cannot be handed to `Diversity.Ecology` as it stands once anything
# deactivates a cell. A destroyed cell is an **all-zero column**, which is not a community at all:
# `simpson`/`pielou` return `NaN` for it and `richness` fails outright with `InexactError: Int64(NaN)`
# — which is exactly how the habitat-loss scenario first died. Dropping the empty columns asks the
# measure about the cells that still hold something, which is the intended question. Normalising at
# the same time also silences Diversity's "abundances not normalised to 1, correcting…" warning,
# which would otherwise be emitted on **every** recorded step.
#
# **Columns only — an empty row must be kept.** The matrix is species × cells, so a zero *column*
# is a cell holding nothing, while a zero *row* is a species extinct across the landscape. Measured:
# the empty row is handled correctly (`richness` reports 2 of 3, no `NaN`), and dropping it would be
# an outright bug — richness would then always equal the number of surviving species by
# construction, which is exactly the quantity a habitat-loss run exists to watch fall.
function _assemblage(eco)
    counts = eco.abundances.matrix
    occupied = counts[:, vec(sum(counts, dims = 1)) .> 0]
    total = sum(occupied)
    return iszero(total) ? occupied : occupied ./ total
end

const MEASURES = (richness = eco -> richness(_assemblage(eco)),
                  shannon = eco -> shannon(_assemblage(eco)),
                  simpson = eco -> simpson(_assemblage(eco)),
                  pielou = eco -> pielou(_assemblage(eco)),
                  norm_alpha = eco -> norm_meta_alpha(eco, [1.0]),
                  raw_alpha = eco -> raw_meta_alpha(eco, [1.0]),
                  norm_beta = eco -> norm_meta_beta(eco, [1.0]),
                  raw_beta = eco -> raw_meta_beta(eco, [1.0]),
                  norm_rho = eco -> norm_meta_rho(eco, [1.0]),
                  raw_rho = eco -> raw_meta_rho(eco, [1.0]),
                  gamma = eco -> meta_gamma(eco, [1.0]))

# The metacommunity measures return a one-row `DataFrame`; the ecological ones one row per
# subcommunity. Either way a single number summarises the landscape, which is what is plotted.
#
# Non-finite entries are skipped rather than propagated. A landscape being eaten away by habitat
# loss will eventually hand some measure a degenerate community, and one `NaN` would otherwise wipe
# out the whole trace from that step onwards — losing the very part of the run the scenario exists to
# show.
_summarise(result) = sum(filter(isfinite, result[!, :diversity]))

"""
    diversity_through_time(eco; times, interval, timestep, intervention = nothing)

Run `eco` for `times`, recording every measure in [`MEASURES`](@ref) every `interval`.

Returns a named tuple of `(times = …, measure = vector, …)` — one vector per measure, over the
recorded time points — so a caller reads `result.gamma` rather than remembering a column index.
"""
function diversity_through_time(eco; times = 10year,
                                interval = 1year,
                                timestep = 1month_mean_duration,
                                intervention = nothing)
    points = length((0 * interval):interval:times)
    recorded = map(_ -> Float64[], MEASURES)
    simulate_action!(eco, times, interval, timestep,
                     intervention = intervention) do _
        for (name, measure) in pairs(MEASURES)
            push!(recorded[name], _summarise(measure(eco)))
        end
    end
    stamps = collect((0 * interval):interval:times)[1:min(points,
                                                          length(recorded.gamma))]
    return merge((times = stamps,), recorded)
end
