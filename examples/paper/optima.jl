# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Figure 7, `OptVarPanel.pdf`: how abundance is distributed over temperature preference, and over
# niche width when the landscape matches every species' optimum and when it is 1 K off. Torus
# topology, as the original.
#
# Each panel is a histogram of individuals over the trait - every species contributes its final
# abundance at its own trait value - so the y axis is abundance, not a count of species. The
# species sit on a grid that puts exactly ten in every bin, so a bar's height is not swayed by how
# many species happened to fall in it.
#
# Included by `examples/paper.jl`, never run directly.

using Plots

# The island's community run for its duration, returning each species' total abundance at the end.
function _optima_run(opts, widths; temperature = island_temperature())
    eco = island_community(island_environment(temperature = temperature), opts,
                           widths)
    simulate!(eco, PAPER.island.duration, PAPER_TIMESTEP)
    return perspecies(eco)
end

# `n` trait values evenly spaced over `[low, high)`, each at the centre of its own slot, so that
# `edges` cutting the range into `n / perbin` bins hold `perbin` species each.
function _traitgrid(low, high, n)
    return low .+ (high - low) .* ((0:(n - 1)) .+ 0.5) ./ n
end

function optima_figure()
    n = PAPER.island.species
    optima = PAPER.optima
    temperature = island_temperature()
    width = PAPER.island.width

    # A: preferences spread evenly either side of the landscape's temperature, every species of the
    # island's width, in `bins` bins. B and C: every optimum at that temperature, widths spread evenly
    # over `optima.widths`, in `bins` bins; a niche of no width is not a plant, hence the floor. B
    # has the landscape at that temperature, C `offset` warmer.
    low, high = temperature - optima.spread, temperature + optima.spread
    spread = _traitgrid(low, high, n)
    graded = _traitgrid(optima.widths..., n)
    runs = [(opts = spread, widths = fill(width, n), temperature = temperature),
        (opts = fill(temperature, n), widths = graded,
         temperature = temperature),
        (opts = fill(temperature, n), widths = graded,
         temperature = temperature + optima.offset)]
    totals = map(r -> _optima_run(r.opts, r.widths,
                                  temperature = r.temperature), runs)

    step = (high - low) / optima.bins
    edges = collect(ustrip(K, low):ustrip(K, step):ustrip(K, high)) .* K
    counts = weighted_histogram(spread, totals[1], edges)
    centres = ustrip.(°C, (edges[1:(end - 1)] .+ edges[2:end]) ./ 2)
    panela = bar(centres, counts, bar_width = ustrip(K, step),
                 xlab = "Temperature preference (C)", ylab = "Abundance",
                 label = "", title = "A", titleloc = :left)

    lowest, widest = optima.widths
    step = (widest - lowest) / optima.bins
    edges = collect(ustrip(K, lowest):ustrip(K, step):ustrip(K, widest)) .* K
    matched = weighted_histogram(graded, totals[2], edges)
    shifted = weighted_histogram(graded, totals[3], edges)
    centres = ustrip.(K, (edges[1:(end - 1)] .+ edges[2:end]) ./ 2)
    # The same y range on both, so the specialist-versus-generalist contrast reads across them.
    top = 1.05 * max(maximum(matched), maximum(shifted))
    panelb = bar(centres, matched, bar_width = ustrip(K, step),
                 xlab = "Niche width (C)",
                 ylab = "Abundance", ylim = (0, top), label = "", title = "B",
                 titleloc = :left)
    panelc = bar(centres, shifted, bar_width = ustrip(K, step),
                 xlab = "Niche width (C)",
                 ylab = "", ylim = (0, top), label = "", title = "C",
                 titleloc = :left)

    return plot(panela, panelb, panelc, layout = @layout([a; b c]),
                size = (1200, 1000), margin = 8Plots.mm, guidefontsize = 12,
                tickfontsize = 10, titlefontsize = 16)
end

timed("OptVarPanel") do
    return savefigure(optima_figure(), "OptVarPanel")
end
