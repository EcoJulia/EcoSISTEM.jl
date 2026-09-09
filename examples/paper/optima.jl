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

# The community of the island configuration run for its years, returning each species' total
# abundance at the end.
function _optima_run(opts, widths; temperature = 298.0K)
    island = island_configuration()
    eco = community(uniform_environment(PAPER_CELLS, PAPER_AREA,
                                        temperature = temperature),
                    opts, widths, individuals = island.individuals,
                    seed = PAPER_SEED)
    simulate!(eco, island.years, PAPER_TIMESTEP)
    return perspecies(eco)
end

# `n` trait values evenly spaced over `[low, high)`, each at the centre of its own slot, so that
# `edges` cutting the range into `n / perbin` bins hold `perbin` species each.
function _traitgrid(low, high, n)
    return low .+ (high - low) .* ((0:(n - 1)) .+ 0.5) ./ n
end

function optima_figure()
    n = island_configuration().numspecies

    # A: preferences spread +/- 3 widths about the landscape's 298 K, all 2 K wide, in 1.2 K bins
    # of ten species. B and C: every optimum at 298 K, widths from 0.1 K to 5.1 K in 0.5 K bins of
    # ten species. The original started at 0.0001 K; the floor of 0.1 K is a realism preference,
    # since a niche a ten-thousandth of a kelvin wide is not a plant. B has the landscape at 298 K,
    # C at 299 K.
    spread = _traitgrid(292.0K, 304.0K, n)
    graded = _traitgrid(0.1K, 5.1K, n)
    runs = [(opts = spread, widths = fill(2.0K, n), temperature = 298.0K),
        (opts = fill(298.0K, n), widths = graded, temperature = 298.0K),
        (opts = fill(298.0K, n), widths = graded, temperature = 299.0K)]
    totals = map(r -> _optima_run(r.opts, r.widths,
                                  temperature = r.temperature), runs)

    edges = collect(292.0:1.2:304.0) .* K
    counts = weighted_histogram(spread, totals[1], edges)
    centres = ustrip.(°C, (edges[1:(end - 1)] .+ edges[2:end]) ./ 2)
    panela = bar(centres, counts, bar_width = 1.2,
                 xlab = "Temperature preference (C)", ylab = "Abundance",
                 label = "", title = "A", titleloc = :left)

    edges = collect(0.1:0.5:5.1) .* K
    matched = weighted_histogram(graded, totals[2], edges)
    shifted = weighted_histogram(graded, totals[3], edges)
    centres = ustrip.(K, (edges[1:(end - 1)] .+ edges[2:end]) ./ 2)
    # The same y range on both, so the specialist-versus-generalist contrast reads across them.
    top = 1.05 * max(maximum(matched), maximum(shifted))
    panelb = bar(centres, matched, bar_width = 0.5, xlab = "Niche width (C)",
                 ylab = "Abundance", ylim = (0, top), label = "", title = "B",
                 titleloc = :left)
    panelc = bar(centres, shifted, bar_width = 0.5, xlab = "Niche width (C)",
                 ylab = "", ylim = (0, top), label = "", title = "C",
                 titleloc = :left)

    return plot(panela, panelb, panelc, layout = @layout([a; b c]),
                size = (1200, 1000), margin = 8Plots.mm, guidefontsize = 12,
                tickfontsize = 10, titlefontsize = 16)
end

timed("OptVarPanel") do
    return savefigure(optima_figure(), "OptVarPanel")
end
