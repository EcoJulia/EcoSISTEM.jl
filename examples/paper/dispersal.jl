# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Figure 6, `DispersalSD.pdf`: two species seeded on opposite edges of an island, after 50 years,
# for mean dispersal distances of 0.5 to 4 km. A wider kernel drains the edges and fills the middle.
#
# Included by `examples/paper.jl`, never run directly.

using Plots

# Total abundance per cell after the run, with species 1 started at 100 in every cell of the first
# column and species 2 in the last (the grid is species x Y x X).
function _dispersal_map(distance, years)
    eco = community(uniform_environment(PAPER_CELLS, PAPER_AREA,
                                        topology = Island()), fill(298.0K, 2),
                    fill(2.0K, 2), dispersal = distance, individuals = 0,
                    seed = PAPER_SEED)
    eco.abundances.grid[1, :, 1] .= 100
    eco.abundances.grid[2, :, end] .= 100
    simulate!(eco, years, PAPER_TIMESTEP)
    return percell(eco)
end

function dispersal_figure()
    years = paper_small() ? 5year : 50year
    # One colour scale across the four panels, at the published 50-year ceiling; the five-year
    # small-scale run is far below it, so there each panel scales itself.
    clim = paper_small() ? :auto : (0, 15_000)
    distances = [0.5, 1.0, 2.0, 4.0] .* km
    maps = map(d -> _dispersal_map(d, years), distances)
    panels = map(zip(maps, ["A", "B", "C", "D"])) do (m, title)
        return heatmap(1:PAPER_CELLS, 1:PAPER_CELLS, m, clim = clim,
                       xlab = "Distance (km)", ylab = "Distance (km)",
                       title = title, titleloc = :left)
    end
    return plot(panels..., layout = @layout([a b; c d]), link = :both,
                size = (1200, 1000), margin = 8Plots.mm, guidefontsize = 12,
                tickfontsize = 10, titlefontsize = 16)
end

timed("DispersalSD") do
    return savefigure(dispersal_figure(), "DispersalSD")
end
