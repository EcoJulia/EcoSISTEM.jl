# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Figure 6, `DispersalSD.pdf`: two species seeded on opposite edges of an island, after 50 years,
# for mean dispersal distances of 0.5 to 4 km. A wider kernel drains the edges and fills the middle.
#
# Included by `examples/paper.jl`, never run directly.

using Plots

# Total abundance per cell after the run, with species 1 started at `founders` in every cell of the
# first column and species 2 in the last (the grid is species x Y x X).
function _dispersal_map(distance)
    founders = PAPER.dispersal.founders
    eco = island_community(island_environment(topology = Island()),
                           fill(island_temperature(), 2),
                           fill(PAPER.island.width, 2), dispersal = distance,
                           individuals = 0)
    eco.abundances.grid[1, :, 1] .= founders
    eco.abundances.grid[2, :, end] .= founders
    simulate!(eco, PAPER.dispersal.duration, PAPER_TIMESTEP)
    return percell(eco)
end

function dispersal_figure()
    # One colour scale across the four panels, at the published 50-year ceiling; the five-year
    # small-scale run is far below it, so there each panel scales itself.
    clim = paper_small() ? :auto : (0, 15_000)
    maps = map(_dispersal_map, PAPER.dispersal.distances)
    cells = PAPER.island.cells
    panels = map(zip(maps, ["A", "B", "C", "D"])) do (m, title)
        return heatmap(1:cells, 1:cells, m, clim = clim,
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
