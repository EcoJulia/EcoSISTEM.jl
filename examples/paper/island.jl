# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Figure 5, `Abundance.pdf`: the model tested on island ecosystems - abundance against the two
# resources together, invariance to grid resolution, scaling with area, and survival in large
# species pools. Island topology throughout.
#
# Included by `examples/paper.jl`, never run directly.

using Distributions
using Random
using Plots

# A: solar supply graded down the rows and water across the columns of one landscape, so every
# cell is a distinct (solar, water) pair and abundance shows which resource binds where. Both
# ranges start at zero, so the first row and column are uninhabitable, deliberately. Returns the
# axes and the per-cell totals; `_surface_panel` draws them.
function _surface_run(island)
    env = uniform_environment(PAPER_CELLS, PAPER_AREA, topology = Island())
    fractions = range(0.0, 1.0, length = PAPER_CELLS)
    fullsun = env.supply.sunlight.matrix[1, 1]
    fullwater = env.supply.water.matrix[1, 1]
    for (i, f) in enumerate(fractions)
        env.supply.sunlight.matrix[i, :] .= f * fullsun     # solar varies by row (Y)
        env.supply.water.matrix[:, i] .= f * fullwater      # water varies by column (X)
    end
    n = island.numspecies
    eco = community(env, fill(298.0K, n), fill(2.0K, n),
                    individuals = island.individuals, seed = PAPER_SEED)
    simulate!(eco, island.years, PAPER_TIMESTEP)
    # A supply layer holds a rate per cell in the resource's canonical unit; the axes show it as
    # the amount a cell receives a month.
    return (water = fractions .*
                    ustrip(Unitful.L, fullwater * 1month_mean_duration),
            sun = fractions .* ustrip(kJ, fullsun * 1month_mean_duration),
            abundance = percell(eco))
end

# `percell` is `(y, x)`, so water is the x axis and solar the y.
function _surface_panel(surface)
    return heatmap(surface.water, surface.sun, surface.abundance,
                   xlab = "Water supply (L/month per cell)",
                   ylab = "Solar supply (kJ/month per cell)",
                   colorbar_title = "Total abundance", title = "A",
                   titleloc = :left)
end

# Total abundance after the run on a uniform island of `cells` a side over `area`.
function _island_total(island, cells, area)
    n = island.numspecies
    eco = community(uniform_environment(cells, area, topology = Island()),
                    fill(298.0K, n), fill(2.0K, n),
                    individuals = island.individuals, seed = PAPER_SEED)
    simulate!(eco, island.years, PAPER_TIMESTEP)
    return sum(eco.abundances.matrix)
end

# The fraction of a heterogeneous pool of `n` species still present after the run: widths uniform
# on 1-5 K, optima spread +/- 3 widths about 298 K, dispersal uniform on 0.6-2.4 km, mortality
# (and birth equal to it) folded normal about 0.15/year, and body size normal about 1 m^2, so
# coexistence cannot be an artefact of every species being identical.
function _survival(island, n, rng)
    widths = rand(rng, Uniform(1.0, 5.0), n) .* K
    opts = 298.0K .+ widths .* range(-3, 3, length = n)
    death = abs.(rand(rng, Normal(0.15, 0.135), n)) ./ year
    eco = community(uniform_environment(PAPER_CELLS, PAPER_AREA,
                                        topology = Island()), opts, widths,
                    dispersal = rand(rng, Uniform(0.6, 2.4), n) .* km,
                    birth = death, death = death,
                    sizes = abs.(rand(rng, Normal(1.0, 0.1), n)) .* m^2,
                    individuals = island.individuals, seed = PAPER_SEED)
    simulate!(eco, island.years, PAPER_TIMESTEP)
    return 100 * count(>(0), perspecies(eco)) / n
end

# The four panels' simulations, then the figure. B is the same 100 km^2 as 1, 4, 25 and 100 cells, where total abundance is a
# property of the area and its supply, so the bars should be equal; C is areas of 10 to 100 km^2
# on a 10 x 10 grid, where supply is per unit area, so the bars should be proportional to it; D is
# the percentage surviving from pools of 100 to 5,000 species, ten replicates, mean and standard
# deviation - the slow panel, ten runs of 5,000 species on 100 cells for ten years.
function island_figure()
    island = island_configuration()
    surface = _surface_run(island)
    sides = [1, 2, 5, 10]
    bysides = map(s -> _island_total(island, s, PAPER_AREA), sides)
    areas = [10.0, 20.0, 50.0, 100.0]
    byarea = map(a -> _island_total(island, PAPER_CELLS, a * km^2),
                 areas)
    pools = paper_small() ? [5, 10, 20] : [100, 500, 1000, 5000]
    replicates = paper_small() ? 1 : 10
    units = [(n, r) for n in pools for r in 1:replicates]
    survived = map(units) do (n, r)
        s = _survival(island, n, Xoshiro(PAPER_SEED + r))
        println("    pool of ", n, ", replicate ", r, ": ",
                round(s, digits = 1), "% survived")
        return s
    end

    panela = _surface_panel(surface)
    panelb = bar(string.(sides .^ 2), bysides, xlab = "Number of grid squares",
                 ylab = "Total abundance", label = "", title = "B",
                 titleloc = :left)
    panelc = bar(string.(areas), byarea, xlab = "Area (km^2)",
                 ylab = "Total abundance", label = "", title = "C",
                 titleloc = :left)
    bypool = [[s for ((n, _), s) in zip(units, survived) if n == pool]
              for pool in pools]
    means = mean.(bypool)
    sds = replicates > 1 ? std.(bypool) : zeros(length(pools))
    paneld = bar(string.(pools), means, yerr = sds, ylim = (0, 100),
                 xlab = "Number of species introduced",
                 ylab = "% Species survived", label = "", title = "D",
                 titleloc = :left)
    return plot(panela, panelb, panelc, paneld, layout = @layout([a b; c d]),
                size = (1600, 1200), margin = 8Plots.mm, guidefontsize = 12,
                tickfontsize = 10, titlefontsize = 16)
end

timed("Abundance") do
    return savefigure(island_figure(), "Abundance")
end
