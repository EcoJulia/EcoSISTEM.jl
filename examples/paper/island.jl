# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Figure 5, `IslandTesting.pdf`: the model tested on island ecosystems - abundance against the two
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
function _surface_run()
    island = PAPER.island
    env = island_environment(topology = Island())
    fractions = range(0.0, 1.0, length = island.cells)
    fullsun = env.supply.sunlight.matrix[1, 1]
    fullwater = env.supply.water.matrix[1, 1]
    for (i, f) in enumerate(fractions)
        env.supply.sunlight.matrix[i, :] .= f * fullsun     # solar varies by row (Y)
        env.supply.water.matrix[:, i] .= f * fullwater      # water varies by column (X)
    end
    n = island.species
    eco = island_community(env, fill(island_temperature(), n),
                           fill(island.width, n))
    simulate!(eco, island.duration, PAPER_TIMESTEP)
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
function _island_total(cells, area)
    island = PAPER.island
    n = island.species
    eco = island_community(island_environment(cells, area, topology = Island()),
                           fill(island_temperature(), n), fill(island.width, n))
    simulate!(eco, island.duration, PAPER_TIMESTEP)
    return sum(eco.abundances.matrix)
end

# The fraction of a heterogeneous pool of `n` species still present after the run, the pool's
# traits drawn as `parameters.toml`'s `testing` set says: widths and dispersal uniform, optima
# spread about the island's temperature by their widths, mortality (and birth equal to it) and body
# size folded normal - so coexistence cannot be an artefact of every species being identical.
function _survival(n, rng)
    pool = PAPER.testing
    widths = rand(rng, Uniform(ustrip.(K, pool.pool_width)...), n) .* K
    opts = island_temperature() .+
           widths .* range(-pool.pool_spread, pool.pool_spread, length = n)
    death = abs.(rand(rng, Normal(ustrip.(year^-1, pool.pool_death)...), n)) ./
            year
    eco = island_community(island_environment(topology = Island()), opts,
                           widths,
                           dispersal = rand(rng,
                                            Uniform(ustrip.(km,
                                                            pool.pool_dispersal)...),
                                            n) .* km,
                           birth = death, death = death,
                           sizes = abs.(rand(rng,
                                             Normal(ustrip.(m^2,
                                                            pool.pool_size)...),
                                             n)) .* m^2)
    simulate!(eco, PAPER.island.duration, PAPER_TIMESTEP)
    return 100 * count(>(0), perspecies(eco)) / n
end

# The four panels' simulations, then the figure. B is the same 100 km^2 as 1, 4, 25 and 100 cells, where total abundance is a
# property of the area and its supply, so the bars should be equal; C is areas of 10 to 100 km^2
# on a 10 x 10 grid, where supply is per unit area, so the bars should be proportional to it; D is
# the percentage surviving from pools of 100 to 5,000 species, ten replicates, mean and standard
# deviation - the slow panel, ten runs of 5,000 species on 100 cells for ten years.
function island_figure()
    island, pool = PAPER.island, PAPER.testing
    surface = _surface_run()
    sides = pool.sides
    bysides = map(s -> _island_total(s, island.area), sides)
    areas = ustrip.(km^2, pool.areas)
    byarea = map(a -> _island_total(island.cells, a * km^2), areas)
    pools = pool.pools
    replicates = pool.replicates
    units = [(n, r) for n in pools for r in 1:replicates]
    survived = map(units) do (n, r)
        s = _survival(n, Xoshiro(PAPER_SEED + r))
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

timed("IslandTesting") do
    return savefigure(island_figure(), "IslandTesting")
end
