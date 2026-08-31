# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The experiments themselves - temperature optima, niche widths, resource supply, area, grid
# resolution, dispersal distance, large species pools, and the climate/habitat-loss/invasion
# scenarios.
#
# Included by `examples/models.jl`, never run directly: it depends on the builders in
# `ecosystems.jl`, the recorder in `simulations.jl`, and the interventions in
# `examples/interventions/`, which that file includes first.
#
# **What distinguishes this from `examples/biodiversity.jl`** is that the same investigations are
# followed *through time* with a full suite of diversity measures, and plotted - where
# `biodiversity.jl` asserts the end state. The ecology is the same; the output is a figure rather
# than a test.

const FIGURES = mktempdir()

# Save under `FIGURES` and report where, so a run leaves something to look at rather than depending
# on an interactive display.
function _savefig(plt, name)
    path = joinpath(FIGURES, name)
    savefig(plt, path)
    println("    wrote ", path)
    return path
end

const CELLS = 10
const AREA = 100.0km^2
const NUMSPECIES = 100

println("Writing figures to ", FIGURES)

# --- 1. different temperature optima -----------------------------------------
# 100 species differing only in preferred temperature, spread ±3 widths around the landscape's own
# 298 K. Those matched to the environment should end up commonest.
println("  1. temperature optima ...")
let widths = fill(2.0K, NUMSPECIES)
    opts = 298.0K .+ widths .* range(-3, stop = 3, length = NUMSPECIES)
    eco = community(uniform_environment(CELLS, AREA, topology = Island()),
                    opts, widths)
    simulate!(eco, 10year, 1month_mean_duration)

    plt = bar(ustrip.(uconvert.(°C, opts)), perspecies(eco),
              xlab = "Temperature preference (°C)", ylab = "Abundance",
              legend = false, size = (1200, 800))
    _savefig(plt, "optima.png")
end

# --- 2. different niche widths, matched ---------------------------------------
# Every species now prefers exactly the landscape's temperature and differs only in how broad a
# range it tolerates. In a constant regime a wide niche buys nothing, so the specialists win.
println("  2. niche widths (matched) ...")
let widths = collect(range(0.0001K, stop = 5.0K, length = NUMSPECIES))
    opts = fill(298.0K, NUMSPECIES)
    eco = community(uniform_environment(CELLS, AREA, topology = Island()),
                    opts, widths)
    simulate!(eco, 10year, 1month_mean_duration)

    plt = bar(ustrip.(widths), perspecies(eco),
              xlab = "Niche width (K)", ylab = "Abundance",
              legend = false, size = (1200, 800))
    _savefig(plt, "widths_matched.png")
end

# --- 3. ...and mismatched -------------------------------------------------------
# The same species against a landscape 1 K off every optimum. Now the narrowest are excluded by the
# mismatch and the widest are still too diffuse, so the advantage moves to intermediate widths - the
# one non-monotonic result in the set, and the reason it is worth plotting beside the last one.
println("  3. niche widths (mismatched) ...")
let widths = collect(range(0.0001K, stop = 5.0K, length = NUMSPECIES))
    opts = fill(298.0K, NUMSPECIES)
    eco = community(uniform_environment(CELLS, AREA, temperature = 299.0K),
                    opts, widths)
    simulate!(eco, 10year, 1month_mean_duration)

    plt = bar(ustrip.(widths), perspecies(eco),
              xlab = "Niche width (K)", ylab = "Abundance",
              legend = false, size = (1200, 800))
    _savefig(plt, "widths_mismatched.png")
end

# --- 4. abundance against the two resources together -------------------------
# **Liebig's law of the minimum, shown as a surface.** Solar is graded down the rows and water
# across the columns of the *same* landscape, so each of the 100 cells is a distinct (solar, water)
# pair and the result is a response surface rather than a line. Abundance should track the **scarcer**
# of the two: rising along the diagonal, and flat along whichever axis is not currently binding.
#
# Grading both at once is the point, not a confound. An earlier rewrite of this experiment graded
# each resource in its own run, on the grounds that varying both "says nothing about either" - that
# was wrong. It would be right for a single cell, where only the minimum is visible; over a grid that
# samples every combination, it is exactly what makes the law legible.
#
# Both ranges start at zero, so the first row and column are uninhabitable. That is deliberate -
# the empty edge is where the surface shows a resource binding absolutely.
println("  4. resource surface (Liebig's law) ...")
let env = uniform_environment(CELLS, AREA, topology = Island()),
    opts = fill(298.0K, NUMSPECIES), widths = fill(2.0K, NUMSPECIES)

    fractions = range(0.0, stop = 1.0, length = CELLS)
    fullsun = env.supply.sunlight.matrix[1, 1]
    fullwater = env.supply.water.matrix[1, 1]
    for (i, f) in enumerate(fractions)
        env.supply.sunlight.matrix[i, :] .= f * fullsun      # solar varies by row (Y)
        env.supply.water.matrix[:, i] .= f * fullwater       # water varies by column (X)
    end

    eco = community(env, opts, widths)
    simulate!(eco, 10year, 1month_mean_duration)

    # `percell` is `(y, x)`, and `heatmap(x, y, z)` wants `z` indexed the same way - so water is
    # the x axis and solar the y, matching how they were written into the layers above.
    plt = heatmap(ustrip.(fractions .* fullwater),
                  ustrip.(fractions .* fullsun),
                  percell(eco),
                  xlab = "Water supply (L/day per cell)",
                  ylab = "Solar supply (kJ/day per cell)",
                  size = (1200, 900))
    _savefig(plt, "resource_surface.png")
end

# --- 5. invariant to grid resolution ------------------------------------------
# The same landscape carved into 1, 4, 25 or 100 cells. Total abundance is a property of the area and
# its resource, so cutting it up more finely must not change the answer.
println("  5. grid resolution ...")
let resolutions = [1, 2, 5, 10]
    totals = map(resolutions) do cells
        eco = community(uniform_environment(cells, AREA,
                                            topology = Island()),
                        fill(298.0K, NUMSPECIES), fill(2.0K, NUMSPECIES))
        simulate!(eco, 10year, 1month_mean_duration)
        return sum(eco.abundances.matrix)
    end
    plt = bar(string.(resolutions .^ 2), totals,
              xlab = "Number of cells", ylab = "Total abundance",
              legend = false, size = (1200, 800))
    _savefig(plt, "resolution.png")
end

# --- 6. abundance scales with area --------------------------------------------
# Resource is given per unit area, so a bigger landscape carries proportionally more of it and must
# support strictly more individuals.
println("  6. area ...")
let areas = [10.0, 20.0, 50.0, 100.0]
    totals = map(areas) do a
        eco = community(uniform_environment(CELLS, a * km^2,
                                            topology = Island()),
                        fill(298.0K, NUMSPECIES), fill(2.0K, NUMSPECIES))
        simulate!(eco, 10year, 1month_mean_duration)
        return sum(eco.abundances.matrix)
    end
    plt = plot(areas, totals, xlab = "Area (km^2)", ylab = "Total abundance",
               legend = false, size = (1200, 800))
    _savefig(plt, "area.png")
end

# --- 7. dispersal distance ----------------------------------------------------
# Two species seeded on opposite edges, with no boundary to wrap around. A wider kernel drains the
# edges and fills the middle, so the edge columns fall as the kernel widens while the centre rises.
println("  7. dispersal ...")
let distances = [0.5, 1.0, 2.0, 4.0] .* km
    profiles = map(distances) do d
        eco = community(uniform_environment(CELLS, AREA,
                                            topology = Island()),
                        fill(298.0K, 2),
                        fill(2.0K, 2), dispersal = d, individuals = 0)
        eco.abundances.grid[1, :, 1] .= 100
        eco.abundances.grid[2, :, end] .= 100
        simulate!(eco, 50year, 1month_mean_duration)
        return vec(sum(percell(eco), dims = 1))
    end
    plt = plot(1:CELLS, profiles, label = reshape(string.(distances), 1, :),
               xlab = "Grid column", ylab = "Abundance", size = (1200, 800))
    _savefig(plt, "dispersal.png")
end

# --- 8. large species pools ---------------------------------------------------
# Heterogeneous species - random widths, optima, dispersal, mortality and body size - so that
# coexistence cannot be an artefact of every species being identical, checked as the pool grows
# tenfold. Every one of these is a plain per-species vector: the builder takes a scalar or a
# length-`n` vector for each, so heterogeneity needs no separate code path.
println("  8. large species pools ...")
let pools = [50, 100, 200, 500]
    fractions = map(pools) do n
        rng = Xoshiro(1)
        widths = rand(rng, Uniform(1.0, 5.0), n) .* K
        opts = 298.0K .+ widths .* range(-3, stop = 3, length = n)
        death = abs.(rand(rng, Normal(0.15, 0.135), n)) ./ year
        eco = community(uniform_environment(CELLS, AREA,
                                            topology = Island()), opts,
                        widths,
                        dispersal = rand(rng, Uniform(0.6, 2.4), n) .* km,
                        birth = death, death = death,
                        sizes = rand(rng, Normal(1.0, 0.05), n) .* m^2)
        simulate!(eco, 10year, 1month_mean_duration)
        return count(>(0), perspecies(eco)) / n
    end
    plt = bar(string.(pools), fractions, ylim = (0, 1),
              xlab = "Species in pool", ylab = "Fraction surviving",
              legend = false, size = (1200, 800))
    _savefig(plt, "pools.png")
end

# --- 9. diversity through time, under each scenario ---------------------------
# **This is where the scenarios live.** They are declarative `Intervention`s, reused verbatim from
# `examples/interventions/`, rather than callbacks mutating the ecosystem from inside the recording
# loop - so the recorder neither knows nor cares that the environment is changing.
#
# Each scenario runs on the ecosystem `examples/interventions/` built for it, because an
# intervention naming a layer (`SetChange(:temperature, ...)`) needs that layer to exist. Comparing
# them against a baseline of the *same* ecosystem with no intervention is the meaningful contrast.
println("  9. diversity through time under each scenario ...")
let scenarios = (baseline = (climate_ecosystem, nothing),
                 warming = (climate_ecosystem, changing_climate()),
                 random_loss = (loss_ecosystem, random_loss(0.05 / year)),
                 clustered_loss = (loss_ecosystem, clustered_loss(0.05 / year)),
                 invasion = (invasion_ecosystem, generalist_invasion()))
    traces = map(scenarios) do (build, intervention)
        return diversity_through_time(build(), times = 10year,
                                      intervention = intervention)
    end

    for measure in (:gamma, :norm_alpha, :norm_beta, :richness, :shannon)
        plt = plot(xlab = "Time (years)", ylab = string(measure),
                   size = (1200, 800))
        for (name, trace) in pairs(traces)
            plot!(plt, ustrip.(uconvert.(year, trace.times)),
                  getproperty(trace, measure), label = string(name))
        end
        _savefig(plt, "through_time_$(measure).png")
    end
end

println("Done - figures are in ", FIGURES)
