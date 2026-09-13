# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Figures 8 and 9: a circular "continent" 8,000 km across on a synthetic grid of 80 km cells, at a
# uniform 25 C with a uniform solar supply, where a specialist is introduced into a population of
# generalists after a burn-in.
#
#   - Figure 8, `InvasionCircle.pdf`: 50,000 species, one of them the specialist -
#     per-cell normalised alpha diversity at q = 0 at the end of the burn-in and 50 and 100 years
#     after the introduction, and at q = 1 at 100 years, both from Diversity.jl. The per-cell maps are saved to `continent.jld2` and the figure is drawn from them,
#     so it can be redrawn without rerunning.
#   - Figure 9, `Invasion.pdf`: a generalist alone, then one specialist of each of six niche widths
#     introduced at the same cell, and how fast each spreads.
#
# Both the seeding and the introduction are `AddAbundance` interventions on a one-cell mask, and
# the maps are taken from the abundance matrix at each target year.
#
# Included by `examples/paper.jl`, never run directly.

using Random
using JLD2
using Plots
using Diversity

# The published continent's 50,000 species, and a 20-cell version of it for the suite. Species and
# years shrink with it; the cell size, the temperature and the supply do not. `numspecies` counts
# the specialist. Supply is 20 MJ/m^2/day, a typical daily insolation for Africa, and a plant
# demands 5 MJ/day, a quarter of a square metre's worth, so a cell carries 2.6e10 plants and the
# disc 2e14 - see `continent_capacity`, which is where the generalists start. From there the
# population settles to its equilibrium within twenty years, so a 25-year burn-in is enough.
function continent_configuration()
    return paper_small() ?
           (side = 20, radius = 800.0km, numspecies = 10,
            burnin = 1year, after = 2year,
            widths = [0.5, 5.0, 50.0] .* K) :
           (side = 100, radius = 4000.0km, numspecies = 50_000,
            burnin = 25year, after = 100year,
            widths = [0.5, 1.0, 5.0, 10.0, 25.0, 50.0] .* K)
end

const PAPER_CELLSIZE = 80.0km
const PAPER_GENERALIST = 50.0K
const PAPER_DEMAND = 5.0MJ / day
# How often the continent run reports, in simulated years: every year of the suite's short run,
# every fifth year of the real one.
const PAPER_REPORT_YEARS = paper_small() ? 1 : 5
# In kelvin, because it is also every species' optimum, and a tolerance's parameter vectors must
# share one unit with its widths.
const PAPER_TEMPERATURE = uconvert(K, 25.0°C)

# The circular continent: uniform temperature and supply, `Island`.
function continent_environment(continent)
    extent = continent.side * PAPER_CELLSIZE
    area = StudyArea(extent = (extent, extent), cellsize = PAPER_CELLSIZE,
                     within = CircleMaskSpec(radius = continent.radius),
                     verbosity = :silent)
    return GridHabitat(regime = UniformSpec(PAPER_TEMPERATURE,
                                            axis = Temperature),
                       supply = UniformSpec(20.0MJ / m^2 / day,
                                            axis = SolarRadiation),
                       area = area, topology = Island())
end

# `n` species of the given niche widths about the continent's own temperature, with the abundances
# given per species, so an invader can start absent. `BirthOnlyMovement` with a 15 km kernel;
# birth and death 0.6/year and survival 0.1 as the original.
function continent_species(widths, abundances)
    n = length(widths)
    return build_species(n, tolerance = (fill(PAPER_TEMPERATURE, n), widths),
                         toleranceaxis = Temperature,
                         demand = PAPER_DEMAND,
                         demandaxis = SolarRadiation, dispersal = 15.0km,
                         movement = BirthOnlyMovement, birth = 0.6 / year,
                         death = 0.6 / year, survival = 0.1,
                         abundance = abundances, seed = PAPER_SEED)
end

# The plants the continent's light can support: every active cell's supply over what one plant
# demands. It is where the generalists start, and it is not where they stay - births balance deaths
# below it, at a fraction of capacity set by the species' match to the environment - so the burn-in
# descends to the equilibrium rather than climbing to it.
function continent_capacity(env)
    supply = EcoSISTEM.totalsupply(env).SolarRadiation
    return round(Int, uconvert(NoUnits, supply / PAPER_DEMAND))
end

# One active cell, drawn from the seeded stream, as a `(Y, X)` mask - the same cell for every run.
function origin_mask(active)
    mask = falses(size(active))
    mask[rand(Xoshiro(PAPER_SEED), findall(active))] = true
    return mask
end

# Per-cell normalised alpha diversity at q = 0 (the species present) and q = 1 (the effective number
# at Shannon's weighting) from a `species x cells` matrix, by Diversity.jl, laid out on the `(Y, X)`
# grid with inactive cells `NaN`. Only the active cells are handed over, as proportions of the whole
# landscape, so an empty inactive column neither warns nor distorts.
function _percellmaps(abundances, active)
    cells = findall(vec(active))
    meta = Metacommunity(abundances[:, cells] ./ sum(abundances))
    map(q -> begin
            m = fill(NaN, size(active))
            m[cells] .= norm_sub_alpha(meta, q)[!, :diversity]
            return m
        end, (alpha0 = 0, alpha1 = 1))
end

# Figure 8's run: the generalists with the specialist present at zero, the burn-in, the
# introduction one step after the first map is taken, and maps at 50 and 100 years after it. The
# action runs every step and takes a map when the clock reaches each target, which keeps the
# recording on the year rather than one step past it.
function continent_run(continent, env)
    n = continent.numspecies
    widths = vcat(fill(PAPER_GENERALIST, n - 1), 0.5K)
    abundances = vcat(fill(div(continent_capacity(env), n - 1), n - 1), 0)
    eco = build_ecosystem(continent_species(widths, abundances), env,
                          seed = PAPER_SEED)
    introduction = continent.burnin + PAPER_TIMESTEP
    invasion = Intervention(AtTime(introduction),
                            CellMask(origin_mask(env.active)),
                            AddAbundance(n, 100))
    targets = [continent.burnin, continent.burnin + continent.after / 2,
        continent.burnin + continent.after]
    maps = Dict{Int, Any}()
    started = time()
    total = ustrip(year, last(targets))
    reported = 0
    simulate_action!(eco, last(targets), PAPER_TIMESTEP, PAPER_TIMESTEP,
                     intervention = invasion) do _
        now = EcoSISTEM.simulationtime(eco)
        # A line every `PAPER_REPORT_YEARS`: the year, the projected finish, and how many
        # individuals the generalists and the specialist hold (Diversity's per-species
        # metaabundance), so a run of hours can be watched.
        elapsed = ustrip(year, uconvert(year, now))
        if floor(Int, elapsed + 1e-6) >= reported + PAPER_REPORT_YEARS
            reported = floor(Int, elapsed + 1e-6)
            counts = getmetaabundance(eco, true)
            hours = (time() - started) / 3600
            togo = hours * (total / elapsed - 1)
            println("    year ", reported, " of ", round(Int, total), ": ",
                    round(hours, digits = 2), " h elapsed, about ",
                    round(togo, digits = 1), " h to go (",
                    Libc.strftime("%H:%M", time() + 3600togo),
                    "); generalists ", sum(counts[1:(n - 1)]),
                    ", specialist ", counts[n], ", total ", sum(counts))
            # Redirected to a file, stdout is buffered until exit, and a line hours late is no
            # use.
            flush(stdout)
        end
        i = findfirst(i -> !haskey(maps, i) && now >= targets[i],
                      eachindex(targets))
        isnothing(i) && return nothing
        maps[i] = _percellmaps(Array(eco.abundances.matrix), env.active)
        println("    continent at ",
                round(typeof(1.0year), uconvert(year, now), digits = 1),
                " after ", round(time() - started, digits = 1), " s")
        flush(stdout)
        return nothing
    end
    return (burnin = maps[1], halfway = maps[2], final = maps[3])
end

# Figure 8 from the saved maps: q = 0 at the three times, q = 1 at the last, where the invasive's
# dominance by abundance reaches well beyond the core it has cleared of generalists.
function continent_figure(maps, nspecies)
    panel(m, title;
          clim) = heatmap(m, c = :algae, clim = clim,
                          aspect_ratio = 1,
                          background_color_inside = :lightblue,
                          ticks = false, title = title,
                          titleloc = :left)
    return plot(panel(maps.burnin.alpha0, "A", clim = (0, nspecies)),
                panel(maps.halfway.alpha0, "B", clim = (0, nspecies)),
                panel(maps.final.alpha0, "C", clim = (0, nspecies)),
                panel(maps.final.alpha1, "D", clim = (0, nspecies)),
                layout = @layout([a b; c d]), size = (1200, 1200),
                margin = 4Plots.mm, titlefontsize = 16)
end

# Figure 9's run for one specialist width: the generalist seeded at the origin on the first step,
# the specialist added there after the burn-in, and the specialist's mean distance from the origin
# over the cells it occupies at the end, per year since its introduction.
function invasion_speed(continent, env, width)
    mask = origin_mask(env.active)
    eco = build_ecosystem(continent_species([PAPER_GENERALIST, width], [0, 0]),
                          env, seed = PAPER_SEED)
    arrivals = InterventionSet(Intervention(AtTime(0year), CellMask(mask),
                                            AddAbundance(1, 100)),
                               Intervention(AtTime(continent.burnin),
                                            CellMask(mask),
                                            AddAbundance(2, 100)))
    total = continent.burnin + continent.after
    simulate_action!(_ -> nothing, eco, total, total, PAPER_TIMESTEP,
                     intervention = arrivals)
    final = Array(eco.abundances.matrix)
    origin = Tuple(findfirst(mask))
    occupied = findall(>(0), reshape(final[2, :], size(env.active)))
    distances = [hypot((Tuple(c) .- origin)...) * PAPER_CELLSIZE
                 for c in occupied]
    return mean(distances) / continent.after
end

timed("InvasionCircle") do
    continent = continent_configuration()
    env = continent_environment(continent)
    println("    ", count(env.active), " active cells of ",
            length(env.active), ", ", continent.numspecies, " species, ",
            continent_capacity(env), " plants at capacity")
    maps = continent_run(continent, env)
    path = joinpath(PAPER_DIR, "continent.jld2")
    # A plain matrix, so the file reads back without the package and its unit types loaded.
    jldsave(path; maps = maps, active = Matrix(env.active),
            numspecies = continent.numspecies)
    println("    wrote ", path)
    return savefigure(continent_figure(maps, continent.numspecies),
                      "InvasionCircle")
end

timed("Invasion") do
    continent = continent_configuration()
    env = continent_environment(continent)
    speeds = map(continent.widths) do width
        speed = invasion_speed(continent, env, width)
        println("    specialist of width ", width, ": ",
                round(typeof(1.0km / year), speed, digits = 1))
        return speed
    end
    advantage = ustrip.(K, PAPER_GENERALIST .- continent.widths)
    plt = plot(advantage, ustrip.(km / year, speeds), marker = :circle,
               xlab = "Selective advantage (C)",
               ylab = "Average invasion speed (km/year)", label = "",
               size = (700, 500), margin = 6Plots.mm)
    return savefigure(plt, "Invasion")
end
