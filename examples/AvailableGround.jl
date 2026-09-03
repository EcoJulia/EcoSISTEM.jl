# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Three runs over one landscape, each taking one more thing into account than the last, so the
# cost of ignoring land cover can be read straight off the comparison.
#
# **1. Land cover ignored.** Every square metre of every cell counted as growable - what a model
#    that never looked at land cover is implicitly assuming.
#
# **2. Available ground.** Not every square metre of a cell is somewhere a plant can grow. EarthEnv's
#    land cover splits each cell into twelve classes; eight of them are ground something can grow on
#    and four are not (built-up, snow and ice, bare rock, open water). Summing the eight gives the
#    **fraction of each cell that is available**, and that fraction is a *space supply*: species
#    compete for room, and a cell that is half water holds half as many plants.
#
# **3. Sunlight on the green bits too.** Sunlight falling on a car park is not sunlight a plant can use.
#    Multiplying the *same* available fraction by the incident radiation gives the light actually
#    reachable - and the second run shows how much difference that makes.
#
# **The point of the pair.** Run 2's layer composes land-cover layers into another
# land-cover-like layer: same axis in, same axis out. Run 3's takes those same cover fractions and
# produces a **solar radiation** layer - a different axis entirely. Nothing about the inputs says so;
# the `axis = SolarRadiation` on the result is the *only* statement of what the layer means. That is
# the rule the whole layer system runs on: **meaning comes from the declared axis, never from the
# values or their units.**
#
# **The green fraction is deliberately used twice** - once as available room, once to derate the
# light - and that is correct rather than double-counting. A half-green cell really does have half
# the room *and* half the usable light; they are two separate facts about the same cell.
#
# The grid is **projected** (British National Grid), which is what a simulation requires: dispersal
# assumes one uniform cell size, and on a lat/long grid a cell's real area changes with latitude.

using EcoSISTEM
using EcoSISTEM.Units
using RasterDataSources
using Rasters: EPSG
using Unitful, Unitful.DefaultSymbols
using DimensionalData: lookup, Y, X
using Plots

# The gate. `runtests.jl` sets `ECOSISTEM_SCALE=small` for every `Pkg.test` run, which shrinks the
# grid and shortens the run so this stays a cheap integration test; run it directly (or set
# `ECOSISTEM_SCALE=large`) for a resolution worth looking at.
const SMALL = get(ENV, "ECOSISTEM_SCALE", "large") == "small"
# **50 km under the suite, not 25.** This file runs *three* whole simulations to make its
# comparison, so it cost 67 s - 40% of `extras_examples` - at 25 km. Halving the resolution quarters
# the cells, and the point being demonstrated (that ignoring land cover overstates capacity) is a
# ratio between the three runs, which does not need a fine grid to show.
const CELLSIZE = SMALL ? 50.0km : 5.0km
const YEARS = SMALL ? 2year : 20year

# --- the green fraction, computed once and used by both examples ----------------
#
# **Summed, never `1 - Σ(unavailable)`.** The twelve bands do not reliably total 100%: measured
# over all 7.5M EarthEnv cells, 89% sum to exactly 100 and the rest run from 55 to 103. In a cell
# summing to 55 the missing 45 is *unknown* coverage, and subtracting the four unavailable classes
# would silently claim all of it as available.
#
# Classes are named, never numbered - the numeric codes are an EarthEnv implementation detail and
# the shipped catalogue (`data/RasterDataSources/LandCover.csv`) is what maps names to them.
const GROWABLE = (:needleleaf_trees, :evergreen_broadleaf_trees,
                  :deciduous_broadleaf_trees, :other_trees, :shrubs,
                  :herbaceous, :cultivated_and_managed, :regularly_flooded)

# Each band is published as a percentage and canonicalised to a **fraction** on read, so summing
# them needs no scaling: the result is already the share of the cell that is available, 0 to 1.
# Each named class arrives as its **own** `ClimateRaster` - a vector of codes gives one layer per
# code, not one stacked array - so the combine takes them as varargs and simply adds them up: a
# raster behaves like its values and stays a raster, so `sum` is the whole of it and no array type is
# named. The result carries no layer code, and does not need one - the codes disagree, so summing
# drops it, and the `axis =` on the `ConstructedRasterSpec` is what declares what the result means.
available = ConstructedRasterSpec(EarthEnv{LandCover}, collect(GROWABLE),
                                  axis = SurfaceArea) do bands...
    return sum(bands)
end

# Incident sunlight, as a flat areal rate - every cell receives the same radiation per square metre.
# (Real solar data would be a second download; the point here is the composition, not the field.)
const INCIDENT = 5.0e5kJ / (m^2 * day)

# **Demands are set from the cell size, not written as literals, and that is what makes the
# comparison mean anything.** Supplies combine by Liebig's minimum, so a resource only shows up in
# the result when it is the one that *binds*: pick demands too small and all three runs sit far below
# capacity, differ by ~2%, and the example demonstrates nothing. Measured while writing this - the
# first attempt had a per-cell capacity of 140 000 against a population of 1 700.
# Scaling them to the cell keeps that true at either resolution, since `ECOSISTEM_SCALE` changes
# the cell size 5-fold and with it every per-cell supply.
const CELLAREA = uconvert(m^2, CELLSIZE^2)
# **And the two capacities are deliberately DIFFERENT, or the third run shows nothing.** With
# equal capacities, space and light are gated by the same fraction and stay tied - so space binds in
# both runs 2 and 3, gating the light costs exactly nothing, and the last panel duplicates the
# middle one. Measured: that first attempt gave 26.4% for both.
# Setting light *more* abundant than space ungated, but scarcer once gated, is what makes each run
# limited by a different thing - which is the honest picture anyway, and lets the comparison show
# all three.
const SPACEPERCELL = 5_000           # a whole cell of ground could hold this many
const LIGHTPERCELL = 1_200           # its full sunlight, rather fewer
const SPACEDEMAND = CELLAREA / SPACEPERCELL
const LIGHTDEMAND = INCIDENT * CELLAREA / LIGHTPERCELL

# **Example 2's layer, and the reason this file exists.** The inputs are land-cover *cover
# fractions* - dimensionless, on `SurfaceArea`. Multiplying by an incident flux gives a solar flux
# density, and `axis = SolarRadiation` is what declares that. Nothing in the values could have said
# it: a fraction times a constant is just a number until an axis gives it a meaning.
usable_light = ConstructedRasterSpec(EarthEnv{LandCover}, collect(GROWABLE),
                                     axis = SolarRadiation) do bands...
    return sum(bands) .* INCIDENT
end

# --- the landscape ---------------------------------------------------------------
# Mainland Scotland, on the British National Grid. A **projected** CRS is required to simulate:
# dispersal is expressed against one cell size, which only a projected grid actually has.
# **`within` is what positions it, and it is not optional.** EarthEnv land cover is a *global*
# dataset, so an area decided from it alone spans the globe - projected into BNG that gives northings
# from -19 000 km, which is nonsense outside Britain. `extent` is a **size**, not a bounding box:
# it says how big, never where. Measured while writing this - without `within` the grid came out
# 953 × 6 cells over mostly ocean, and every number in the comparison was meaningless.
# `boundingbox` reads the shipped `data/bounding_boxes.csv`, so naming a region costs no download.
area = StudyArea(supply = available,
                 within = EcoSISTEM.boundingbox("Scotland",
                                                coverage = LargestLandmass()),
                 crs = EPSG(27700), cellsize = CELLSIZE, verbosity = :silent)

# The regime and the tolerances are deliberately dull - this example is about the supplies, so
# everything else is a uniform temperature and a niche centred on it.
# They are **synthetic** layers on a **positioned** grid, which the builder refused until step 13
# was written: a generated layer needs shape and orientation, never coordinates, which is exactly why
# the *supply* side had always allowed it.
# **`build_habitat` supplies both**, announcing each as it fills it in - which is the point of
# it, and why `verbosity = :silent` exists for a caller who has read them once.
#
# **Each run is built from the PREVIOUS HABITAT, not from `area`, and that is what makes the
# comparison honest.** Passing a habitat as the first argument takes everything not named from it -
# here the regime and, crucially, **the grid it was actually built on**. So each run inherits its
# predecessor's active cells and can only take more away: the three landscapes are nested by
# construction, and a later, more constrained run can never gain ground an earlier one had already
# dropped. Building all three from `area` instead would let a data layer that costs cells in run 2
# leave run 3 simulating ground run 2 had ruled out - silently, since the grid *shape* is fixed
# either way and only the activity differs. The headline numbers are ratios between the three, so
# that would be a defect in the argument rather than in the code.
# Only the **supply** is named each time; the regime and topology carry across untouched, which is
# exactly the one-thing-changed reading the comparison wants.

# --- example 1: land cover ignored entirely (the baseline) ------------------------
# `SurfaceSpec()` is the whole cell - every square metre counted as growable, which is what a
# model that never looked at land cover is implicitly assuming.
env_plain = build_habitat(verbosity = :silent,
                          supply = (space = SurfaceSpec(),
                                    light = UniformSpec(INCIDENT,
                                                        axis = SolarRadiation)),
                          area = area)

# --- example 2: available ground only, sunlight still full ------------------------
env_space = build_habitat(env_plain, verbosity = :silent,
                          supply = (space = available,
                                    light = UniformSpec(INCIDENT,
                                                        axis = SolarRadiation)))

# --- example 3: available ground, and sunlight gated by it too --------------------
env_both = build_habitat(env_space, verbosity = :silent,
                         supply = (space = available, light = usable_light))

# --- one species set, run on each ------------------------------------------------
# Both demands are per individual: `m^2` of ground and `kJ/day` of light. Supply / demand is a
# count either way - space is a *stock* and light a *flow*, and the model needs only that each axis
# is consistent about which it is.
function community(env)
    species = build_species(DefaultEcosystem(), numspecies = 20,
                            verbosity = :silent,
                            demand = (space = SPACEDEMAND, light = LIGHTDEMAND),
                            demandaxis = (space = SurfaceArea,
                                          light = SolarRadiation),
                            dispersal = 2CELLSIZE, abundance = 500_000,
                            seed = 1)
    eco = build_ecosystem(species, env, seed = 1)
    simulate!(eco, YEARS, 1month_mean_duration)
    return eco
end

plain = community(env_plain)
space = community(env_space)
both = community(env_both)

function percell(eco)
    return dropdims(sum(eco.abundances.dimgrid, dims = :species),
                    dims = :species)
end
total(eco) = sum(eco.abundances.matrix)

greenfrac = sum(env_space.supply.space.matrix) /
            sum(env_plain.supply.space.matrix)
println("Average available (green) fraction of a cell : ",
        round(100 * greenfrac, digits = 1), "%")
println()
println("1. land cover ignored entirely  : ", total(plain))
println("2. available ground only        : ", total(space),
        "   (", round(100 * (1 - total(space) / total(plain)), digits = 1),
        "% fewer)")
println("3. available ground + its light : ", total(both),
        "   (", round(100 * (1 - total(both) / total(plain)), digits = 1),
        "% fewer)")

# **The losses are real but not proportional, and that is the interesting part.** Roughly a third
# of Scotland is not growable ground, yet run 2 does not lose a third of its plants, and gating the
# light as well does not halve them again. Supplies combine by **Liebig's law of the minimum**:
# whichever is scarcest binds, so cutting a resource that was not binding costs nothing and cutting
# one that was costs the difference. A reader expecting the population to track the fraction is
# being taught something wrong - which is exactly why both corrections are worth making explicitly
# rather than assumed to be small.

# --- side by side by side ----------------------------------------------------------
northing = lookup(plain.abundances.dimgrid, Y)
easting = lookup(plain.abundances.dimgrid, X)

maps = [(plain, "1. land cover ignored"), (space, "2. available ground"),
    (both, "3. available ground + light")]
panels = map(maps) do (eco, title)
    return heatmap(easting, northing, Array(percell(eco)),
                   title = title, xlabel = "Easting (BNG)",
                   ylabel = "Northing (BNG)", colorbar_title = "Individuals",
                   aspect_ratio = 1)
end
display(plot(panels..., layout = (1, 3), size = (1500, 450)))

# Where the point actually lands: the loss is **not uniform**. Cells that are mostly built-up or
# water lose nearly everything; cells that were already green barely move. A single headline number
# hides that, and the map does not.
display(heatmap(easting, northing, Array(percell(plain) .- percell(both)),
                title = "Individuals lost by accounting for land cover",
                xlabel = "Easting (BNG)", ylabel = "Northing (BNG)",
                colorbar_title = "Difference", aspect_ratio = 1))
