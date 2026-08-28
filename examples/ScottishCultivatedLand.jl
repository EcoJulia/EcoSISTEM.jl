# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A real-data ecosystem over the whole of Scotland, built entirely from lazy `*Spec` objects —
# no manually-read/instantiated raster or shapefile is ever passed around directly. The regime is
# EarthEnv land cover's own `cultivated_and_managed` class (its raw per-cell % cover, 0-100),
# the active area is Scotland's real coastline (a shapefile, not its rectangular bounding box),
# and the resource is a synthetic north-south solar gradient laid over the real grid. Species are
# split between "lovers" (do best where % cultivated cover is high) and "haters" (do best where
# it's low), with random solar demand and random initial abundances, then simulated for a few
# months. Two heatmaps (raw, then masked to Scotland's coastline) are plotted directly from the
# built `habitat`'s own fields — no separate raw read is needed to see the land-cover data or the
# effect of the shapefile mask.
#
# Everything is seeded, so re-running this script reproduces exactly the same result. The first
# run downloads real data (a NatureScot shapefile, an EarthEnv land-cover raster) into
# `EcoSISTEM.assetdir()`; later runs reuse the cache.

using EcoSISTEM
using EcoSISTEM.Units
using RasterDataSources
using Rasters: EPSG
using Unitful, Unitful.DefaultSymbols
using Distributions
using Random
using DimensionalData: lookup, Y, X
using Plots

# The regime: EarthEnv land cover class 7 (`cultivated_and_managed`, see
# `data/RasterDataSources/LandCover.csv`) read as a lazy `SourceSpec` — its own per-cell % cover,
# not the full multi-class winning-code collapse `landcoverhabitat` uses elsewhere (EarthEnv's own
# per-class bands are continuous % cover, not categorical — `iscategorical` treats them as such).
cultivated = SourceSpec(EarthEnv{LandCover}, :cultivated_and_managed)

# NatureScot's "Landscape Map of Scotland" (79 landscape character features, British National
# Grid). `ShapeSpec` downloads and caches URLs lazily (via `EcoSISTEM.CachedAsset`, under its
# own `EcoSISTEM.assetdir(owner = ShapeSpec)` subdirectory) the first time it is
# materialised, so re-runs skip the download.
active_shape = ShapeSpec("https://gis-downloads.nature.scot/LSCMAP_SCOTLAND_SHP_27700.zip")

# The synthetic resource: a solar supply that increases heading south (real solar intensity
# increases toward the equator) — `orientation = 180°` (south) puts the higher value at the
# southern edge, the lower one at the northern edge.
sunlight = GradientSpec(3000.0kJ / km^2 / day, 8000.0kJ / km^2 / day,
                        axis = SolarRadiation, orientation = 180°)

# --- decide the grid, and look at it before committing -------------------------
# The grid is chosen once, up front, and can be argued with before anything is built on it.
# `investigate_study_area` runs exactly the analysis `StudyArea` would and hands back a report
# instead of an area, so the CRS, cell size, extent, per-layer cost and warnings are all visible
# first — nothing about the grid is decided inside `GridHabitat`.
#
#   * `within = active_shape` — the coastline both restricts the simulated cells *and sets the
#     extent*, so no separate bounding box is needed: the mask leads, and the grid comes out
#     Scotland-shaped rather than global. (A shipped bounding box such as
#     `boundingbox("Scotland")` would say nothing the coastline's own envelope does not.)
#   * `crs = EPSG(27700)` — British National Grid, Scotland's own national grid and the CRS the
#     LSCMAP shapefile is published in, so mask and grid agree without reprojection. A projected
#     CRS is *required* for a physical `cellsize`: on a geographic (°) grid a cell's real size
#     varies with latitude, and simulation assumes uniform cells.
#   * `cellsize = 5km` — square 5 km cells, exactly.
#
# The report warns that EarthEnv's WGS84 land cover has to be resampled onto that BNG grid (no
# layer here is natively in BNG, so none can be kept exactly) — real information, worth seeing.
report = investigate_study_area(regime = cultivated, within = active_shape,
                                crs = EPSG(27700), cellsize = 5km)
display(report)

# Commit to exactly the grid just displayed. Passing the report back as the `base` reuses its cache
# of raw reads, so the (global) land-cover layer is not read a second time; `:silent` because the
# report has already said everything there is to say about this area.
area = StudyArea(report, verbosity = :silent)

# Now build on it. `GridHabitat` chooses nothing — it samples the layers onto the decided grid.
habitat = GridHabitat(regime = cultivated, supply = sunlight,
                      area = area)

# --- heatmaps: raw cultivated-land % cover, then masked to Scotland's real coastline ----------
# Both plotted directly from the already-built `habitat`'s own fields (`regime.matrix`, `active`)
# — no separate raw `read(...)` is needed, since `GridHabitat` already materialised the
# layer onto a real `(Y, X)`-dimensioned grid. Those coordinates are British National Grid
# eastings/northings in metres (not degrees), because of the `crs = EPSG(27700)` above.
northing = lookup(habitat.regime.matrix, Y)
easting = lookup(habitat.regime.matrix, X)
heatmap(easting, northing, Array(habitat.regime.matrix),
        xlabel = "Easting (BNG)", ylabel = "Northing (BNG)",
        title = "EarthEnv land cover class 7 (cultivated_and_managed) — Scotland",
        colorbar_title = "% cover", aspect_ratio = 1)

# The same data with everything outside the shapefile's coastline blanked out (`NaN`) — the plot
# visibly changes from a rectangle to Scotland's actual outline, showing the effect of
# `active = active_shape` above.
masked = Float64.(Array(habitat.regime.matrix))
masked[.!Array(habitat.active)] .= NaN
heatmap(easting, northing, masked,
        xlabel = "Easting (BNG)", ylabel = "Northing (BNG)",
        title = "Cultivated land restricted to Scotland's real coastline",
        colorbar_title = "% cover", aspect_ratio = 1)

# --- 20 species, half loving cultivated land, half hating it -------------------
seed = 1
Random.seed!(seed)

numspecies = 20
halfspecies = numspecies ÷ 2

# "Lovers" have a niche mean near 100% cultivated cover; "haters" near 0%. Both groups get a
# random width, so the two groups aren't otherwise identical.
#
# **Cover is a fraction, 0–1, not a percentage.** The shipped land-cover layers are published on
# the 0–100 scale — which their `Units = percent` cell records — and that is divided out once when
# the layer is built, so a tolerance written here must be in the same 0–1 frame the regime is in.
# Getting this wrong does not error: a mean of 70 against values in 0–1 simply puts every species
# far outside its niche, and the run dies quietly.
niche_means = vcat(rand(Uniform(0.7, 1.0), halfspecies),
                   rand(Uniform(0.0, 0.3), halfspecies))
niche_widths = rand(Uniform(0.1, 0.25), numspecies)

# Random solar-radiation demand and random dispersal distance per species.
demand = rand(Uniform(2.0, 10.0), numspecies) * kJ / day
dispersal = rand(Uniform(3.0, 15.0), numspecies) * km

species = build_species(numspecies, tolerance = (niche_means, niche_widths),
                        toleranceaxis = EcoSISTEM.SurfaceArea,
                        dispersal = dispersal,
                        demand = demand,
                        demandaxis = SolarRadiation,
                        # A total of 20000 individuals, split at random across species —
                        # genuinely random initial abundances, reproducible via `seed`.
                        abundance = 20_000, seed = seed)

# --- assemble and simulate ----------------------------------------------------
eco = build_ecosystem(species, habitat, seed = seed)

times = 6month_mean_duration
timestep = 1month_mean_duration

startabun = sum(eco.abundances.matrix)
println("Per-species initial abundance: ",
        sum(eco.abundances.matrix, dims = 2)[:, 1])

simulate!(eco, times, timestep)
endabun = sum(eco.abundances.matrix)

println("Simulated $numspecies species over a $(size(eco.abundances.grid, 2)) × " *
        "$(size(eco.abundances.grid, 3)) grid ($(count(habitat.active)) active cells " *
        "within Scotland's coastline) for $times.")
println("Total abundance: $startabun -> $endabun.")
println("Per-species final abundance (first $halfspecies = lovers, " *
        "last $halfspecies = haters): ",
        sum(eco.abundances.matrix, dims = 2)[:, 1])

# --- final spatial distribution: lovers vs haters ------------------------------
# `eco.abundances.grid` is a `(species, Y, X)` `DimArray`; selecting each half by name and
# summing over species leaves a `(Y, X)` `DimArray` — still real coordinates, so `heatmap` plots
# it directly, just like the habitat heatmaps above.
lovers_total = dropdims(sum(eco.abundances.grid[species = 1:halfspecies],
                            dims = :species), dims = :species)
haters_total = dropdims(sum(eco.abundances.grid[species = (halfspecies + 1):numspecies],
                            dims = :species), dims = :species)
heatmap(easting, northing, Array(lovers_total),
        xlabel = "Easting (BNG)", ylabel = "Northing (BNG)",
        title = "Final population: cultivated-land \"lovers\"",
        colorbar_title = "Total abundance", aspect_ratio = 1)
heatmap(easting, northing, Array(haters_total),
        xlabel = "Easting (BNG)", ylabel = "Northing (BNG)",
        title = "Final population: cultivated-land \"haters\"",
        colorbar_title = "Total abundance", aspect_ratio = 1)
