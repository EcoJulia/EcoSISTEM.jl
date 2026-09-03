# SPDX-License-Identifier: LGPL-3.0-or-later
#
#### SINGLE SPECIES ####
# One species across Africa, limited by water, on WorldClim bioclim data.
#
# **The grid must be decided before anything is built on it, and projected.** Building the layers by
# hand - reading rasters, rewrapping their coordinates, and handing the result straight to
# `Ecosystem` - decides no grid, so the run stays on WorldClim's own **geographic** (°) one. A degree
# cell's physical size changes with latitude while dispersal assumes a single uniform cell size, so
# the simulator refuses it, correctly.
#
# **This is not cosmetic.** A `StudyArea`
# resolves the CRS, extent and cell size up front; `GridHabitat` only samples the layers onto
# it. That is also what lets `within` do the positioning: WorldClim is global, and an area decided
# from it alone would span the world.

using EcoSISTEM
using EcoSISTEM.Units
using RasterDataSources
using Rasters: EPSG
using Unitful
using Unitful.DefaultSymbols
using Plots

# The gate. `runtests.jl` sets `ECOSISTEM_SCALE=small` for every `Pkg.test` run; run it directly
# (or set `ECOSISTEM_SCALE=large`) for a resolution worth looking at.
const SMALL = get(ENV, "ECOSISTEM_SCALE", "large") == "small"
const CELLSIZE = SMALL ? 200.0km : 50.0km
const YEARS = SMALL ? 2year : 10year

# The two layers, named rather than read: bioclim 1 is annual mean temperature and bioclim 12 annual
# precipitation. Neither the unit nor the axis is written here - the shipped catalogue
# (`data/RasterDataSources/BioClim.csv`) says bio1 is a `Temperature` in °C and bio12 a
# `Precipitation`, and the layer is canonicalised on that authority rather than on a guess.
temperature = SourceSpec(WorldClim{BioClim}, :bio1)
rainfall = SourceSpec(WorldClim{BioClim}, :bio12)

# --- decide the grid, before anything is built on it ------------------------------
# **A projected CRS is required to simulate**, and this is the whole reason the old version of
# this file could not run. `EPSG:10592` (WGS 84 / GLANCE Africa) is the package's *own* advice:
# ask for this extent without a `crs` and the report names it in the warning it raises.
# `within` positions the area and is not optional - WorldClim is a global dataset, so without it
# the grid is the world. `boundingbox` reads the shipped table, so naming a region costs no download.
area = StudyArea(regime = temperature, supply = rainfall,
                 within = EcoSISTEM.boundingbox("Africa", level = "CONTINENT"),
                 crs = EPSG(10592), cellsize = CELLSIZE,
                 verbosity = :silent)

# `GridHabitat` chooses nothing - it samples the named layers onto the grid just decided.
env = GridHabitat(regime = temperature, supply = rainfall, area = area)

# --- one species, tolerant of warmth and limited by water -------------------------
# The tolerance is a `Temperature` and the demand a `Precipitation`, and each says so: meaning
# comes from the declared axis, never from the values or their units.
species = build_species(1, tolerance = (280.0K, 10.0K),
                        toleranceaxis = Temperature,
                        demand = 0.1Unitful.L / day,
                        demandaxis = Precipitation,
                        dispersal = 2CELLSIZE, abundance = 20_000, seed = 1)

eco = build_ecosystem(species, env, seed = 1)

# --- run, recording monthly -------------------------------------------------------
const INTERVAL = 1month_mean_duration
lensim = length((0year):INTERVAL:YEARS)
abuns = zeros(Int64, 1, length(env.active), lensim)
@time simulate_record!(abuns, eco, YEARS, INTERVAL, INTERVAL)

# --- what happened ----------------------------------------------------------------
ny, nx = size(env.active)
frames = reshape(abuns[1, :, :], ny, nx, lensim)
masked(i) = replace(Float64.(frames[:, :, i]), 0.0 => NaN)

println("Africa on a $(ny) × $(nx) grid of $(CELLSIZE) cells ",
        "($(count(env.active)) active) over $(YEARS).")
println("Total abundance: ", sum(frames[:, :, 1]), " -> ",
        sum(frames[:, :, end]), ".")

display(heatmap(masked(1), title = "start", background_color = :lightblue,
                grid = false, color = cgrad(:algae, scale = :exp)))
display(heatmap(masked(lensim), title = "end", background_color = :lightblue,
                grid = false, color = cgrad(:algae, scale = :exp)))
