# SPDX-License-Identifier: LGPL-3.0-or-later
#
#### SINGLE SPECIES ####
# One species across Africa that lives on vegetated land cover and is limited by water.
#
# **The grid is decided before anything is built on it**, for the same reason as its sibling
# `Worldclim_africa.jl`: layers built by hand decide no grid, leaving the run on the sources' own
# geographic (°) coordinates, which cannot be simulated.
#
# **And a second constraint, which only becomes visible once the first is met**: pairing an EarthEnv
# land-cover grid with a WorldClim layer *upresolutioned ×2* pairs two grids that are genuinely
# different grids - the builder refused them with *"`supply`'s grid does not match the regime's
# grid"*. The hand-rolled `upresolution` existed precisely to force the two onto a common grid, and
# a `StudyArea` is the thing that does that properly: both layers are now sampled onto one decided
# grid, so no resampling has to be arranged by hand and nothing can silently disagree.

using EcoSISTEM
using EcoSISTEM.Units
using RasterDataSources
using Rasters: EPSG
using Unitful
using Unitful.DefaultSymbols
using Plots

# The gate - see `Worldclim_africa.jl`.
const SMALL = get(ENV, "ECOSISTEM_SCALE", "large") == "small"
const CELLSIZE = SMALL ? 200.0km : 50.0km
const YEARS = SMALL ? 2year : 10year

# --- the regime: which land-cover class wins each cell -----------------------------
# `compress_landcover` collapses EarthEnv's twelve per-class cover fractions into the single
# **winning class code** per cell. A bare dataset argument (no code) hands the combine the whole
# multi-band raster, which is exactly what it needs.
#
# **`LandCoverTypology`, and the axis is what makes the layer categorical** - not the fact that the
# values happen to be integers. Declaring `SurfaceArea` here (a *continuous* axis) canonicalises them
# to `Float64` and the builder then refuses the pairing outright: *"layer `tolerance` is Int64 in the
# species tolerance but Float64 in the environment regime"*. That refusal is the axes rule working:
# meaning comes from the declared axis, never from the values, so the axis has to say "these are class
# labels" for a `SimpleCategoricalTolerance` to match them.
landcover = ConstructedSpec(compress_landcover, EarthEnv{LandCover},
                            axis = LandCoverTypology)

# Bioclim 13 is precipitation of the wettest month. Its unit and axis come from the shipped
# catalogue rather than being written here.
rainfall = SourceSpec(WorldClim{BioClim}, :bio13)

# --- decide the grid, before anything is built on it ------------------------------
# A **projected** CRS is required to simulate. `EPSG:10592` (WGS 84 / GLANCE Africa) is the
# package's own advice for this extent. `within` positions the area: both sources are global.
area = StudyArea(regime = landcover, supply = rainfall,
                 within = EcoSISTEM.boundingbox("Africa", level = "CONTINENT"),
                 crs = EPSG(10592), cellsize = CELLSIZE,
                 verbosity = :silent)

env = GridHabitat(regime = landcover, supply = rainfall, area = area)

# --- one species, living on the vegetated classes ---------------------------------
# **Classes by name, never by number.** The numeric codes are an EarthEnv implementation detail;
# the shipped catalogue (`data/RasterDataSources/LandCover.csv`) maps names to them, and
# `EcoSISTEM.landcoverclass` is the lookup (`public`, not exported). The old version wrote `collect(1:8)`, which silently depended on
# those eight codes staying in that order.
const LIVEABLE = (:needleleaf_trees, :evergreen_broadleaf_trees,
                  :deciduous_broadleaf_trees, :other_trees, :shrubs,
                  :herbaceous, :cultivated_and_managed, :regularly_flooded)

# A **pre-built** tolerance, which `build_species` accepts as it stands: a land-cover tolerance is
# a set of acceptable classes, not a mean and a width, so there is no `(mean, width)` form for it and
# no `toleranceaxis` to give - the tolerance already carries its own meaning. It does still declare
# its `axis`, which is what pairs it with the `LandCoverTypology` regime above.
tolerance = SimpleCategoricalTolerance([collect(EcoSISTEM.landcoverclass.(LIVEABLE))],
                                       axis = LandCoverTypology)

species = build_species(1, tolerance = tolerance,
                        demand = 10.0Unitful.L / day,
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

display(heatmap(Array(env.regime.matrix), title = "winning land-cover class",
                grid = false))
display(heatmap(masked(lensim), title = "end", background_color = :lightblue,
                grid = false, color = cgrad(:algae, scale = :exp)))
