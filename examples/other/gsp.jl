# SPDX-License-Identifier: LGPL-3.0-or-later
#
# **The same layer in both roles.** `gsp` — CHELSA's growing-season precipitation — appears twice
# in this model, and means something different each time:
#
#   * as a **regime** (a *condition*) it is the water that fell over the growing season, an amount:
#     `L·m⁻²`. Species are matched to it by a tolerance — some like a wet season, some a dry one.
#   * as a **supply** (a *resource*) it is water *per day*, a rate: `L·m⁻²·d⁻¹`. Species consume it
#     against a demand, and compete for it.
#
# Nothing in the script converts between the two. The catalogue records that `gsp` accumulates over
# `gsl` (growing-season length, which varies by cell), and `GridHabitat` desugars the **supply**
# occurrence into `gsp ÷ gsl` automatically — reading, gridding and windowing `gsl` alongside `gsp`,
# and dividing on the *source* grid because division is cell-wise but nonlinear and so does not
# commute with regridding. The regime occurrence is left exactly as it is.
#
# **Why the second layer of each pair is what it is.** Solar radiation is a genuine resource: finite,
# and rival by shading. Temperature is a genuine *condition*: nothing consumes it, and one plant being
# warm does not make another colder — so it is a regime, and the type system enforces that
# (`supplytype(Temperature)` is `nothing`, so it cannot be built as a supply even by mistake).
#
# The first run downloads CHELSA BioClimPlus (`gsp`, `gsl`, `bio1`) and WorldClim solar radiation into
# `EcoSISTEM.assetdir()`; later runs reuse the cache. Everything is seeded, so re-running reproduces
# the same result exactly.

using EcoSISTEM
using EcoSISTEM.Units
using RasterDataSources
using Rasters: EPSG
using Extents: Extent
using Unitful, Unitful.DefaultSymbols
using Distributions
using Random
using DimensionalData
using Plots

# --- the layers ---------------------------------------------------------------
# One spec, used twice below. It is the *role* it is passed in — `regime =` or `supply =` — that
# decides whether it stays an amount or becomes a rate; the spec itself says nothing about which.
gsp = SourceSpec(CHELSA{BioClimPlus}, :gsp)

# Annual mean temperature, from the same CHELSA product so both regime layers share a native grid.
# CHELSA ships it in °C; a regime is always held in an absolute unit, so it arrives as K.
temperature = SourceSpec(CHELSA{BioClimPlus}, :bio1)

# January solar radiation only (`month = 1`), so this is a *static* supply rather than a
# twelve-slice series — it keeps the example about the two roles of `gsp` rather than about
# time-varying layers. Drop the keyword to get all twelve months and a series-driven supply.
sunlight = SourceSpec(WorldClim{Climate}, :srad, month = 1)

# --- the study area: Ireland --------------------------------------------------
# `crs = EPSG(2157)` — Irish Transverse Mercator. A *projected* CRS is required for a physical
# `cellsize`: on a geographic (°) grid a cell's real size shrinks with latitude, and the simulation
# assumes uniform cells.
#
# Only the regime layers are named here, so they decide the grid; `sunlight` is sampled onto it later
# and can mark cells inactive but not move or resize it. `gsl` comes along automatically with the
# `gsp` supply, so it too is read on this grid and within this window.
ireland = Extent(Y = (51.4°, 55.5°), X = (-10.6°, -5.4°))
report = investigate_study_area(regime = (water = gsp,
                                          temperature = temperature),
                                within = ireland, crs = EPSG(2157),
                                cellsize = 10km)
display(report)

# Commit to exactly the grid just displayed, reusing the report's cache so nothing is read twice.
area = StudyArea(report, verbosity = :silent)

# --- the environment ----------------------------------------------------------
# `gsp` appears in **both** arguments. As a regime it stays `L·m⁻²`; as a supply it is silently
# rewritten to `gsp ÷ gsl` and arrives as `L·d⁻¹` per cell. The names line up with the species'
# tolerances and demands below, which is what `build_ecosystem` checks.
habitat = GridHabitat(regime = (water = gsp, temperature = temperature),
                      supply = (water = gsp, solar = sunlight),
                      area = area)

northing = lookup(habitat.regime.water.matrix, Y)
easting = lookup(habitat.regime.water.matrix, X)
active = Array(habitat.active)

# The two readings side by side — the whole point of the example. Same source layer, same cells,
# different quantity: an amount over the season, and a rate per day.
seasonwater = Float64.(ustrip.(Array(habitat.regime.water.matrix)))
seasonwater[.!active] .= NaN
heatmap(easting, northing, seasonwater,
        xlabel = "Easting (ITM)", ylabel = "Northing (ITM)",
        title = "gsp as a REGIME: growing-season water (L/m²)",
        colorbar_title = "L/m²", aspect_ratio = 1)

dailywater = Float64.(ustrip.(Array(habitat.supply.water.matrix)))
dailywater[.!active] .= NaN
heatmap(easting, northing, dailywater,
        xlabel = "Easting (ITM)", ylabel = "Northing (ITM)",
        title = "gsp as a SUPPLY: water per day (gsp ÷ gsl)",
        colorbar_title = "L/day per cell", aspect_ratio = 1)

# --- species ------------------------------------------------------------------
seed = 1
Random.seed!(seed)
numspecies = 20
half = numspecies ÷ 2

# Half prefer a wet growing season, half a dry one; every species has its own temperature optimum
# drawn across the range Ireland actually spans. The two tolerances are named to match the two
# regime layers — `build_ecosystem` refuses a mismatch rather than pairing them by position.
wetlovers = rand(Uniform(700.0, 1100.0), half)
drylovers = rand(Uniform(300.0, 600.0), half)
watermeans = vcat(wetlovers, drylovers)
waterwidths = rand(Uniform(80.0, 200.0), numspecies)

tempmeans = rand(Uniform(281.0, 285.0), numspecies)      # ≈ 8–12 °C, Ireland's real range
tempwidths = rand(Uniform(1.0, 3.0), numspecies)

# One demand per resource, named to match the supply. Water is consumed against the daily rate the
# `gsp ÷ gsl` division produced; solar against January's flux.
waterdemand = rand(Uniform(500.0, 1500.0), numspecies) * Unitful.L / day
solardemand = rand(Uniform(0.5e6, 1.5e6), numspecies) * kJ / day

species = build_species(numspecies,
                        tolerance = (water = (watermeans, waterwidths),
                                     temperature = (tempmeans, tempwidths)),
                        toleranceaxis = (water = EcoSISTEM.GrowingSeasonPrecipitation,
                                         temperature = Temperature),
                        demand = (water = waterdemand, solar = solardemand),
                        demandaxis = (water = Precipitation,
                                      solar = SolarRadiation),
                        dispersal = rand(Uniform(5.0, 20.0), numspecies) * km,
                        abundance = 200_000, seed = seed)

# --- simulate -----------------------------------------------------------------
eco = build_ecosystem(species, habitat, seed = seed)

# The starting distribution, before anything has moved: `populate!` places individuals in
# proportion to supply, so this already reflects where the water and light are.
startgrid = dropdims(sum(eco.abundances.dimgrid, dims = :species),
                     dims = :species)
startmap = Float64.(Array(startgrid))
startmap[.!active] .= NaN
heatmap(easting, northing, startmap,
        xlabel = "Easting (ITM)", ylabel = "Northing (ITM)",
        title = "Total abundance at the start",
        colorbar_title = "individuals", aspect_ratio = 1)

startabun = sum(eco.abundances.matrix)
times = 10year
timestep = 1month_mean_duration
simulate!(eco, times, timestep)
endabun = sum(eco.abundances.matrix)

endgrid = dropdims(sum(eco.abundances.dimgrid, dims = :species),
                   dims = :species)
endmap = Float64.(Array(endgrid))
endmap[.!active] .= NaN
heatmap(easting, northing, endmap,
        xlabel = "Easting (ITM)", ylabel = "Northing (ITM)",
        title = "Total abundance after $times",
        colorbar_title = "individuals", aspect_ratio = 1)

println("Ireland on a $(size(eco.abundances.grid, 2))×$(size(eco.abundances.grid, 3)) grid, " *
        "$(count(active)) active cells.")
println("Total abundance: $startabun -> $endabun after $times.")
println("Species surviving: ",
        count(>(0), sum(eco.abundances.matrix, dims = 2)), " of $numspecies.")
println("Wet-season specialists (1-$half) vs dry ($(half + 1)-$numspecies): ",
        sum(sum(eco.abundances.matrix, dims = 2)[1:half]), " vs ",
        sum(sum(eco.abundances.matrix, dims = 2)[(half + 1):numspecies]))
