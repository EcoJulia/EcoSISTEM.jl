# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Download EarthEnv land-cover layer 9 with EcoSISTEM, cut it to the bounding box of
# mainland Scotland, and plot it. Then use a real shapefile (NatureScot's "Landscape Map of
# Scotland") to restrict a simulation's active area to Scotland's actual coastline — rather
# than its rectangular bounding box — and demonstrate the effect on the land-cover heatmap.

using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using RasterDataSources
using Unitful, Unitful.DefaultSymbols
using Downloads: download
using Logging
using Plots

# --- mainland-Scotland bounding box, from the shipped data/bounding_boxes.csv table ----------
# `boundingbox` returns a `(lat, long)` NamedTuple of ° intervals ready to pass as `cut`/`region`.
# Pass `islands = true` for the island-inclusive extent, or e.g. `round = 5°` to snap outwards.
scotland = EcoSISTEM.ClimatePref.boundingbox("Scotland")

# --- download + read EarthEnv land-cover layer 9, cut to Scotland ----------
landcover = read(EarthEnv{LandCover}, 9; cut = scotland)

# --- plot ------------------------------------------------------------------
# `landcover.array` is a lat × long `AxisArray` (dim 1 = latitude, dim 2 = longitude).
lat = ustrip.(landcover.array.axes[1].val)
long = ustrip.(landcover.array.axes[2].val)
heatmap(long, lat, Array(landcover.array);
        xlabel = "Longitude (°)", ylabel = "Latitude (°)",
        title = "EarthEnv land cover class 9 — mainland Scotland",
        colorbar_title = "% cover", aspect_ratio = 1)

# --- shapefile-derived active area ------------------------------------------
# NatureScot's "Landscape Map of Scotland" (79 landscape character features, British National
# Grid). Downloaded once into EcoSISTEM's own Scratch.jl cache (`assetdir`) — re-runs (and a
# warm `~/.julia/scratchspaces` cache in CI) skip the download entirely.
shapezip = joinpath(EcoSISTEM.assetdir(), "LSCMAP_SCOTLAND_SHP_27700.zip")
isfile(shapezip) ||
    download("https://gis-downloads.nature.scot/LSCMAP_SCOTLAND_SHP_27700.zip",
             shapezip)
scotlandshape = shapemask(shapezip)

# The mask restricted to Scotland's actual coastline (not the rectangular bounding box) —
# resolved onto the same grid as the land-cover raster above. Scotland's bounding box spans
# enough latitude for `build_environment` to log an informational note on meridian convergence
# (`@info`, not a warning); silenced locally so this stays usable as a `@test_nowarn` example —
# the notebooks harness already carries the same carve-out for this exact message.
landcover_env = with_logger(NullLogger()) do
    return build_environment(landcover; active = scotlandshape)
end

# Re-plot the land-cover heatmap with everything outside the shapefile's coastline blanked out
# (`NaN`) — the plot visibly changes from a rectangle to Scotland's actual outline, showing the
# effect of `active = shapemask(...)`.
masked = Float64.(Array(landcover.array))
masked[.!landcover_env.active] .= NaN
heatmap(long, lat, masked;
        xlabel = "Longitude (°)", ylabel = "Latitude (°)",
        title = "Land cover restricted to the shapefile active area",
        colorbar_title = "% cover", aspect_ratio = 1)
