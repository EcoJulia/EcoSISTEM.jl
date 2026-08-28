# SPDX-License-Identifier: LGPL-3.0-or-later

# Shared synthetic raster fixtures for the grid/study-area tests. `include`d into each test module
# rather than duplicated, so the two files cannot drift apart. Not named `test_*.jl`, so `runtests.jl`
# neither runs it directly nor expects a matching `src/` file.

# Small synthetic climate rasters (no downloads) for the data-driven builders — a genuine
# `Projected` (CRS-bearing, `Regular`-span) lookup, matching what every real `ClimateRaster` (via
# `datasetread.jl`) carries, since the grid engine/`materialise` etc. reproject via `Rasters.resample`,
# which requires a real CRS + a `Regular` span on both source and target. `crs` defaults to WGS84
# (degree coordinates); pass a projected CRS (with matching length-valued `lat`/`long`) for a
# metric-grid fixture — see `_bngraster`.
#
# **`Intervals(Start)`, because that is what a reader hands back.** Measured 2026-08-13 on real
# reads (`WorldClim{Climate}`, `WorldClim{BioClim}`, `EarthEnv{LandCover}`): `Y`/`X` come back
# `Projected`/**`Intervals(Start)`**/`Regular`, a monthly `Ti` comes back `Sampled`/`Points`/
# `Irregular` holding **1-based** `(1:n) .* month_mean_duration`, and a layer axis comes back
# `Categorical` holding **`Symbol`** codes. `_locus` preserves whatever the source file declares and
# refuses `Points` outright, so nothing shipped ever arrives as `Center`.
#
# **This declared `Center` until 2026-08-13, while the comment above claimed it matched real
# reads.** The difference is half a cell — a `Center` cell labelled 20° spans `[10°, 30°)`, a `Start`
# one spans `[20°, 40°)` — so a query at 19.9° landed in a different cell in the tests than in
# production, ~9 km for WorldClim. That is the mechanism behind four separate hidden bugs
# (`[S3-F3]`, `[S3-F5]`, the `Dim{:var}` band fixture, and the 0-based month axis): **the fixture was
# unrepresentative and the code was calibrated against the fixture.** `locus` stays a keyword so a
# test can still ask for `Center` deliberately, but it must now be asked for.
function _testraster(T, values; lat = (0.0:1.0:4.0) .* °,
                     long = (0.0:1.0:4.0) .* °, crs = Rasters.EPSG(4326),
                     locus = DimensionalData.Lookups.Start(),
                     time = nothing, code = nothing)
    centre = DimensionalData.Lookups.Intervals(locus)
    latstep = length(lat) > 1 ? lat[2] - lat[1] : oneunit(eltype(lat))
    longstep = length(long) > 1 ? long[2] - long[1] : oneunit(eltype(long))
    # **The order follows the values.** It was hardcoded `ForwardOrdered` while the span was
    # computed *signed*, so a **descending** axis — the normal orientation of a north-up GeoTIFF —
    # was declared incoherently: values falling, span negative, order claiming to rise. Nothing
    # noticed while consumers read the raw coordinates and took `abs` of their difference; asking the
    # lookup for its cells' own intervals (`[LOCUS-BLIND]`) is what made the fixture's own
    # declaration load-bearing.
    ord(step) = step < zero(step) ?
                DimensionalData.Lookups.ReverseOrdered() :
                DimensionalData.Lookups.ForwardOrdered()
    yd = Y(Rasters.Projected(collect(lat), sampling = centre, crs = crs,
                             order = ord(latstep),
                             span = DimensionalData.Lookups.Regular(latstep)))
    xd = X(Rasters.Projected(collect(long), sampling = centre, crs = crs,
                             order = ord(longstep),
                             span = DimensionalData.Lookups.Regular(longstep)))
    third = _thirdaxis(time, code)
    return ClimateRaster(T,
                         DimArray(values,
                                  isnothing(third) ? (yd, xd) : (yd, xd, third)))
end

# The optional third axis, built exactly as `datasetread.jl` builds it — `_mkstackaxis`/`_stackcoords`
# for a monthly `Ti`, `_stacklayers`' `Categorical` for a layer axis. Naming which one is wanted is
# how the fixture avoids the ambiguity a bare `1:n` third dimension had: `code = 2` is a *label* on a
# `Categorical` axis and a *coordinate* on a numeric one, and those coincide only when the labels
# happen to be `1:n`.
_thirdaxis(::Nothing, ::Nothing) = nothing
function _thirdaxis(time, code)
    return error("a fixture raster has at most one third axis; got both `time` and `code`.")
end
# `Points`/`Irregular`, which is what a real monthly read carries — *not* `Intervals`. The `Ti`
# half of the `Intervals(Start)` rule is `[TIME-AXIS]`, deliberately not part of this step.
function _thirdaxis(time, ::Nothing)
    return DimensionalData.Ti(DimensionalData.Lookups.Sampled(collect(time),
                                                              order = DimensionalData.Lookups.ForwardOrdered(),
                                                              span = DimensionalData.Lookups.Irregular(),
                                                              sampling = DimensionalData.Lookups.Points()))
end
function _thirdaxis(::Nothing, code)
    return DimensionalData.Dim{:layer}(DimensionalData.Lookups.Categorical(collect(code)))
end

# A projected (British National Grid, EPSG:27700) fixture over central Scotland — eastings
# 245–265 km, northings 640–660 km in 2.5 km cells. The metric counterpart of `_testraster`'s
# default WGS84 grid, for the projected-target paths.
function _bngraster(T, values; east = (245000.0:2500.0:265000.0) .* m,
                    north = (640000.0:2500.0:660000.0) .* m)
    return _testraster(T, values, lat = north, long = east,
                       crs = Rasters.EPSG(27700))
end

# Wrap an already-in-memory `ClimateRaster` as a layer spec.
#
# **A thin alias for the package's own `in_memory_raster`, deliberately** — never a hand-rolled
# `ConstructedSpec(() -> r; axis = axis)`, which is the shape that lets a fixture drift away from
# what the package actually does. Calling the real function means these ~97 uses
# exercise it rather than a parallel copy of it. The short name is kept because it appears that often.
_reg(r; axis = EcoSISTEM.NicheAxis) = EcoSISTEM.in_memory_raster(r, axis = axis)
