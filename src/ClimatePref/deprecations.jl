# SPDX-License-Identifier: LGPL-3.0-or-later

# ===========================================================================
# Deprecations — `ClimatePref` submodule
#
# (One live deprecation cannot be moved here: the `xmin`/`xmax`/`ymin`/`ymax`
# keyword branch inside `readfile(file; …)` in `ReadData.jl` is part of that
# live method's body.)
# ===========================================================================

# ---------------------------------------------------------------------------
# Climate types: per-source wrapper constructors → `ClimateRaster{<:RasterDataSource}`
#
# Prior to the unified `read`/`ClimateRaster` API these sources each had their own wrapper type. They are
# retained as deprecated constructors that forward to the equivalent `ClimateRaster` so existing user code
# keeps working (with a warning). The monthly series (`Worldclim_monthly`/`CHELSA_monthly`) are
# `WorldClim{Climate}` / `CHELSA{Climate}`; the 12-month check they used to enforce is dropped (the readers
# that build them already assemble a 12-month stack).
# ---------------------------------------------------------------------------
@deprecate Worldclim_bioclim(array::AxisArray) ClimateRaster(RDS.WorldClim{RDS.BioClim},
                                                             array)
@deprecate CHELSA_bioclim(array::AxisArray) ClimateRaster(RDS.CHELSA{RDS.BioClim},
                                                          array)
@deprecate Landcover(array::AxisArray) ClimateRaster(RDS.EarthEnv{RDS.LandCover},
                                                     array)
@deprecate Worldclim_monthly(array::AxisArray) ClimateRaster(RDS.WorldClim{RDS.Climate},
                                                             array)
@deprecate CHELSA_monthly(array::AxisArray) ClimateRaster(RDS.CHELSA{RDS.Climate},
                                                          array)

# ---------------------------------------------------------------------------
# Readers → the unified `read`/`cut` API
# ---------------------------------------------------------------------------
"""
    readworldclim(T::Type{WorldClim{Climate}}, files; cut = nothing)

Deprecated — use `read(WorldClim{Climate}, layers; …)` instead, which downloads via `getraster`
and reads through the same machinery. Retained to read an already-resolved set of monthly climate
raster file paths.
"""
function readworldclim(T::Type{WorldClim{Climate}}, files; cut = nothing)
    Base.depwarn("`readworldclim` is deprecated; use `read(WorldClim{Climate}, layers; …)`.",
                 :readworldclim)
    return _readsource(T, _filelist(files); cut = cut)
end

# Deprecated positional-extent form `readfile(file, xmin, xmax, ymin, ymax)`; convert to `cut` and forward
# to the keyword `readfile(file; cut)` method in `ReadData.jl`.
function readfile(file::String, xmin, xmax, ymin, ymax)
    Base.depwarn("`readfile(file, xmin, xmax, ymin, ymax)` is deprecated; pass " *
                 "`cut = LatLong(ymin .. ymax, xmin .. xmax)` (e.g. from `boundingbox`) " *
                 "instead.", :readfile)
    return readfile(file; cut = LatLong(ymin .. ymax, xmin .. xmax))
end

# ---------------------------------------------------------------------------
# `compressLC` → `compressLandCover` (v0.4.0 rename; the `LC` land-cover abbreviation is expanded, matching
# RasterDataSources' `LandCover`).
# ---------------------------------------------------------------------------
@deprecate compressLC compressLandCover
