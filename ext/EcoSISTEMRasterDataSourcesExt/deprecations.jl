# SPDX-License-Identifier: LGPL-3.0-or-later

# --- Deprecations whose signatures name a dataset --------------------------------------------------
#
# Both halves of the package's RDS-typed deprecations — `ClimatePref`'s per-source wrapper
# constructors and readers, and the parent's `*habitat` shims. A `@deprecate` or a method that
# dispatches on `ClimateRaster{WorldClim{…}}` cannot be defined without `RasterDataSources` loaded,
# so the methods are here; the *names* stay exported by the parent, and nothing a user can write
# disappears — it simply has no method until the package that gives it meaning is loaded.

# ---------------------------------------------------------------------------
# Climate types: per-source wrapper constructors → `ClimateRaster{<:RasterDataSource}`
#
# Prior to the unified `read`/`ClimateRaster` API these sources each had their own wrapper type. They are
# retained as deprecated constructors that forward to the equivalent `ClimateRaster` so existing user code
# keeps working (with a warning). The monthly series (`Worldclim_monthly`/`CHELSA_monthly`) are
# `WorldClim{Climate}` / `CHELSA{Climate}`, and enforce no 12-month check of their own: the readers
# that build them already assemble a 12-month stack.
# ---------------------------------------------------------------------------
# These took an `AxisArray` and forwarded it **verbatim** to `ClimateRaster`, which since the
# DimensionalData migration accepts only an `AbstractDimArray` — so all five were dead, failing with a
# bare `MethodError` rather than warning and working, and that was hidden inside the suite's expected
# error count. They now take a `DimArray`: what is deprecated here is the per-source *wrapper type*,
# not the array type it wraps, and `AxisArrays` is on its way out of the package entirely.
# `export_old = false` (the trailing `false`) on each: `@deprecate` cannot export a qualified name,
# and it must not — `ClimatePref` already exports all five, so the macro would be re-exporting them
# from an extension module nobody `using`s.
@deprecate ClimatePref.Worldclim_bioclim(array::DimensionalData.AbstractDimArray) ClimateRaster(RDS.WorldClim{RDS.BioClim},
                                                                                                array) false
@deprecate ClimatePref.CHELSA_bioclim(array::DimensionalData.AbstractDimArray) ClimateRaster(RDS.CHELSA{RDS.BioClim},
                                                                                             array) false
@deprecate ClimatePref.Landcover(array::DimensionalData.AbstractDimArray) ClimateRaster(RDS.EarthEnv{RDS.LandCover},
                                                                                        array) false
@deprecate ClimatePref.Worldclim_monthly(array::DimensionalData.AbstractDimArray) ClimateRaster(RDS.WorldClim{RDS.Climate},
                                                                                                array) false
@deprecate ClimatePref.CHELSA_monthly(array::DimensionalData.AbstractDimArray) ClimateRaster(RDS.CHELSA{RDS.Climate},
                                                                                             array) false

function ClimatePref.readworldclim(T::Type{WorldClim{Climate}}, files;
                                   cut = nothing)
    Base.depwarn("`readworldclim` is deprecated; use `read(WorldClim{Climate}, layers; …)`.",
                 :readworldclim)
    return _readsource(T, _filelist(files); cut = cut)
end

function ClimatePref.readCHELSA_monthly(dir::String, var_name::String;
                                        scale = 1, fn = mean,
                                        cut = nothing)
    Base.depwarn("`readCHELSA_monthly` is deprecated; use `read(CHELSA{Climate}, dir, " *
                 "var_name; scale, fn, cut)`.", :readCHELSA_monthly)
    return read(RDS.CHELSA{RDS.Climate}, dir, var_name; scale = scale, fn = fn,
                cut = cut)
end

# **There are deliberately no data-backed `*habitat` shims here.** Such a shim built a regime
# straight from `raster.array` with no resampling, and `GridHabitat` has no parts-taking constructor
# to assemble one from — a habitat is built from specs plus a `StudyArea`. Their released **names**
# live in `src/deprecations.jl`, raising an error that names both replacements (`SourceSpec(...)` and
# `in_memory_raster(...)`); see `_removedbuilder` there. That needs nothing from this extension,
# which is why this file does not mention them.

# Re-runs the parent's axis loop for the one method whose signature names a `RasterDataSources`
# type; the rest of that loop stays in `src/deprecations.jl`. Keyed on the **axis**, not on a
# supply name: `SolarTimeBudget`/`WaterTimeBudget` are bindings onto
# `Supply{SolarRadiation}`/`Supply{Precipitation}`, so this is the method those bindings resolve to
# for a monthly climate raster.
for axis in (:SolarRadiation, :Precipitation)
    # The `ClimateRaster` forms read their stack from a monthly climate source.
    # **Corrected with the live path, deliberately.** This derives its cell area from the *grid*
    # — unlike `simplehabitat`/`raingradhabitat`, whose area is a caller-supplied total divided
    # evenly and so has no latitude in it — and WorldClim is geographic, so a scalar area was wrong
    # here for the same reason it was wrong everywhere else. The v0.4.0 contract these shims
    # reproduce is *which supply type a unit selects*, not how big a cell is, so correcting the area
    # is a bug fix rather than a change of contract; and its live sibling
    # `Supply{A}(::ClimateRaster{WorldClim{BioClim}})` is corrected, so leaving this one would make
    # two adjacent constructors on the same data source disagree.
    @eval function Supply{$axis}(worldclim::ClimateRaster{WorldClim{Climate}},
                                 time::Integer)
        cellareas = _cellareas(worldclim.array)
        return Supply{$axis}(cancel.(worldclim.array, cellareas, $axis),
                             time)
    end
end
