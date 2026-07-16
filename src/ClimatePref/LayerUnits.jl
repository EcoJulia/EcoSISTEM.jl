# SPDX-License-Identifier: LGPL-3.0-or-later

using CSV
using Unitful
using EcoSISTEM.Units
using RasterDataSources
# `const RDS = RasterDataSources` is already defined for this module in ClimateTypes.jl.

# The dataset subtype wrapped by a `RasterDataSource`, e.g. `WorldClim{BioClim}` →
# `BioClim` and `EarthEnv{LandCover}` → `LandCover`. Sources without a dataset
# parameter (such as `AWAP`) are returned unchanged.
function _datasettype(::Type{T}) where {T <: RDS.RasterDataSource}
    params = Base.unwrap_unionall(T).parameters
    return isempty(params) ? T : first(params)
end

# The layer→unit table shipped alongside this file for each RasterDataSources dataset
# type. `LandCover` maps to `Landcover.csv` (note the case). Topography also has a
# table but no RasterDataSources type, so it is not reachable through this map.
const _LAYERUNIT_FILES = Dict{Type, String}(RDS.BioClim => "BioClim.csv",
                                            RDS.BioClimPlus => "BioClimPlus.csv",
                                            RDS.Climate => "Climate.csv",
                                            RDS.LandCover => "Landcover.csv",
                                            RDS.Elevation => "Elevation.csv",
                                            RDS.HabitatHeterogeneity => "HabitatHeterogeneity.csv")

# Cache of parsed `Code`→unit-string maps, keyed by CSV basename. The unit strings are
# `uparse`d lazily (in `layerunit`) so a malformed entry in one table can only affect
# the layer that actually uses it.
const _LAYERUNIT_CACHE = Dict{String, Dict{String, String}}()

function _layerunit_file(T::Type)
    ds = _datasettype(T)
    haskey(_LAYERUNIT_FILES, ds) ||
        return error("No layer-unit table is known for raster source `$T`")
    return _LAYERUNIT_FILES[ds]
end

function _layerunit_strings(file::String)
    return get!(_LAYERUNIT_CACHE, file) do
        table = CSV.File(joinpath(@__DIR__, file); normalizenames = true)
        strings = Dict{String, String}()
        for row in table
            (ismissing(row.Code) || ismissing(row.Units)) && continue
            strings[string(row.Code)] = String(strip(String(row.Units)))
        end
        return strings
    end
end

"""
    layerunit(T::Type{<:RasterDataSources.RasterDataSource}, code)

Return the physical unit of the layer identified by `code` in the raster dataset
`T` (e.g. `layerunit(WorldClim{BioClim}, 1)` is `K`, `layerunit(WorldClim{Climate},
:srad)` is `kJ m⁻² day⁻¹`). The unit is looked up in the layer table shipped for `T`
in this directory and parsed with `Unitful.uparse` (so `NoUnits` becomes the
dimensionless unit). `code` is matched by its string form, so integer layer numbers
and `Symbol`/`String` layer keys are both accepted.
"""
function layerunit(T::Type{<:RDS.RasterDataSource}, code)
    strings = _layerunit_strings(_layerunit_file(T))
    key = string(code)
    haskey(strings, key) ||
        return error("Layer `$code` is not in the unit table for raster source `$T`")
    return uparse(strings[key], unit_context = [Unitful, Units])
end
