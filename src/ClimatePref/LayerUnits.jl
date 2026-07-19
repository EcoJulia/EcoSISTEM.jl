# SPDX-License-Identifier: LGPL-3.0-or-later

using CSV
using Unitful
using EcoSISTEM.Units
using EcoSISTEM: NicheAxis
using RasterDataSources
using InteractiveUtils
# `const RDS = RasterDataSources` is already defined for this module in ClimateTypes.jl.

# The dataset subtype wrapped by a `RasterDataSource`, e.g. `WorldClim{BioClim}` →
# `BioClim` and `EarthEnv{LandCover}` → `LandCover`. Sources without a dataset
# parameter (such as `AWAP`) are returned unchanged.
function _datasettype(::Type{T}) where {T <: RDS.RasterDataSource}
    params = Base.unwrap_unionall(T).parameters
    return isempty(params) ? T : first(params)
end

# The layer table shipped in the package `data/RasterDataSources/` directory, named by convention
# after the dataset type (`WorldClim{BioClim}` → `data/RasterDataSources/BioClim.csv`,
# `EarthEnv{LandCover}` → `data/RasterDataSources/LandCover.csv`); the file's presence is what
# makes a source supported. (They live under `RasterDataSources/` to keep them distinct from other
# shipped data such as `data/bounding_boxes.csv`.)
function _layerfile(T::Type)
    path = pkgdir(@__MODULE__, "data", "RasterDataSources",
                  "$(nameof(_datasettype(T))).csv")
    isfile(path) ||
        return error("No layer table for raster source `$T` (expected $(basename(path)) in data/)")
    return path
end

# Cache of parsed `Code` → (units, axis) strings, keyed by file path. Blank cells are stored
# as `""` and interpreted lazily — a blank `Units` ⇒ dimensionless, a blank/absent `Axis` ⇒
# unclassified — so a malformed entry can only affect the layer that actually uses it.
const _LAYER_CACHE = Dict{String,
                          Dict{String,
                               @NamedTuple{units::String, axis::String}}}()

function _layertable(path::String)
    return get!(_LAYER_CACHE, path) do
        table = CSV.File(path; normalizenames = true)
        cols = propertynames(table)
        function cell(row, col)
            return (col in cols && !ismissing(getproperty(row, col))) ?
                   String(strip(String(getproperty(row, col)))) : ""
        end
        rows = Dict{String, @NamedTuple{units::String, axis::String}}()
        for row in table
            ismissing(row.Code) && continue
            rows[string(row.Code)] = (units = cell(row, :Units),
                                      axis = cell(row, :Axis))
        end
        return rows
    end
end

# All concrete `NicheAxis` types. `subtypes` only returns direct children, so recurse
# through the abstract intermediates (e.g. `AbstractTemperature`) to reach the leaf axes.
function _leafaxes(T = NicheAxis)
    return isabstracttype(T) ?
           mapreduce(_leafaxes, vcat, subtypes(T); init = Type[]) : [T]
end

# Resolve an axis name from a table to its `NicheAxis` type by autodiscovery — no registry:
# any loaded concrete `NicheAxis` with that name works (build-time only).
function _resolve_axis(name::AbstractString)
    matches = filter(A -> string(nameof(A)) == name, _leafaxes())
    length(matches) == 1 ||
        return error("`$name` does not name exactly one loaded `NicheAxis` (found $(length(matches)))")
    return only(matches)
end

function _layerrow(T::Type{<:RDS.RasterDataSource}, code)
    tbl = _layertable(_layerfile(T))
    key = string(code)
    haskey(tbl, key) ||
        return error("Layer `$code` is not in the table for raster source `$T`")
    return tbl[key]
end

"""
    layerunit(T::Type{<:RasterDataSources.RasterDataSource}, code)

Return the physical unit of layer `code` in raster dataset `T` (e.g.
`layerunit(WorldClim{BioClim}, 1)` is `K`, `layerunit(WorldClim{Climate}, :srad)` is
`kJ m⁻² day⁻¹`). Looked up in the dataset's shipped `data/` table and parsed with
`Unitful.uparse`; a **blank** `Units` cell means **dimensionless** (`NoUnits`). `code` is
matched by its string form, so integer layer numbers and `Symbol`/`String` keys both work.
"""
function layerunit(T::Type{<:RDS.RasterDataSource}, code)
    u = _layerrow(T, code).units
    return isempty(u) ? NoUnits : uparse(u, unit_context = [Unitful, Units])
end

"""
    layeraxis(T::Type{<:RasterDataSources.RasterDataSource}, code)

Return the [`NicheAxis`](@ref) type declared for layer `code` of dataset `T` in its shipped
`data/` table, or `nothing` if the `Axis` cell is blank (or absent). A blank axis is an
*unclassified* layer — documented and unit-bearing, but not modelled as a niche axis; to
use it, pass an explicit axis when building the layer.
"""
function layeraxis(T::Type{<:RDS.RasterDataSource}, code)
    a = _layerrow(T, code).axis
    return isempty(a) ? nothing : _resolve_axis(a)
end
