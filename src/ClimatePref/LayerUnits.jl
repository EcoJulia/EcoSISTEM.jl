# SPDX-License-Identifier: LGPL-3.0-or-later

using CSV
using Unitful
using EcoSISTEM.Units
using EcoSISTEM: NicheAxis, canonicalunit
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
            meta = (units = cell(row, :Units), axis = cell(row, :Axis))
            # `Code` may be a semicolon-separated list of equivalent RDS aliases for one layer (e.g. `1;bio1`,
            # or `Contrast;contrast`); register each alias as a key so any accepted form resolves.
            for code in split(string(row.Code), ";")
                key = strip(code)
                isempty(key) && continue
                rows[String(key)] = meta
            end
        end
        return rows
    end
end

# All concrete `NicheAxis` types. `subtypes` only returns direct children, so recurse
# through the abstract intermediates (e.g. `TemperatureAxis`) to reach the leaf axes.
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

# ---------------------------------------------------------------------------
# Catalogue + discovery helpers (public, not exported — declared in ClimatePref.jl)
# ---------------------------------------------------------------------------

"""
    LayerRecord

A single layer's catalogue entry: its `dataset` (the table basename, e.g. `:BioClim`), code `aliases` (the
`;`-separated `Code`, e.g. `["1", "bio1"]`), `name`, `definition`, physical `unit` (`NoUnits` if
dimensionless), `axis` (a [`NicheAxis`](@ref) type, or `nothing` if unclassified), the `sources` that provide
it (the full list of supporting `RasterDataSources`), the `numlayers` the code expands to, and its `temporal`
native reporting period — a genuine `Unitful.Units` time unit (`day`/`month`/`quarter`/`year`), or `nothing`
for a static/non-time-varying value. Returned by [`layerinfo`](@ref) / [`layersbyaxis`](@ref).
"""
struct LayerRecord
    dataset::Symbol
    aliases::Vector{String}
    name::String
    definition::String
    unit::Unitful.Units
    axis::Union{DataType, Nothing}
    sources::Vector{String}
    numlayers::Int
    temporal::Union{Unitful.Units, Nothing}
end

# Full per-layer metadata for the discovery helpers, read ONCE from every shipped table (the lean
# `_layertable` keeps only units/axis, keyed per alias). One record per CSV row; the `;`-separated `Code`
# becomes `aliases`. Cached in `_CATALOGUE`.
const _CATALOGUE = LayerRecord[]

# Parse a shipped CSV's `Temporal Resolution` cell (day/month/quarter/year, or blank for "no native
# period" — a static/non-time-varying value) into the unit itself, or `nothing`. Validated as a genuine
# time unit at parse time, since a non-time value here would silently break any later use as a native
# reporting period.
function _parsetemporal(cell)
    isempty(cell) && return nothing
    u = uparse(cell, unit_context = [Unitful, Units])
    dimension(u) == Unitful.𝐓 ||
        return error("Temporal Resolution `$cell` is not a time unit (got dimension $(dimension(u)))")
    return u
end

function _catalogue()
    isempty(_CATALOGUE) || return _CATALOGUE
    datadir = pkgdir(@__MODULE__, "data", "RasterDataSources")
    for f in sort(filter(endswith(".csv"), readdir(datadir)))
        dataset = Symbol(first(splitext(f)))
        table = CSV.File(joinpath(datadir, f); normalizenames = true)
        cols = propertynames(table)
        # `string` (not `String`) so numeric columns like `NumLayers` (an Int) convert too
        cell(row, col) = (col in cols && !ismissing(getproperty(row, col))) ?
                         String(strip(string(getproperty(row, col)))) : ""
        splitsemis(s) = filter(!isempty, String.(strip.(split(s, ";"))))
        for row in table
            ismissing(row.Code) && continue
            aliases = splitsemis(string(row.Code))
            isempty(aliases) && continue
            u = cell(row, :Units)
            a = cell(row, :Axis)
            nl = cell(row, :NumLayers)
            push!(_CATALOGUE,
                  LayerRecord(dataset, aliases, cell(row, :Name),
                              cell(row, :Definition),
                              isempty(u) ? NoUnits :
                              uparse(u, unit_context = [Unitful, Units]),
                              isempty(a) ? nothing : _resolve_axis(a),
                              splitsemis(cell(row, :Source)),
                              isempty(nl) ? 1 : parse(Int, nl),
                              _parsetemporal(cell(row, :Temporal_Resolution))))
        end
    end
    return _CATALOGUE
end

# The abstract-supertype chain of a leaf axis, up to (but excluding) `NicheAxis` — e.g.
# `Temperature` → `[TemperatureAxis]`; a direct leaf like `SolarRadiation` → `[]`.
function _axischain(A::DataType)
    chain = DataType[]
    S = supertype(A)
    while S !== NicheAxis && S !== Any
        push!(chain, S)
        S = supertype(S)
    end
    return chain
end

_unitstr(u) = u == NoUnits ? "dimensionless" : string(u)

# compact one-liner — used when a `LayerRecord` is an element of a printed `Vector` (e.g. `layersbyaxis`)
function Base.show(io::IO, r::LayerRecord)
    return print(io, join(r.aliases, ";"), " — ", r.name, " [",
                 _unitstr(r.unit),
                 ", ", r.dataset, "]")
end

# detailed block — used at the REPL for a single `LayerRecord` (e.g. `layerinfo(T, code)`)
function Base.show(io::IO, ::MIME"text/plain", r::LayerRecord)
    println(io, r.name)
    println(io, "  dataset    : ", r.dataset)
    println(io, "  code(s)    : ", join(r.aliases, ", "))
    println(io, "  sources    : ",
            isempty(r.sources) ? "—" : join(r.sources, ", "))
    println(io, "  unit       : ", _unitstr(r.unit))
    if r.axis === nothing
        println(io, "  axis       : unclassified")
    else
        chain = _axischain(r.axis)
        print(io, "  axis       : ", nameof(r.axis))
        isempty(chain) || print(io, " (⊂ ", join(nameof.(chain), " ⊂ "), ")")
        cu = canonicalunit(r.axis())
        cu === nothing || print(io, "  [canonical unit ", _unitstr(cu), "]")
        println(io)
    end
    if r.temporal !== nothing || r.numlayers > 1
        print(io, "  temporal   : ",
              r.temporal === nothing ? "—" : string(r.temporal))
        r.numlayers > 1 &&
            print(io, " (", r.numlayers, " layers) — select with `month=`")
        println(io)
    end
    return print(io, "  definition : ", r.definition)
end

"""
    layerinfo(T::Type{<:RasterDataSources.RasterDataSource}, code)
    layerinfo(code)

Return the catalogue [`LayerRecord`](@ref) for a layer, with its name, definition, unit, niche axis and
sources. The two-argument form looks up `code` in dataset `T`'s table (mirrors `layerunit(T, code)`) and
returns one record. The single-argument form (an `Integer`/`Symbol`/`String` code) searches EVERY shipped
table and returns a `Vector` of all matches — the same code can appear in several datasets (e.g. `bio1` is in
both `BioClim` and `BioClimPlus`).
"""
function layerinfo(T::Type{<:RDS.RasterDataSource}, code)
    ds = nameof(_datasettype(T))
    key = string(code)
    for r in _catalogue()
        (r.dataset === ds && key in r.aliases) && return r
    end
    return error("Layer `$code` is not in the `$ds` table.")
end

function layerinfo(code::Union{Integer, Symbol, AbstractString})
    key = string(code)
    recs = filter(r -> key in r.aliases, _catalogue())
    isempty(recs) &&
        return error("No layer with code `$code` in any shipped table.")
    return recs
end

"""
    layersbyaxis(A::Type{<:NicheAxis})

Return a `Vector` of [`LayerRecord`](@ref)s for every layer (across all shipped tables) whose axis is `A` or a
concrete leaf beneath it. Passing an abstract group spans all its axes — `layersbyaxis(TemperatureAxis)`
returns every temperature layer, `layersbyaxis(Temperature)` just that one axis. Printed as a compact list
of codes with short summaries.
"""
function layersbyaxis(A::Type{<:NicheAxis})
    leaves = Set{DataType}(_leafaxes(A))
    return filter(r -> r.axis !== nothing && r.axis in leaves, _catalogue())
end

"""
    AxisNode

One node of the [`layeraxes`](@ref) niche-axis tree: the `axis` type itself, the `names` of shipped layers
that use it directly (an abstract grouping axis never carries a layer itself — only a concrete leaf does, so
this is empty for every non-leaf node), and its `children` — the axis types immediately below it, each an
`AxisNode` built the same way, recursively down to the concrete leaves.
"""
struct AxisNode
    axis::Type{<:NicheAxis}
    names::Vector{String}
    children::Vector{AxisNode}
end

# Shipped layer names (the CSV `Name` column, via `LayerRecord.name`) that use axis `A` directly.
function _axisnames(A::Type{<:NicheAxis})
    return [r.name for r in _catalogue() if r.axis === A]
end

# Build one `AxisNode` for `A`, recursing through `subtypes` down to the concrete leaves.
function _axisnode(A::Type{<:NicheAxis})
    children = sort(_axisnode.(subtypes(A)); by = c -> string(nameof(c.axis)))
    return AxisNode(A, _axisnames(A), children)
end

function _showtree(io::IO, node::AxisNode, depth::Int)
    println(io, "  "^depth, nameof(node.axis),
            isempty(node.names) ? "" : " — " * join(node.names, ", "))
    for child in node.children
        _showtree(io, child, depth + 1)
    end
end

function Base.show(io::IO, node::AxisNode)
    return print(io, nameof(node.axis), " (", length(node.children),
                 " children, ",
                 length(node.names), " layers)")
end

Base.show(io::IO, ::MIME"text/plain", node::AxisNode) = _showtree(io, node, 0)

"""
    layeraxes(A::Type{<:NicheAxis} = NicheAxis)

Return the niche-axis hierarchy at and below `A` (the whole tree by default) as a nested [`AxisNode`](@ref):
each node carries the shipped layer names that use its axis directly (only a concrete leaf ever has any —
an abstract grouping node's own `names` is always empty) plus its child axis nodes, recursively down to the
concrete leaves. Use it to discover which axes exist and what shipped layers use them, then drill into a
group with [`layersbyaxis`](@ref).
"""
layeraxes(A::Type{<:NicheAxis} = NicheAxis) = _axisnode(A)
