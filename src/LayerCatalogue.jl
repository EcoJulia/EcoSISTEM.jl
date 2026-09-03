# SPDX-License-Identifier: LGPL-3.0-or-later
#
# One row of the shipped layer catalogue, and a node of the axis tree the discovery helpers print.

using Unitful

using CSV

using .Units

using InteractiveUtils

"""
    LayerRecord

A single layer's catalogue entry: its `dataset` (the table basename, e.g. `:BioClim`), code `aliases` (the
`;`-separated `Code`, e.g. `["1", "bio1"]`), `name`, `definition`, physical `unit` (`NoUnits` if
dimensionless), `axis` (a [`NicheAxis`](@ref) type, or `nothing` if unclassified), the `sources` that provide
it (the full list of supporting `RasterDataSources`), the `numslices` the code expands to, its `temporal`
native reporting period - a genuine `Unitful.Units` time unit (`day`/`month_mean_duration`/`year`), or `nothing`
for a static/non-time-varying value - its `valuetype`, and its `category`. Returned by
[`layerinfo`](@ref) / [`layersbyaxis`](@ref).

`category` says what *kind* of quantity the layer is - `:instantaneous`, `:range`, `:rate`, `:count`,
`:stock`, `:balance` or `:categorical` - which is a different question from `valuetype` and not
derivable from it. It is the only machine-readable record of the rate-vs-level distinction, the one
`canonicalunit(::Type{<:Precipitation})` got wrong.

`publishedscale` is the factor a **provider** ships the layer at when that is not the documented
one - `100` for `bio4` and `10 000` for the fourteen HabitatHeterogeneity layers, blank everywhere else.
It is divided out on read, so values arrive at the order of magnitude the source's own documentation
claims. It lives in the catalogue rather than in code deliberately: the correction is then visible
to `layerinfo` as documentation *and* applied as behaviour, instead of being buried in a reader.
The test guarding it is **inverted** - it fires when the upstream data is *repaired*, since a test
asserting the data is still broken is what tells us the workaround can go.

The remaining fields carry the rest of the shipped table so that **no column is dead data**:
`officialunit` (the source's own documented unit string, against which `unit` is our Unitful
reading), `notes`, and the file-specific `group` (BioClimPlus) and `order` (HabitatHeterogeneity),
each `nothing` where the table has no such column. `UnitDimension` is deliberately **validated
rather than stored**, because it duplicates `dimension(unit)` and could only drift.

Two columns have been **removed from the shipped tables outright**, on the same reasoning: a
column that cannot disagree with what we already know is not a check, it is a second place to get
wrong. `IsRate` said exactly `category === :rate` in all 139 rows. LandCover's
`Filename (full version)`/`(reduced version)` held `consensus_full_class_N.tif` /
`Consensus_reduced_class_N.tif` - precisely what RasterDataSources' own `rastername` computes from
the class number, down to its odd lower/upper-case asymmetry, so we were duplicating another
package's convention and would not have noticed it changing.

`period` is the interval the value accumulated over - an [`AbstractAccumulationPeriod`](@ref), or
`nothing` where none applies. It is a **different question** from `temporal`: `temporal` is how often
the layer is *sampled*, twelve monthly slices, while `period` is what each value *accumulated over*.
For `prec` those are both months and the difference does not show; for `srad` they differ, sampled
monthly and accumulated per day; and for a single-slice layer such as `bio13` there is no sampling at
all.

`valuetype` is `:categorical`, `:discrete` or `:continuous`, and the distinction that matters is
**categorical against the other two**, not the names' everyday senses: a `:discrete` layer (a
day-of-year, or a count of growing-degree days) is an ordinary number that may be averaged, while a
`:categorical` one is a class code - Köppen-Geiger and friends, `kg0`-`kg5`, the only six in the
shipped tables - where the arithmetic mean of two classes is meaningless. Anything resampling a layer
must therefore pick nearest-class (`:mode`) for `:categorical` and may interpolate the rest.

So `:discrete` is **descriptive metadata, not a switch**: every consumer in the package tests
`=== :categorical`, and nothing anywhere reads `:discrete`. It is kept because it is honest about
the source data - `gsl` really does count whole days - but a three-way column is a two-way split in
practice, and a layer marked `:discrete` becomes an ordinary `ContinuousLayer` that is interpolated
like any other number. Do not add behaviour keyed on it without deciding that question first.
"""
struct LayerRecord
    dataset::Symbol
    aliases::Vector{String}
    name::String
    definition::String
    unit::Unitful.Units
    axis::Union{DataType, Nothing}
    sources::Vector{String}
    numslices::Int
    temporal::Union{Unitful.Units, Nothing}
    period::Union{AbstractAccumulationPeriod, Nothing}
    valuetype::Symbol
    category::Symbol
    officialunit::String
    notes::String
    group::Union{String, Nothing}
    order::Union{String, Nothing}
    # Published values are this multiple of the layer's **documented** unit, or `nothing` where
    # they agree (the overwhelming majority). Recorded per layer because it is a property of that
    # provider's data, and discoverable through `layerinfo` rather than buried in code.
    publishedscale::Union{Rational{Int}, Nothing}
end

"""
    AxisNode

One node of the [`layeraxes`](@ref) niche-axis tree: the `axis` type itself, the `names` of shipped layers
that use it directly (an abstract grouping axis never carries a layer itself - only a concrete leaf does, so
this is empty for every non-leaf node), and its `children` - the axis types immediately below it, each an
`AxisNode` built the same way, recursively down to the concrete leaves.
"""
struct AxisNode
    axis::Type{<:NicheAxis}
    names::Vector{String}
    children::Vector{AxisNode}
end

# Cache of parsed `Code` -> (units, axis) strings, keyed by file path. Blank cells are stored
# as `""` and interpreted lazily - a blank `Units` => dimensionless, a blank/absent `Axis` =>
# unclassified - so a malformed entry can only affect the layer that actually uses it.
const _LAYER_CACHE = Dict{String,
                          Dict{String,
                               @NamedTuple{units::String, axis::String}}}()

# Name -> axis type, memoised.
#
# **Without this cache a catalogue build costs 83 seconds.** `_resolveaxis` runs once per row, 139 of
# them, and each call walks `_leafaxes()`, which is `subtypes()` all the way down. `subtypes` scans
# every loaded module, so it is around 0.6 s a call rather than free, and every axis being abstract
# means the walk recurses through some fifty types rather than stopping at fifteen concrete leaves.
# Nothing fails without the cache; the gates simply get slower, which is why such a cost can go
# unnoticed.
const _AXIS_BY_NAME = Dict{String, Type}()

# ---------------------------------------------------------------------------
# Catalogue + discovery helpers (public, not exported - declared in ClimatePref.jl)
# ---------------------------------------------------------------------------

# The closed vocabulary of `perslice=` right-hand sides. Only one today; named rather than assumed so
# that a typo is an error instead of silently becoming the sole option.
const _PERSLICE_PERIODS = ("calendar_month",)

# The cross-column invariant, and the reason ~100 blank cells mean something: a period is
# *required* exactly where a quantity accumulates, and *forbidden* everywhere else. Without it a
# blank is indistinguishable from "nobody filled this in" - and it was this check that exposed `gsl`
# being miscatalogued as a `stock` when it is a span between two dates.
const _PERIOD_CATEGORIES = (:rate, :stock, :balance)

# `:range` is deliberately in *neither* list: a range **inherits** whether it has a period from the
# quantity it ranges over. The range of an instantaneous value (`bio2`, a diurnal temperature range)
# has none; the range of an accumulated one (`cmi_range`, the difference of two monthly totals) has
# its siblings'. So a per-row check cannot decide it, and `_checkrangerows` below does it against the
# siblings instead - a stronger check than the blanket "must be blank" this replaces.
#
# Allowing a period on a `range` cannot make one supply-eligible: §4b moved supply-eligibility to
# the axis's `Resource`-role `canonicalunit`, and no `*Range` axis declares one. That leak was the
# whole reason for forbidding it, and it was closed elsewhere.
const _PERIOD_INHERITED_CATEGORIES = (:range,)

# Full per-layer metadata for the discovery helpers, read ONCE from every shipped table (the lean
# `_layertable` keeps only units/axis, keyed per alias). One record per CSV row; the `;`-separated `Code`
# becomes `aliases`. Cached in `_CATALOGUE`.
const _CATALOGUE = LayerRecord[]

# The `ValueType` column as a checked `Symbol`. Errors rather than defaulting on an unrecognised
# value: a typo silently becoming "safe to interpolate" would mangle a class-code layer without
# anyone noticing, which is exactly the failure this column exists to prevent. A blank cell is a
# genuine "not stated" and stays `:continuous`, matching how every other optional column reads.
const _VALUETYPES = (:categorical, :discrete, :continuous)

# The `Category` column as a checked `Symbol` - what *kind* of quantity a layer holds, which is a
# different question from `ValueType` (how its values behave) and not derivable from it: `count` and
# one `stock` are `:discrete`, the other five `stock`s and every `rate`/`instantaneous`/`range`/
# `balance` are `:continuous`. It is the only machine-readable record of the rate-vs-level
# distinction, which is exactly what `canonicalunit(::Type{<:Precipitation})` got wrong.
const _CATEGORIES = (:instantaneous, :range, :rate, :count, :stock, :balance,
                     :categorical)

# **There is no separate "is a rate" column, deliberately**: `Category === :rate` already says it, and
# two columns holding one fact is how they drift apart. Nothing infers rate-ness from the **unit**
# either - `WindSpeed` is `m s^-1`, dimensionally `𝐋𝐓^-1`, and yet honestly `instantaneous` - so the
# dimensional check below keys on `dimension` and never on whether something is a rate.

# The two notations for a dimension reduced to one comparable form. Needed because the table and
# Unitful write the same thing differently - the CSV has `𝐋 * 𝐓^-1`, Unitful shows `𝐋 𝐓^-1` - so a
# string comparison reports every dimensioned layer as a mismatch. Superscripts become `^n`, the
# separators are normalised, and the factors are sorted so ordering cannot matter either.
const _SUPERSCRIPTS = Dict('⁻' => '-', '⁰' => '0', '¹' => '1', '²' => '2',
                           '³' => '3', '⁴' => '4', '⁵' => '5', '⁶' => '6',
                           '⁷' => '7', '⁸' => '8', '⁹' => '9')

# The exact column set every shipped table must have, plus the extra columns each file is allowed.
#
# The guard that stops a column going quietly dead. The readers fetch cells through a permissive
# helper that answers `""` for an absent column, so a renamed or misspelled header degrades to blank
# rather than failing, and a whole column can end up shipped but unread with nothing to show for it.
# Extra columns are refused too: an unexpected header is either a typo for a real one, or data nothing
# will ever read.
const _REQUIRED_COLUMNS = (:Code, :Axis, :OfficialUnit, :Units, :UnitDimension,
                           :Category, :ValueType, :Name, :Definition,
                           :Source, :NumSlices, :Temporal_Resolution,
                           :AccumulationPeriod, :Notes)

# `PublishedScaleFactor` records a measured discrepancy between what a provider *publishes* and
# what it *documents* - see `LayerRecord.publishedscale`. It is optional and lives only in the tables
# that need it, so an empty cell is not an invitation to guess: it means nothing is known.
# Add it to any table where such a defect is found, and record the evidence in `Notes`.
const _OPTIONAL_COLUMNS = Dict(:BioClimPlus => (:Group,),
                               :HabitatHeterogeneity => (:Order,
                                :PublishedScaleFactor),
                               :BioClim => (:PublishedScaleFactor,))

# == The shipped layer catalogue ====================================================================
# compact one-liner - used when a `LayerRecord` is an element of a printed `Vector` (e.g. `layersbyaxis`)
function Base.show(io::IO, r::LayerRecord)
    return print(io, join(r.aliases, ";"), " - ", r.name, " [",
                 _unitstr(r.unit),
                 ", ", r.dataset, "]")
end

# detailed block - used at the REPL for a single `LayerRecord` (e.g. `layerinfo(T, code)`)
function Base.show(io::IO, ::MIME"text/plain", r::LayerRecord)
    println(io, r.name)
    println(io, "  dataset    : ", r.dataset)
    println(io, "  code(s)    : ", join(r.aliases, ", "))
    println(io, "  sources    : ",
            isempty(r.sources) ? "-" : join(r.sources, ", "))
    println(io, "  unit       : ", _unitstr(r.unit))
    if isnothing(r.axis)
        println(io, "  axis       : unclassified")
    else
        chain = _axischain(r.axis)
        print(io, "  axis       : ", nameof(r.axis))
        isempty(chain) || print(io, " (⊂ ", join(nameof.(chain), " ⊂ "), ")")
        cu = canonicalunit(r.axis)
        isnothing(cu) || print(io, "  [canonical unit ", _unitstr(cu), "]")
        println(io)
    end
    # Sampling and accumulation are shown on separate lines because they are separate questions:
    # `srad` is *sampled* monthly but *accumulates* per day.
    if !isnothing(r.temporal) || r.numslices > 1
        print(io, "  sampled    : ", string(something(r.temporal, "-")))
        r.numslices > 1 &&
            print(io, " (", r.numslices, " slices) - select with `month=`")
        println(io)
    end
    isnothing(r.period) ||
        println(io, "  accumulated: ", _periodphrase(r.period))
    return print(io, "  definition : ", r.definition)
end

function Base.show(io::IO, node::AxisNode)
    return print(io, nameof(node.axis), " (", length(node.children),
                 " children, ",
                 length(node.names), " layers)")
end

Base.show(io::IO, ::MIME"text/plain", node::AxisNode) = _showtree(io, node, 0)

"""
    layerunit(T::Type, code)

Return the physical unit of layer `code` in raster dataset `T` (e.g.
`layerunit(WorldClim{BioClim}, 1)` is `K`, `layerunit(WorldClim{Climate}, :srad)` is
`kJ m^-2 day^-1`). Looked up in the dataset's shipped `data/` table and parsed with
`Unitful.uparse`; a **blank** `Units` cell means **dimensionless** (`NoUnits`). `code` is
matched by its string form, so integer layer numbers and `Symbol`/`String` keys both work.

# Arguments

  - `T`: the `RasterDataSources` dataset type, e.g. `WorldClim{BioClim}`.
  - `code`: which layer of it - an `Integer` layer number, or a `Symbol`/`String` alias
    (`:bio1`, `"bio1"`, `1` all reach the same row).
"""
function layerunit(T::Type, code)
    u = _layerrow(T, code).units
    return isempty(u) ? NoUnits : uparse(u, unit_context = [Unitful, Units])
end

"""
    layeraxis(T::Type, code)

Return the [`NicheAxis`](@ref) type declared for layer `code` of dataset `T` in its shipped
`data/` table, or `nothing` if the `Axis` cell is blank (or absent). A blank axis is an
*unclassified* layer - documented and unit-bearing, but not modelled as a niche axis; to
use it, pass an explicit axis when building the layer.

# Arguments

  - `T`: the `RasterDataSources` dataset type.
  - `code`: which layer of it, in any of the forms [`layerunit`](@ref) accepts.
"""
function layeraxis(T::Type, code)
    a = _layerrow(T, code).axis
    return isempty(a) ? nothing : _resolveaxis(a)
end

"""
    layerrate(unit, period, axis)

Return the unit a layer's values actually carry once its accumulation `period` has been applied - the
declared amount, or that amount per day where the layer's canonical reading is a rate.

# Arguments

  - `unit`: the layer's declared unit, as [`layerunit`](@ref) returns it.
  - `period`: what each value accumulated over - an [`AbstractAccumulationPeriod`](@ref), or
    `nothing` for a layer that accumulates nothing.
  - `axis`: the layer's [`NicheAxis`](@ref), which is what decides whether the canonical reading is
    a rate at all. Not the period - see below.

**Two different questions, kept apart**: [`layerunit`](@ref) answers *what the shipped table
declares*, this answers *what materialising it yields*. Monthly precipitation is catalogued as
`L*m^-2`, an amount, and read as `L*m^-2*d^-1`, a rate.

**Having a period does NOT mean a layer is read as a rate - the axis decides that, not the
period.** `gdd0` accumulates over a year, but `CumulativeHeat`'s canonical unit is `d*K`: the heat
*sum* is the reading that means something, since heat is a condition rather than a consumable
resource. Dividing it would produce a mean daily temperature excess nobody asked for.
So the period says *what interval this accumulated over*, and the axis says *which reading is
canonical*; both are needed and they answer different questions.

Per **day** when it does divide, uniformly, whatever the period - that is what makes twelve slices
divided by twelve different month durations share one unit, and it matches the canonical units the
axes already declare (`mm/day`, `kJ/day`). The per-slice divisor changes the *value*, never the unit.

A [`PerCellAccumulationPeriod`](@ref) is deliberately never divided here: its stock-and-rate
readings are role-dependent and the role is unknown at this point, so it keeps the declared amount.

An axis may be **resource-only** - `CarbonFlux` (`npp`) declares no `Condition` unit at all,
deliberately, because potential productivity is a resource species compete for and not a condition
they are matched against. Its `Resource`-role unit stands in as the canonical reading, so `npp` is
divided like any other rate layer rather than left as an annual total.
"""
layerrate(unit, ::Nothing, axis) = unit

layerrate(unit, ::PerCellAccumulationPeriod, axis) = unit

function layerrate(unit, ::AbstractAccumulationPeriod, axis)
    isnothing(axis) && return unit / Units.day
    # The Condition reading is the canonical one wherever it exists; a resource-only axis has none,
    # so its Resource-role unit answers instead (`canonicalunit(Resource, A)` falls back to the
    # 1-argument form, so this is the Condition unit again for every axis that declares one).
    cu = something(canonicalunit(axis), canonicalunit(Resource, axis),
                   Some(nothing))
    # No canonical unit to consult at all, or the axis wants the accumulated amount itself (a heat
    # sum): leave it alone. Otherwise the canonical reading is a rate, so divide.
    (isnothing(cu) || dimension(unit) == dimension(cu)) && return unit
    return unit / Units.day
end

"""
    layerinfo(T::Type, code)
    layerinfo(code)

Return the catalogue [`LayerRecord`](@ref) for a layer, with its name, definition, unit, niche axis and
sources. The two-argument form looks up `code` in dataset `T`'s table (mirrors `layerunit(T, code)`) and
returns one record. The single-argument form (an `Integer`/`Symbol`/`String` code) searches EVERY shipped
table and returns a `Vector` of all matches - the same code can appear in several datasets (e.g. `bio1` is in
both `BioClim` and `BioClimPlus`).

# Arguments

  - `T`: the dataset to look in. Omit it to search every shipped table instead, which is what
    turns one record into a vector of them.
  - `code`: which layer - an `Integer`/`Symbol`/`String`, as for [`layerunit`](@ref).

The two forms differ in **return type**, not just in scope: name the dataset when you want *the*
record, and omit it when you want to see where a code appears.
"""
function layerinfo(T::Type, code)
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
    layersbyaxis(::Nothing)
    layersbyaxis()

Return a `Vector` of [`LayerRecord`](@ref)s for every layer (across all shipped tables) whose axis is `A` or a
concrete leaf beneath it. Passing an abstract group spans all its axes - `layersbyaxis(TemperatureAxis)`
returns every temperature layer, `layersbyaxis(Temperature)` just that one axis. Printed as a compact list
of codes with short summaries.

`layersbyaxis(nothing)` returns the **unclassified** layers instead - those whose `Axis` cell is blank, which
[`layeraxis`](@ref) also reports as `nothing`. They are documented and unit-bearing but not modelled as a
niche axis, so no axis type can reach them.

`layersbyaxis()` returns **every** layer, classified or not: the whole catalogue.

`layersbyaxis()` exists because `layersbyaxis(NicheAxis)` looks complete and is not - it spans every axis,
but an unclassified layer has none, so it would be silently omitted. No shipped layer is unclassified today,
which is exactly what makes that trap easy to fall into. Iterate the catalogue with the no-argument form.

# Arguments

  - `A`: the axis to search under - a concrete leaf for that axis alone, or an abstract `⋯Axis`
    group to span every axis beneath it. Pass `nothing` for the unclassified layers, or omit it
    entirely for the whole catalogue.
"""
function layersbyaxis(A::Type{<:NicheAxis})
    return filter(r -> !isnothing(r.axis) && r.axis <: A, _catalogue())
end

layersbyaxis(::Nothing) = filter(r -> isnothing(r.axis), _catalogue())

# `copy`, not the cached vector itself: `_catalogue()` hands back the module-level cache, and a caller
# that sorted or filtered it in place would corrupt every later lookup.
layersbyaxis() = copy(_catalogue())

# == Functions ==================================================================================

# **This file names no dataset package, deliberately, and takes `::Type` rather than a
# `RasterDataSources` bound.** The shipped tables are the *package's* knowledge about the data - half
# of what is here (`_leafaxes`, `_resolveaxis`, the accumulation-period family) is about EcoSISTEM's
# own niche axes and not about any source at all - so the catalogue stays in the parent while the
# readers become an extension.
#
# **Why `::Type` and not the `IsRasterData` trait**, which is what `SourceSpec` and `ClimateRaster`
# use: those are **doors**, where an unmarked type would otherwise be stored and surface later; these
# are **lookups**, and `_layerfile` already refuses an unknown type with a message naming the table it
# looked for. Guard the door, not every read - the trait buys nothing here and would cost every public
# name in this file a `@traitfn` rewrite.

# The dataset subtype wrapped by a `RasterDataSource`, e.g. `WorldClim{BioClim}` ->
# `BioClim` and `EarthEnv{LandCover}` -> `LandCover`. Sources without a dataset
# parameter (such as `AWAP`) are returned unchanged.
function _datasettype(::Type{T}) where {T}
    params = Base.unwrap_unionall(T).parameters
    return isempty(params) ? T : first(params)
end

# The layer table shipped in the package `data/RasterDataSources/` directory, named by convention
# after the dataset type (`WorldClim{BioClim}` -> `data/RasterDataSources/BioClim.csv`,
# `EarthEnv{LandCover}` -> `data/RasterDataSources/LandCover.csv`); the file's presence is what
# makes a source supported. (They live under `RasterDataSources/` to keep them distinct from other
# shipped data such as `data/bounding_boxes.csv`.)
function _layerfile(T::Type)
    path = pkgdir(@__MODULE__, "data", "RasterDataSources",
                  "$(nameof(_datasettype(T))).csv")
    isfile(path) ||
        return error("No layer table for raster source `$T` (expected $(basename(path)) in data/)")
    return path
end

# Every row of one shipped catalogue CSV, keyed by layer code, read once and memoised per path.
# Memoised because a catalogue is consulted per layer per build, and re-parsing the file each time
# is what made a catalogue build cost seconds rather than milliseconds.
function _layertable(path::String)
    return get!(_LAYER_CACHE, path) do
        table = CSV.File(path, normalizenames = true)
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

# The leaf `NicheAxis` types - those with no subtypes of their own. `subtypes` only returns direct
# children, so recurse through the intermediate groupings (e.g. `TemperatureAxis`) to reach them.
# The test is "has no children", not `isconcretetype`: *every* axis is an abstract type, so that a
# leaf can be subdivided later without disturbing the hierarchy or any dispatch written against it.
function _leafaxes(T = NicheAxis)
    children = subtypes(T)
    return isempty(children) ? Type[T] :
           mapreduce(_leafaxes, vcat, children, init = Type[])
end

# Resolve an axis name from a table to its `NicheAxis` type by autodiscovery - no registry:
# any loaded leaf `NicheAxis` with that name works (build-time only).
#
# **One walk populates every name, not one walk per name.** Filtering `_leafaxes()` per lookup was
# still ~0.6 s × 35 distinct names even with the results cached, because the cost is the walk and not
# the filter. Building the whole map at once makes it one walk per session.
#
# An ambiguous name is deliberately left **out** of the map, so it still errors on lookup rather
# than silently resolving to whichever came first - the check survives the caching.
function _refreshaxisnames!()
    byname = Dict{String, Vector{Type}}()
    for A in _leafaxes()
        push!(get!(() -> Type[], byname, string(nameof(A))), A)
    end
    empty!(_AXIS_BY_NAME)
    for (n, As) in byname
        length(As) == 1 && (_AXIS_BY_NAME[n] = only(As))
    end
    return _AXIS_BY_NAME
end

# A miss re-walks before giving up, so an axis declared by a downstream package after the map was
# first built is still found; only a genuinely absent or ambiguous name reaches the error.
function _resolveaxis(name::AbstractString)
    key = String(name)
    haskey(_AXIS_BY_NAME, key) && return _AXIS_BY_NAME[key]
    _refreshaxisnames!()
    haskey(_AXIS_BY_NAME, key) && return _AXIS_BY_NAME[key]
    n = count(A -> string(nameof(A)) == key, _leafaxes())
    return error("`$key` does not name exactly one loaded `NicheAxis` (found $n)")
end

# One layer's row from `T`'s catalogue, or an error naming both. `code` is stringified first, so a
# `Symbol`, an `Int` and a `String` all reach the same row - the three spellings the package accepts.
function _layerrow(T::Type, code)
    tbl = _layertable(_layerfile(T))
    key = string(code)
    haskey(tbl, key) ||
        return error("Layer `$code` is not in the table for raster source `$T`")
    return tbl[key]
end

# Whether layer `code` of dataset `S` holds class labels, by asking the axis the catalogue declares
# for it. `iscategorical` is one function throughout - an axis declares it, a dataset and a code look
# it up here, a stack of codes agrees on it or is refused, and a raster asks about its own - and
# every method ends at the axis declaration emitted by `@nicheaxis`, the only place it is stated.
#
# It joins `layerunit`/`layeraxis`/`layerinfo` in taking a bare `::Type` rather than an
# `IsRasterData` bound, for the reason given at the head of this group: these are lookups, not doors,
# and `_layerfile` already refuses an unknown type by naming the table it looked for. Passing an axis
# type as `S` reaches that same refusal, exactly as it does for the rest of the family.
iscategorical(S::Type, code) = iscategorical(layerinfo(S, code).axis)

# The categorical-ness `codes` agree on, or an error naming the ones that disagree. A mixed stack
# has no answer rather than a `false` one: every caller chooses a single behaviour for the whole
# array -- `:mode` against `:bilinear`, `CategoricalRegime` against `ContinuousRegime` -- so either
# answer would be wrong for half the layers and would be applied silently. Interpolating class codes
# is meaningless and reducing measurements to a nearest class discards them, so such a stack is
# refused rather than resolved.
#
# This is the one method of `iscategorical` that can throw, and it is the reason the whole family
# is one function rather than a predicate plus a separately named check: the rule and its message
# then exist once, where two copies had already drifted apart in wording. Reading a multi-layer
# `SourceSpec` calls it for that refusal before downloading anything (`_read`), and asking a raster
# built some other way reaches it through `iscategorical(::ClimateRaster)`.
function iscategorical(S::Type, codes::AbstractVector)
    isempty(codes) && return false
    allequal(iscategorical(S, c) for c in codes) &&
        return iscategorical(S, first(codes))
    cats = filter(c -> iscategorical(S, c), codes)   # built only to name the offenders
    return error("layers $(join(cats, ", ")) hold class codes but " *
                 "$(join(setdiff(codes, cats), ", ")) hold measurements, so they have no one " *
                 "resampling method and no one kind of regime: class codes must not be " *
                 "interpolated, and measurements should not be reduced to a nearest class. Name " *
                 "the layers you want instead - one code for a single layer, or pass a codes " *
                 "vector to `ConstructedRasterSpec`, which reads each on its own terms.")
end

# ---------------------------------------------------------------------------
# Known upstream scaling defects - and a check that fires when one is FIXED
# ---------------------------------------------------------------------------
# **Some providers publish values that are not in the unit they document**, and there is no way to
# tell from the file: the anomaly is in the provider's own metadata, not in ours.
#
# **Why the check is inverted.** These are not defects we can fix once and forget: if the provider
# corrects the data, any compensation *we* apply silently becomes a second error in the opposite
# direction - and a silent factor of 10^4 is exactly what this whole area exists to prevent. So each
# entry records what the values look like **today, while still wrong**, and the check fires when they
# stop looking that way. A failure here is good news that needs acting on, not a bug.
#
# **The discriminator is whether the GeoTIFF carries a scale tag**, not the provider or the dtype -
# checked with `ArchGDAL.getscale`. CHELSA's rasters do carry one and GDAL applies it on read
# (verified: `CHELSA_bio1` returns °C, not raw `UInt16`), so CHELSA needs **no** entry here. EarthEnv's
# do not, which is the whole of the problem below. Following "RasterDataSources applies scaling"
# uniformly is what let this through.
# Parse a `PublishedScaleFactor` cell. Blank => `nothing` (published and documented agree, which is
# the case for all but a handful of layers). Read as a `Rational` so `10000` stays exact and the
# reported factor is never `9999.999...`.
function _parsepublishedscale(s::AbstractString, code)
    isempty(s) && return nothing
    v = tryparse(Int, s)
    isnothing(v) &&
        return error("layer `$code` has `PublishedScaleFactor` = `$s`, which is not an integer. " *
                     "The cell records how many times larger the published values are than the " *
                     "layer's documented unit (e.g. `100`, `10000`); leave it blank if they agree.")
    v > 0 ||
        return error("layer `$code` has a non-positive `PublishedScaleFactor` ($v).")
    return v // 1
end

# The recorded discrepancy for a layer, or `nothing`. Reads the catalogue, so it is visible through
# `layerinfo` rather than hidden here - which is the point of it being a column at all.
# Tolerant of an unknown code: this runs on every read, and a caller may legitimately name a layer
# that has no catalogue row (a hand-built raster). A missing row means nothing is known, not an error.
function _publishedscale(T::Type, code)
    ds = nameof(_datasettype(T))
    key = string(code)
    for r in _catalogue()
        (r.dataset === ds && key in r.aliases) && return r.publishedscale
    end
    return nothing
end

# Check a layer's values still show the defect recorded for it, and say so loudly if they do not.
#
# **A warning, not an error, and deliberately**: a provider fixing its data must not make the
# package unusable until we catch up. The *test* is what makes this loud (`test_LayerCatalogue.jl`);
# this is the belt for anyone reading a dataset the suite never downloads.
# Sampled rather than swept - these rasters run to ~10^9 cells and the check runs on every read.
function _checkupstreamscale(T::Type, code, values)
    factor = _publishedscale(T, code)
    isnothing(factor) && return nothing
    n = denominator(factor) == 1 ? numerator(factor) : factor
    finite = Iterators.filter(isfinite,
                              Iterators.map(Float64, _stridedsample(values)))
    top = reduce(max, finite, init = -Inf)
    isfinite(top) || return nothing
    # **The threshold is the geometric midpoint**, `√factor × ceiling`, which is what makes this
    # sensitive to an order-of-magnitude change and blind to an ordinary data revision. Published
    # values sit at about `factor × ceiling`; corrected ones would sit at about `ceiling`. Splitting
    # the difference on a log scale allows the data to move by up to `√factor` - a factor of **100**
    # for the heterogeneity layers, **10** for `bio4` - before this says anything, while leaving a
    # full rescale nowhere to hide.
    top > sqrt(factor) * _documentedceiling(T, code) && return nothing   # still as recorded
    return @warn "$T layer `$code` no longer looks $(n)× larger than its documented unit - its " *
                 "values now peak at $top. `PublishedScaleFactor` in the shipped table records " *
                 "that this provider publishes values $(n)× too large; if that has been corrected " *
                 "upstream, clear the cell (and remove any compensation keyed on it) - otherwise " *
                 "this package will now be wrong by $(n)× in the other direction."
    # No `maxlog`: it is per **call site**, so one layer reporting would silence every other -
    # and these are distinct, actionable facts about different datasets.
end

# The largest value the layer's documented unit could reach, used only to tell "still scaled" from
# "now corrected". Deliberately crude: only layers with a genuine documented ceiling can be
# checked at all, and `PublishedScaleFactor` is only set on those.
_documentedceiling(::Type, code) = 1.0

# Every ~1000th element, so the cost is independent of raster size.
function _stridedsample(values)
    n = length(values)
    return n <= 1000 ? vec(collect(values)) :
           @view(vec(values)[1:cld(n, 1000):end])
end

# Parse an `AccumulationPeriod` cell. `=` is the discriminator and needs no escaping, because a
# Unitful unit can never contain one: a bare cell is a unit, anything with `=` names how to resolve a
# period that is not constant.
function _parseperiod(cell)
    isempty(cell) && return nothing
    if occursin('=', cell)
        kind, _, rhs = partition_eq(cell)
        kind == "perslice" && return _persliceperiod(rhs, cell)
        kind == "percell" && return PerCellAccumulationPeriod(_periodcode(rhs))
        return error("unrecognised AccumulationPeriod `$cell`: the part before `=` must be " *
                     "`perslice` or `percell`.")
    end
    u = uparse(cell, unit_context = [Unitful, Units])
    dimension(u) == Unitful.𝐓 ||
        return error("AccumulationPeriod `$cell` is not a time unit (got dimension " *
                     "$(dimension(u))) - a period is an interval.")
    return ConstantAccumulationPeriod(u)
end

# `cell` split at its first `=`, so a right-hand side may itself contain one without ambiguity.
function partition_eq(cell)
    i = findfirst('=', cell)
    return cell[1:(i - 1)], '=', cell[(i + 1):end]
end

# `perslice=` right-hand side, checked against the closed vocabulary above.
function _persliceperiod(rhs, cell)
    rhs in _PERSLICE_PERIODS ||
        return error("unrecognised AccumulationPeriod `$cell`: `perslice=` takes one of " *
                     join(_PERSLICE_PERIODS, ", ") * ".")
    return PerSliceAccumulationPeriod()
end

# A `percell=` right-hand side as a layer code: an all-digit reference is a numeric code, anything
# else a symbol, matching how `Code` cells are read everywhere else.
function _periodcode(rhs)
    isempty(rhs) &&
        return error("AccumulationPeriod `percell=` needs a layer code.")
    return all(isdigit, rhs) ? parse(Int, rhs) : Symbol(rhs)
end

# A layer that accumulates over an interval must say which interval, and one that does not must not
# claim one. Both directions are checked here because the `Category` and `AccumulationPeriod`
# columns are written by hand and a mismatch between them is invisible until a rate comes out wrong.
function _checkperiod(period, category, code)
    category in _PERIOD_INHERITED_CATEGORIES && return nothing
    needed = category in _PERIOD_CATEGORIES
    if needed && isnothing(period)
        return error("layer `$code` has Category = :$category, which accumulates over an " *
                     "interval, but its AccumulationPeriod is blank - say what interval, or " *
                     "correct the Category.")
    elseif !needed && !isnothing(period)
        return error("layer `$code` has Category = :$category, which does not accumulate, but " *
                     "declares an AccumulationPeriod - a period there would mean nothing.")
    end
    return nothing
end

# Parse a shipped CSV's `Temporal Resolution` cell (day/month_mean_duration/year, or blank for "no native
# period" - a static/non-time-varying value) into the unit itself, or `nothing`. Validated as a genuine
# time unit at parse time, since a non-time value here would silently break any later use as a native
# reporting period.
function _parsetemporal(cell)
    isempty(cell) && return nothing
    u = uparse(cell, unit_context = [Unitful, Units])
    dimension(u) == Unitful.𝐓 ||
        return error("Temporal Resolution `$cell` is not a time unit (got dimension $(dimension(u)))")
    return u
end

# The `ValueType` cell as a checked `Symbol`, defaulting to `:continuous` when blank. An
# unrecognised spelling is refused rather than passed through, since it would otherwise reach the
# resampler as a silently wrong method.
function _parsevaluetype(cell)
    isempty(cell) && return :continuous
    v = Symbol(lowercase(cell))
    v in _VALUETYPES ||
        return error("unrecognised ValueType \"$cell\"; expected one of $(join(_VALUETYPES, ", "))")
    return v
end

# The `Category` cell as a checked `Symbol`, defaulting to `:instantaneous` when blank. Refuses an
# unrecognised spelling, as `_parsevaluetype` does and for the same reason.
function _parsecategory(cell)
    isempty(cell) && return :instantaneous
    v = Symbol(lowercase(cell))
    v in _CATEGORIES ||
        return error("unrecognised Category \"$cell\"; expected one of $(join(_CATEGORIES, ", "))")
    return v
end

# A unit string reduced to a sorted list of its factors, so two spellings of one unit compare equal.
# Superscripts are rewritten to `^n` first because the catalogue is written with them and Julia's
# parser is not: `kJ d⁻¹` and `d^-1*kJ` are the same dimension and must not read as different ones.
function _normdim(s::AbstractString)
    plain = replace(s,
                    r"[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+" =>
                        m -> "^" * map(c -> _SUPERSCRIPTS[c], m))
    factors = filter(!isempty, strip.(split(plain, r"[*\s]+")))
    return sort(String.(factors))
end

# `UnitDimension` is a hand-maintained copy of `dimension(Units)`, so it can only ever drift from it -
# checked rather than read.
function _checkunitdimension(cell, u, code)
    isempty(cell) && return nothing
    _normdim(cell) == _normdim(string(dimension(u))) ||
        return error("layer `$code` declares UnitDimension `$cell` but its Units `$u` have " *
                     "dimension $(dimension(u)).")
    return nothing
end

# How many days each slice of a read accumulated over, or `nothing` where no division applies.
#
# Returns a **vector** for a per-slice period - one divisor per slice - which is the whole reason
# the conversion cannot be a unit multiply: `_attachunit` puts one scalar unit on a whole array, but
# a twelve-slice monthly read needs twelve different divisors.
#
# `months` comes from the caller's own `readkw`, never from the raster's `Ti` axis: a partial read
# (`month = 2:4`) builds an axis labelled 1:3, so the axis cannot say which real months these are.
# See the deferred partial-read bug in the layer-units plan.
_perioddivisors(::Nothing, months) = nothing

_perioddivisors(::PerCellAccumulationPeriod, months) = nothing

function _perioddivisors(p::ConstantAccumulationPeriod, months)
    return ustrip(Units.day, 1.0 * p.duration)
end

function _perioddivisors(::PerSliceAccumulationPeriod, months)
    return [ustrip(Units.day, 1.0 * Units.month_duration(m)) for m in months]
end

# The divisors a `SourceSpec`-style read must apply, given the layer's record and the months asked
# for - `nothing` when the layer's canonical reading is the accumulated amount itself, so that a heat
# sum is never quietly turned into a temperature.
function _readdivisors(rec::LayerRecord, months)
    layerrate(rec.unit, rec.period, rec.axis) === rec.unit && return nothing
    return _perioddivisors(rec.period, months)
end

# Does a layer's unit agree dimensionally with the canonical unit of the axis it is assigned to?
#
# This is the check whose absence let `canonicalunit(::Type{<:Precipitation}) = mm` (a depth) sit against
# 15 rate-valued precipitation layers, so every one of them threw a bare `DimensionError` from
# `_tocanonu` when built as a regime. Caught here, at catalogue read, it names the layer and both
# units instead. Deliberately worded like the matching check at the trait end
# (`_checksupport`, in `src/Tolerance.jl`), so the two failures read alike.
#
# Skipped when the axis declares no canonical unit (`nothing` - most axes), since there is then
# nothing to disagree with.
#
# Compared against the **rate** the layer resolves to, not the amount the table declares. Since the
# accumulation period was split out of `Units`, a monthly precipitation layer is catalogued `L*m^-2`
# (𝐋) while `Precipitation`'s canonical unit is `mm/day` (𝐋𝐓^-1) - dimensionally different, yet both
# right. It is `layerrate` that reconciles them, so that is what has to be checked; comparing the raw
# cell would now reject every accumulating layer in the tables.
function _checkaxisunit(u, period, axis, code)
    isnothing(axis) && return nothing
    # **Deliberately the Condition unit only, NOT `layerrate`'s Resource fallback** - I tried that
    # and it is wrong. A catalogue row states an **areal flux density** (`npp` is `g*m^-2`), while a
    # Resource-role canonical unit is a **per-cell rate** (`g*day^-1`); the two differ by the cell area
    # on purpose, so comparing them reports a dimension mismatch on a correct pair. `layerrate` can
    # consult the Resource unit because it only asks *whether the dimensions differ* in order to
    # decide about dividing by time; this check asks whether they **agree**, which is a stricter
    # question the Resource unit cannot answer.
    cu = canonicalunit(axis)
    isnothing(cu) && return nothing
    r = layerrate(u, period, axis)
    dimension(r) == dimension(cu) ||
        return error("layer `$code` resolves to unit $r but its axis $(nameof(axis))'s canonical " *
                     "unit $cu has a different dimension ($(dimension(r)) vs $(dimension(cu))); " *
                     "one of the two is wrong." *
                     (r === u ? "" :
                      " (The table declares $u, which its accumulation period turns into $r.)"))
    return nothing
end

# Every shipped CSV must have exactly the columns expected of it - no missing ones, and no extras.
# Both halves matter: a renamed header reads as blank everywhere it is used rather than failing, and
# an unexpected column is usually a rename that half-happened.
function _checkschema(dataset::Symbol, cols)
    allowed = (_REQUIRED_COLUMNS..., get(_OPTIONAL_COLUMNS, dataset, ())...)
    missingcols = setdiff(_REQUIRED_COLUMNS, cols)
    isempty(missingcols) ||
        return error("$dataset.csv is missing required column(s) $(join(missingcols, ", ")); " *
                     "a renamed header would otherwise read as blank everywhere it is used.")
    extra = setdiff(cols, allowed)
    isempty(extra) ||
        return error("$dataset.csv has unexpected column(s) $(join(extra, ", ")); " *
                     "add them to `_REQUIRED_COLUMNS`/`_OPTIONAL_COLUMNS` and read them, or " *
                     "remove them - an unread column is data nothing honours.")
    return nothing
end

# Do two accumulation periods say the same thing? One method per pair of kinds, so a `Constant` is
# never silently "equal" to a `PerSlice` that happens to average the same length.
_periodsagree(::Nothing, ::Nothing) = true

function _periodsagree(a::ConstantAccumulationPeriod,
                       b::ConstantAccumulationPeriod)
    return a.duration == b.duration
end

_periodsagree(::PerSliceAccumulationPeriod, ::PerSliceAccumulationPeriod) = true

function _periodsagree(a::PerCellAccumulationPeriod,
                       b::PerCellAccumulationPeriod)
    return a.code == b.code
end

_periodsagree(_, _) = false

# How a period reads when a check has to name one, including the blank case `_periodphrase` never sees.
_periodwording(::Nothing) = "no period at all"

_periodwording(p) = _periodphrase(p)

# The whole-catalogue invariant for `*Range` axes - the one check `_checkperiod` cannot make alone,
# because a range **inherits** its period rather than having one of its own. A range is the difference
# of two values of the quantity it ranges over, and subtraction preserves both unit and accumulation
# period, so a range row must carry `absoluteunit` of its siblings' unit and exactly their period.
#
# `absoluteunit`, not the unit itself, because a difference of *positions* is an *interval*: that is
# why `bio2` (diurnal temperature range) is `K` where `Temperature` is `°C`. Unitful does this for
# free - `°C -> K`, `°F -> Ra`, identity on everything non-affine.
#
# Resolved against **same-dataset** siblings, not the axis alone, and that distinction bites:
# `SolarRadiation` carries both `MJ*m^-2` (CHELSA `rsds`) and `kJ*m^-2` (WorldClim `srad`), so an
# axis-wide lookup could hand `rsds_range` the wrong unit. Where siblings disagree among themselves
# there is nothing to check against and the row is skipped rather than guessed at.
# Do all the shipped rows on one axis resolve to the same dimension?
#
# **The check that needs nothing declared**, which is precisely where `_checkaxisunit` gives up: that
# one compares a row against its axis's canonical unit and returns early when there is none, so an
# axis nobody has declared is unchecked from *both* ends. Comparing the rows to **each other** needs
# no declaration, so it covers exactly that gap - and it keeps covering it for an axis a downstream
# package declares loosely.
#
# Iterates the **rows**, not the axes, so an axis with no shipped rows - derived from the grid, or
# belonging to another package - simply never appears, rather than being reported as an empty
# disagreement.
#
# Grouped across **all** datasets, unlike `_checkrangerows` which groups within one: reconciling
# WorldClim's `kJ*m^-2` with CHELSA's `MJ*m^-2` on one axis is the point, and a per-dataset grouping
# would never see the pair. **Dimensions, not units** - differing *scale* on one axis is legal, and
# reconciling it is what `canonicalunit` is for. Compared after `layerrate`, for the same reason
# `_checkaxisunit` is: a monthly total and a daily rate are both right and differ dimensionally.
function _checkaxishomogeneity(catalogue)
    first_on = Dict{Type, LayerRecord}()
    for r in catalogue
        isnothing(r.axis) && continue
        f = get!(first_on, r.axis, r)
        f === r && continue
        rate, frate = layerrate(r.unit, r.period, r.axis),
                      layerrate(f.unit, f.period, f.axis)
        dimension(rate) == dimension(frate) ||
            error("layers `$(first(f.aliases))` and `$(first(r.aliases))` are both on axis " *
                  "$(nameof(r.axis)) but resolve to different dimensions - $frate " *
                  "($(dimension(frate))) and $rate ($(dimension(rate))). One row's Units, " *
                  "AccumulationPeriod or Axis cell is wrong: rows sharing an axis may differ in " *
                  "scale, never in dimension.")
    end
    return nothing
end

# A `range` layer must agree with its siblings on the same axis: it carries their unit made absolute,
# and their accumulation period. Checked across the whole catalogue rather than per row, because the
# constraint is between rows - a range with no siblings has nothing to disagree with.
function _checkrangerows(catalogue)
    axisname(r) = isnothing(r.axis) ? "" : String(nameof(r.axis))
    family = Dict{Tuple{Symbol, String}, Vector{LayerRecord}}()
    for r in catalogue
        key = (r.dataset, axisname(r))
        push!(get!(() -> LayerRecord[], family, key), r)
    end
    for r in catalogue
        ax = axisname(r)
        endswith(ax, "Range") || continue
        base = ax[1:(end - length("Range"))]
        sibs = get(family, (r.dataset, base), LayerRecord[])
        isempty(sibs) && continue
        units = unique(s.unit for s in sibs)
        length(units) == 1 || continue
        all(s -> _periodsagree(s.period, first(sibs).period), sibs) || continue
        code = first(r.aliases)
        want = Unitful.absoluteunit(only(units))
        r.unit == want ||
            error("layer `$code` is a range on $ax, so it is the difference of two $base " *
                  "values and must carry their unit made absolute - `$want` - but its Units " *
                  "cell says `$(r.unit)`. (Absolute because a difference of positions is an " *
                  "interval: that is why a temperature range is K where a temperature is °C.)")
        _periodsagree(r.period, first(sibs).period) ||
            error("layer `$code` is a range on $ax, so it spans the same interval as the $base " *
                  "values it is the difference of - which accumulate $(_periodwording(first(sibs).period)) " *
                  "- but it declares $(_periodwording(r.period)).")
    end
    return nothing
end

# Every record from every shipped CSV, built once and memoised. This is the single place the data
# directory is walked, so a dataset is in the catalogue exactly when its file is there.
function _catalogue()
    isempty(_CATALOGUE) || return _CATALOGUE
    datadir = pkgdir(@__MODULE__, "data", "RasterDataSources")
    for f in sort(filter(endswith(".csv"), readdir(datadir)))
        dataset = Symbol(first(splitext(f)))
        table = CSV.File(joinpath(datadir, f), normalizenames = true)
        cols = propertynames(table)
        _checkschema(dataset, cols)
        # `string` (not `String`) so numeric columns like `NumSlices` (an Int) convert too
        cell(row, col) = (col in cols && !ismissing(getproperty(row, col))) ?
                         String(strip(string(getproperty(row, col)))) : ""
        splitsemis(s) = filter(!isempty, String.(strip.(split(s, ";"))))
        for row in table
            ismissing(row.Code) && continue
            aliases = splitsemis(string(row.Code))
            isempty(aliases) && continue
            u = cell(row, :Units)
            a = cell(row, :Axis)
            ns = cell(row, :NumSlices)
            code = first(aliases)
            unit = isempty(u) ? NoUnits :
                   uparse(u, unit_context = [Unitful, Units])
            axis = isempty(a) ? nothing : _resolveaxis(a)
            category = _parsecategory(cell(row, :Category))
            period = _parseperiod(cell(row, :AccumulationPeriod))
            # The three cross-checks: a period is declared exactly where one applies, the duplicated
            # dimension column must agree with the unit, and the unit must be usable on its axis.
            _checkperiod(period, category, code)
            _checkunitdimension(cell(row, :UnitDimension), unit, code)
            _checkaxisunit(unit, period, axis, code)
            optional(col) = col in cols ? cell(row, col) : nothing
            push!(_CATALOGUE,
                  LayerRecord(dataset, aliases, cell(row, :Name),
                              cell(row, :Definition), unit, axis,
                              splitsemis(cell(row, :Source)),
                              isempty(ns) ? 1 : parse(Int, ns),
                              _parsetemporal(cell(row, :Temporal_Resolution)),
                              period,
                              _parsevaluetype(cell(row, :ValueType)), category,
                              cell(row, :OfficialUnit), cell(row, :Notes),
                              optional(:Group), optional(:Order),
                              _parsepublishedscale(cell(row,
                                                        :PublishedScaleFactor),
                                                   code)))
        end
    end
    # After the loop, not inside it: a range row is checked against its siblings, which may sit
    # anywhere in the same table (or be read after it).
    _checkrangerows(_CATALOGUE)
    _checkaxishomogeneity(_CATALOGUE)
    return _CATALOGUE
end

# The abstract-supertype chain of a leaf axis, up to (but excluding) `NicheAxis` - e.g.
# `Temperature` -> `[TemperatureAxis]`; a direct leaf like `SolarRadiation` -> `[]`.
function _axischain(A::DataType)
    chain = DataType[]
    S = supertype(A)
    while S !== NicheAxis && S !== Any
        push!(chain, S)
        S = supertype(S)
    end
    return chain
end

# How a unit reads in a record's detail block. `NoUnits` is spelled out rather than printed, since an
# empty string beside a layer's name reads as missing data rather than as a deliberate absence.
_unitstr(u) = u == NoUnits ? "dimensionless" : string(u)

# How an accumulation period reads in a record's detail block. One method per kind, so a period that
# is a unit and one that is a layer reference each say what they actually mean rather than printing a
# struct.
_periodphrase(p::ConstantAccumulationPeriod) = "over 1 " * string(p.duration)

function _periodphrase(::PerSliceAccumulationPeriod)
    return "over each slice's own calendar month (31 d for January, 28.25 d for a mean February ...)"
end

function _periodphrase(p::PerCellAccumulationPeriod)
    return "over the `" * string(p.code) * "` layer, which varies by cell"
end

# Shipped layer names (the CSV `Name` column, via `LayerRecord.name`) that use axis `A` directly.
function _axisnames(A::Type{<:NicheAxis})
    return [r.name for r in _catalogue() if r.axis === A]
end

# Build one `AxisNode` for `A`, recursing through `subtypes` down to the concrete leaves.
function _axisnode(A::Type{<:NicheAxis})
    children = sort(_axisnode.(subtypes(A)), by = c -> string(nameof(c.axis)))
    return AxisNode(A, _axisnames(A), children)
end

# One node of the axis tree and everything under it, indented by depth. Recursive, and the only
# printer for the tree - `layerinfo`'s axis view is this function all the way down.
function _showtree(io::IO, node::AxisNode, depth::Int)
    println(io, "  "^depth, nameof(node.axis),
            isempty(node.names) ? "" : " - " * join(node.names, ", "))
    for child in node.children
        _showtree(io, child, depth + 1)
    end
end

"""
    layeraxes(A::Type{<:NicheAxis} = NicheAxis)

Return the niche-axis hierarchy at and below `A` (the whole tree by default) as a nested [`AxisNode`](@ref):
each node carries the shipped layer names that use its axis directly (only a concrete leaf ever has any -
an abstract grouping node's own `names` is always empty) plus its child axis nodes, recursively down to the
concrete leaves. Use it to discover which axes exist and what shipped layers use them, then drill into a
group with [`layersbyaxis`](@ref).

# Arguments

  - `A`: the axis to root the tree at. Defaults to `NicheAxis`, the whole hierarchy; pass a group
    such as `TemperatureAxis` to see only that branch.
"""
layeraxes(A::Type{<:NicheAxis} = NicheAxis) = _axisnode(A)
