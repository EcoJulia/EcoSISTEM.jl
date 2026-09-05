# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The data-backed layer recipes: name a catalogued source, a vector file, or a function combining
# other specs. None holds any data - each is resolved against the target grid at build time.

using DimensionalData

"""
    AbstractLazySpec <: AbstractSpec

Abstract supertype of the lazy, data-backed / derived specs. Resolved against the target grid at
build time and usable in *either* role - a regime/supply layer or an active mask: [`SourceSpec`](@ref)
(read a data source), `ShapeSpec` (a vector file), `ConstructedRasterSpec` (combine child specs by a
function).
"""
abstract type AbstractLazySpec <: AbstractSpec end

"""
    LayerSpec

Type-union of everything accepted as a regime or supply layer: a synthetic layer spec
([`AbstractSyntheticLayerSpec`](@ref)) or any lazy data-backed spec ([`AbstractLazySpec`](@ref)).
"""
const LayerSpec = Union{AbstractSyntheticLayerSpec, AbstractLazySpec}

"""
    LayerInput

Type-union of everything a `regime` or `supply` keyword accepts: one [`AbstractSpec`](@ref) - which
includes a [`Varying`](@ref) wrapping one - or a `Tuple`/`NamedTuple` of them for a multi-layer
environment, a named tuple keeping the caller's names.

**Written into the builders' signatures**, not merely documented, so that `methods(GridHabitat)` and
the rendered docs both *show* what a builder accepts, and a wrong kind of argument is rejected where
it was passed rather than several calls later.

Deliberately **not** the same as [`LayerSpec`](@ref), which is the union of spec *types*: this also
admits the tuple forms and the [`Varying`](@ref) wrapper, because those are things a caller writes at
a keyword rather than kinds of layer recipe. A `Tuple` here always means **several layers**, one spec
per member; name a single data layer with a [`SourceSpec`](@ref), which is where its axis, unit and
read options live in any case.

Two costs, both accepted. Julia does not dispatch on keyword types, so a wrong argument is refused by
a `TypeError` naming the keyword and printing this union expanded rather than by a message suggesting
a remedy - and no fallback method can improve on that, since a second method with the same positional
signature replaces the first instead of adding to it. And a signature cannot see **inside** a
container, so an element of a tuple that is not a spec is caught later, by the resolvers.
"""
const LayerInput = Union{AbstractSpec, Tuple, NamedTuple}

"""
    MaskSpec

Type-union of everything accepted as an `active` mask: a synthetic mask spec
([`AbstractSyntheticMaskSpec`](@ref)) or any lazy data-backed spec ([`AbstractLazySpec`](@ref)).
"""
const MaskSpec = Union{AbstractSyntheticMaskSpec, AbstractLazySpec}

"""
    SourceSpec{A <: NicheAxis, U, K <: NamedTuple}

Name a layer of a catalogued data source, without reading it. It holds **no** grid array: the read,
the cut and the resample happen only when it is materialised onto a decided grid.

Passing **no** `code` describes the *whole* dataset - every layer read into one multi-band raster,
which is the form [`ConstructedRasterSpec`](@ref) uses for a bare dataset, such as all the land-cover class
bands for `compress_landcover`.

# Arguments

  - `source`: the dataset type to read from, `WorldClim{BioClim}` and the like.
  - `code`: which layer of it, as one [`CODE_TYPE`](@ref) or a vector of them. Omit it for the whole
    dataset.
  - `unit`: the physical unit to attach on read. Defaults to the layer's own, from the shipped table.
  - `axis`: the niche axis, which is what matches the layer to species niches. Defaults from the
    shipped table, and to [`NicheAxis`](@ref) where the table names none.
  - any other keyword: kept as a pass-through argument for the eventual `read`, so
    `SourceSpec(WorldClim{Climate}, :wind, month = 1:12)` reads the twelve monthly layers and
    `month = 1` just the one. A `SourceSpec` nested inside a [`ConstructedRasterSpec`](@ref) can therefore
    carry its own read options.

    Two of those keywords decide how much is read, and are worth knowing about together. `cut`
    windows the read, so only the cells inside a bounding box come off disk. `scale` coarsens by an
    integer factor, and on its own it coarsens the *whole* source file however small a result is
    wanted, because the aggregated form is memoised per file rather than per window: for a global
    dataset that first read can need many gigabytes. Give a `scale` a `cut` as well where the whole
    world is not needed, and the memo is skipped along with the cost.

# Fields

  - `source`, `code`, `unit`, `readkw`: as above. `code` is never `nothing` - a whole-dataset spec
    resolves the dataset's own code list at construction, so every layer's identity is known before
    anything is read and each can keep its own unit. `unit` falls back to `NoUnits` as a neutral
    placeholder where a multi-layer spec's layers disagree.

# Type parameters

  - `A`: the niche axis. A type parameter rather than a field because it is dispatched on, which is
    what lets a layer's meaning be checked at compile time rather than looked up.
  - `U`, `K`: the types of `unit` and `readkw`.
"""
struct SourceSpec{A <: NicheAxis, U, K <: NamedTuple} <:
       EcoSISTEM.AbstractLazySpec
    source::Type
    # One layer, or several. Never `nothing`: `SourceSpec(dataset)` resolves the dataset's own code
    # list here rather than carrying a late-bound "all of them" sentinel, so there is one code shape
    # instead of three, and every layer's identity is known before anything is read - which is what
    # lets each keep its *own* unit.
    code::Union{CODE_TYPE, Vector{CODE_TYPE}}
    unit::U
    readkw::K   # extra keywords forwarded to `read` at materialise time (e.g. `month = 1:12`)
    # The sole constructor - defined with the type's own name and taking `axis` as a runtime *type*
    # argument (see `GradientSpec`'s inner constructor comment), so there is exactly one way to build
    # one. Omitting `code` gives the whole dataset. `unit`/`axis` are resolved in the *body* rather
    # than as signature defaults because their defaults are shipped-table lookups keyed on `code`.
    # Any keyword other than `axis` is captured as a pass-through read keyword, so a `SourceSpec`
    # nested inside a `ConstructedRasterSpec` can specify e.g. its own `month`.
    #
    # A multi-layer spec does **not** error here when its layers disagree on unit or axis, even
    # though it cannot then honestly claim one. Four of the seven shipped datasets are heterogeneous -
    # including `WorldClim{BioClim}` (6 units) and `CHELSA{BioClimPlus}` (13 units, 29 axes) - so
    # refusing them at construction would rule out the flagship sources. Its real use is inside a
    # `ConstructedRasterSpec`, where `_parselayers` expands it to one correctly-united spec per layer and
    # the disagreement never arises. The error belongs where a *single array* is genuinely required -
    # materialising it directly as a regime or supply - and lives in `_read` accordingly.
    #
    # **Trait-gated on `IsRasterData`, not bounded by `RasterDataSources.RasterDataSource`** -
    # the same treatment `ClimateRaster`'s sole constructor already has, and for the same reason: a
    # `<:` bound names one package, whereas the trait asks the question the code actually cares
    # about and a third party's raster type can answer it in one `@traitimpl` line. It is also what
    # lets this struct stay here while `RasterDataSources` is a weak dependency.
    @traitfn function SourceSpec(::Type{S},
                                 code::Union{CODE_TYPE,
                                             AbstractVector{<:CODE_TYPE},
                                             Nothing} = nothing,
                                 unit = nothing;
                                 axis::Union{Type{<:NicheAxis}, Nothing} = nothing,
                                 readkw...) where {S; IsRasterData{S}}
        source = S
        c = isnothing(code) ? EcoSISTEM._alllayercodes(source) :
            code isa AbstractVector ? collect(CODE_TYPE, code) : code
        u = !isnothing(unit) ? unit : _sharedunit(source, c)
        a = !isnothing(axis) ? axis : _sharedaxis(source, c)
        kw = NamedTuple(readkw)
        return new{a, typeof(u), typeof(kw)}(source, c, u, kw)
    end
end

# Must be defined here rather than beside `ERA`/`CERA` in `Climate.jl`: `SourceSpec` is declared in
# this file, and a method written above its own type's file is silently discarded -- the package
# still loads, and the call fails later with a bare `MethodError`.
# Worth the extra method, exactly as on `ClimateRaster`: without it an unmarked source fails with a
# bare `MethodError` naming `SimpleTraits.Not{IsRasterData{...}}`, which leaks the trait machinery and
# names no remedy. It also covers the case this file now cares most about - a user who has not
# loaded `RasterDataSources` and so cannot name a dataset at all.
@traitfn function SourceSpec(::Type{S}, args...;
                             kw...) where {S; !IsRasterData{S}}
    return error("`$S` cannot name a data source. Load the package that defines the dataset " *
                 "(`using RasterDataSources` for the shipped ones), or mark your own raster " *
                 "type with `@traitimpl EcoSISTEM.IsRasterData{$S}`.")
end

"""
    AbstractShapeSpec

A mask that is a piece of **ground** rather than data - [`ShapeSpec`](@ref) for one given as a vector
file, [`NaturalEarthSpec`](@ref) for one given by name, [`ConstructedShapeSpec`](@ref) for several
combined or one transformed.

What these share, and what the abstract type is for, is that they resolve to **geometry** before any
grid exists. That is what lets them be combined exactly and at no resolution - a study area of your
own, unioned with a country you named - where combining rasterised masks would have to fix a
resolution before the study grid had been decided. It is the vector mirror of
[`ConstructedRasterSpec`](@ref), which composes rasters and so does need a grid.
"""
abstract type AbstractShapeSpec <: EcoSISTEM.AbstractLazySpec end

"""
    ShapeSpec(path::AbstractString; layer = 0)

Name an active-area mask taken from the polygons of a vector file, without reading it. It holds
**no** geometry: the read, any download, the per-feature reprojection into the target grid's own CRS
and the cell-membership test all happen when it is materialised onto a decided grid, as for
[`SourceSpec`](@ref).

# Arguments

  - `path`: a shapefile, GeoJSON, GeoPackage, or any format GDAL reads. A path ending in `.zip` is
    read directly, with no need to unzip first, and a URL is downloaded into
    `EcoSISTEM.assetdir(owner = ShapeSpec)` as an [`EcoSISTEM.CachedAsset`](@ref) the first time it
    is needed.
  - `layer`: which layer of the file, 0-indexed. Every polygon feature in it is used.
"""
struct ShapeSpec <: AbstractShapeSpec
    path::Union{String, EcoSISTEM.CachedAsset}
    layer::Int
    # A leading URL scheme (`scheme://...`) marks `path` as a download, deferred to a `CachedAsset`;
    # anything else is taken to be an already-local path, used as-is.
    function ShapeSpec(path::AbstractString; layer::Integer = 0)
        p = occursin(r"^[a-zA-Z][a-zA-Z0-9+.-]*://", path) ?
            EcoSISTEM.CachedAsset(ShapeSpec, path) : String(path)
        return new(p, Int(layer))
    end
end

"""
    NaturalEarthSpec(name::AbstractString; level = nothing, coverage = AllTerritories(),
                     outline = true)

Name an active-area mask as a **named region** - a country, a continent, an island - without
reading anything. The polygons are fetched and cut to the grid when the spec is materialised, as for
[`ShapeSpec`](@ref).

The name is resolved against the shipped region table at construction, so a name that does not exist
is an error where it was written rather than minutes later mid-build. It is resolved by exactly the
rule [`boundingbox`](@ref) uses, which is what makes the box that function reports the box this
spec's shape actually has.

# Arguments

  - `name`: the region's name, matched case-insensitively but otherwise as Natural Earth spells it.
  - `level`: which kind of region the name means - `"ADMIN"` for a country, `"Physical Island"` for a
    landmass; `EcoSISTEM.naturalearth_levels()` lists them. Only needed where a name means genuinely
    different ground at different levels, and the error says so when it does.
  - `coverage`: how much of what the name covers to take - [`AllTerritories`](@ref), the default and
    what the source itself means by the name, or [`LargestLandmass`](@ref) for the principal landmass
    alone.
  - `outline`: `true`, the default, activates only the cells whose centres fall inside the region.
    `false` activates every cell in the region's bounding box instead, which is the cheaper thing to
    want when the region is only being used to say *where* to work rather than to mask a coastline.
"""
struct NaturalEarthSpec{C <: EcoSISTEM.AbstractCoverage} <: AbstractShapeSpec
    level::String
    name::String
    coverage::C
    outline::Bool

    function NaturalEarthSpec(name::AbstractString; level = nothing,
                              coverage::EcoSISTEM.AbstractCoverage = AllTerritories(),
                              outline::Bool = true)
        lvl = isnothing(level) ? EcoSISTEM._resolvelevel(name, coverage) :
              EcoSISTEM._checklevel(level).name
        row = EcoSISTEM._regionrow(lvl, name)
        isnothing(row) &&
            error("No region named \"$name\" at level \"$lvl\". " *
                  "`EcoSISTEM.naturalearth_levels()` lists the levels.")
        # The source's own spelling is stored, not the caller's: the lookup is case-insensitive, and
        # what is kept should be what the data says so that `show` and any later report agree with it.
        return new{typeof(coverage)}(lvl, row.Name, coverage, outline)
    end
end

"""
    NaturalEarthSpec(match::EcoSISTEM.RegionMatch; coverage = AllTerritories(), outline = true)

Turn one match from [`investigate_regions`](@ref) into a spec, without naming it again.

A match already carries the level and the name, which is the whole of a spec's identity, so nothing
is re-derived and the shape agrees with the box the report displayed.

A *report* cannot be converted, because it may hold several regions. Pick one first - `only(report)`
asserts there was exactly one, `first(report)` takes the best by the report's own ordering, and
`report[i]` takes a chosen one.
"""
function NaturalEarthSpec(match::EcoSISTEM.RegionMatch;
                          coverage::EcoSISTEM.AbstractCoverage = AllTerritories(),
                          outline::Bool = true)
    return NaturalEarthSpec(match.name, level = match.level.name,
                            coverage = coverage, outline = outline)
end

# A report is ambiguous by construction, so converting one would have to pick silently. `first` is
# meaningful under `Encloses`, whose order is smallest-enclosing-first, and not under the other two -
# it would be right a third of the time. An error naming the three ways to choose beats that, and a
# `MethodError` would name none of them.
function NaturalEarthSpec(report::EcoSISTEM.RegionReport; kw...)
    return error("A `RegionReport` holds $(length(report)) region" *
                 (length(report) == 1 ? "" : "s") *
                 ", so it does not name one spec. Choose: `only(report)` asserts there was exactly " *
                 "one, `first(report)` takes the best by the report's own ordering, `report[i]` " *
                 "takes the one you want.")
end

function Base.show(io::IO, s::NaturalEarthSpec)
    print(io, "NaturalEarthSpec(\"", s.name, "\", level = \"", s.level, "\"")
    EcoSISTEM._isdefaultcoverage(s.coverage) ||
        print(io, ", coverage = ", s.coverage)
    s.outline || print(io, ", outline = false")
    return print(io, ")")
end

"""
    ConstructedShapeSpec(operation, members...; coverage = AllTerritories(), outline = true)

Combine several named regions into one mask - the union of the United Kingdom, Ireland and the Isle
of Man, or a country with an island group cut out of it.

Regions combine as **geometry**, so the result is exact and carries no resolution of its own: the
grid is still decided afterwards, and nothing is rasterised twice.

```julia
# The British Isles, including Shetland - which Natural Earth's own polygon of that name omits
ConstructedShapeSpec(ShapeUnion(),
                   NaturalEarthSpec("United Kingdom", coverage = AllTerritories()),
                   NaturalEarthSpec("Ireland", level = "ADMIN"),
                   NaturalEarthSpec("Isle of Man", level = "ADMIN"),
                   coverage = LandmassesAbove(1km^2))     # ...and without Rockall
```

# Arguments

  - `operation`: how they combine - [`ShapeUnion`](@ref), [`ShapeIntersection`](@ref) or
    [`ShapeDifference`](@ref), the last taking every later member away from the first.
  - `members`: two or more region specs, either [`NaturalEarthSpec`](@ref)s or nested
    `ConstructedShapeSpec`s.
  - `coverage`: which components of the *result* to keep, applied after the operation -
    [`AllTerritories`](@ref) by default, since a combination usually means all of what it built.
  - `outline`: as [`NaturalEarthSpec`](@ref) - `false` activates the result's bounding box instead of
    its outline.
"""
struct ConstructedShapeSpec{O, M <: Tuple,
                            C <: EcoSISTEM.AbstractCoverage} <:
       AbstractShapeSpec
    operation::O
    members::M
    coverage::C
    outline::Bool

    function ConstructedShapeSpec(operation::Union{EcoSISTEM.AbstractShapeOperation,
                                                   Function},
                                  members::AbstractShapeSpec...;
                                  coverage::EcoSISTEM.AbstractCoverage = AllTerritories(),
                                  outline::Bool = true)
        least = EcoSISTEM._minmembers(operation)
        length(members) >= least ||
            throw(ArgumentError("`$operation` needs at least $least shape" *
                                (least == 1 ? "" : "s") *
                                "; it was given $(length(members))."))
        return new{typeof(operation), typeof(members), typeof(coverage)}(operation,
                                                                         members,
                                                                         coverage,
                                                                         outline)
    end
end

function Base.show(io::IO, s::ConstructedShapeSpec)
    print(io, "ConstructedShapeSpec(", s.operation, ", ")
    join(io, s.members, ", ")
    EcoSISTEM._isdefaultcoverage(s.coverage) ||
        print(io, ", coverage = ", s.coverage)
    s.outline || print(io, ", outline = false")
    return print(io, ")")
end

"""
    ConstructedRasterSpec(combine, layers...)

The universal lazy escape hatch: read each of `layers` onto the working grid, then apply `combine`
to the resulting rasters. Because `combine` is the **first** argument it can be written as a
do-block:

```julia
ConstructedRasterSpec(EarthEnv{LandCover}, :open_water) do water
    water .< 50   # a mask of cells less than half open water
end
```

`layers` are given as alternating `dataset, code(s)` - a `RasterDataSources` type followed by an
`Int`/`Symbol` code, a vector/tuple of codes (several layers), or **no** code (a bare dataset =>
*all* its layers, e.g. every land-cover class band, passed to `combine` as one multi-band raster);
a pre-built spec is also accepted directly - including a **synthetic** one, so a combine may mix
generated layers with read ones. A bare dataset becomes a whole-dataset spec. With **no** layers,
`combine` is a nullary thunk that produces the layer itself (reading a source directly, or wrapping
a literal in-memory array).

**`combine` takes rasters and must return a raster** - one contract, whichever way the spec is later
used, because that decision is made elsewhere and cannot be known where the combine is written. A
**mask is simply a raster whose element type is `Bool`**, so the element type still distinguishes the
two; only the container is fixed.

Nothing has to be wrapped or unwrapped to satisfy that: a raster broadcasts and yields a raster, as
the example above shows, and `sum(bands)` over a multi-band combine does the same. So a combine names
**no array type at all**, which is the point - the array type is an implementation detail, and a
combine is user code.

Usable as a regime/supply layer or an active mask. Covers what the other specs don't:
bespoke data, derived layers (anomalies, blends) and thresholded/combined masks. See `hasdata`
and `landcoverclass` for ready-made combine building blocks. No data is read or downloaded at
construction - only when materialised onto a grid (`GridHabitat`), mirroring
[`SourceSpec`](@ref); each layer's unit/axis is resolved from the shipped table eagerly, so an
invalid code errors here rather than at materialise time.

`combinestage` says *when* `combine` runs - after its layers are put on the study grid
([`CombineOnTargetGrid`](@ref), the default) or before ([`CombineOnSourceGrid`](@ref)). See
[`AbstractCombineStage`](@ref) for which a given combine needs.

# Fields

  - `combine`: the rule itself - a function taking one raster per layer (none, for a thunk) and
    returning a raster.
  - `layers`: the child specs to materialise and hand it, already normalised; empty for a thunk.
  - `combinestage`: when the combine runs, as above.

`axis` is a **type parameter** rather than a field, as on [`SourceSpec`](@ref) - and on this type it
is the *only* statement of what the result means, since a derived raster has no layer code to resolve
one from. That includes whether the result holds class codes: there is no `valuetype` keyword,
because a `TypologyAxis` says the values are class labels and so must be resampled by nearest class,
while any other axis says they may be interpolated. A separate declaration could only agree with the
axis or contradict it.
"""
struct ConstructedRasterSpec{A <: NicheAxis, F} <: EcoSISTEM.AbstractLazySpec
    axis::Type{A}  # the niche axis of the produced layer (matched to species tolerances); mask => ignored
    combine::F
    # **`AbstractSpec`, not `Vector{SourceSpec}`** - two things at once. It is what lets this type
    # live outside `ClimatePref` (a `SourceSpec` is defined *after* this file, so naming it here would
    # be a cycle), and it is what lets a combine take a **synthetic** layer, which the argument parser
    # used to refuse outright. Abstractly typed, which costs nothing here: layers are walked once
    # per materialisation, never in a hot loop.
    layers::Vector{AbstractSpec}
    # A runtime field rather than a type parameter, for the same reason `ClimateRaster`'s `code`
    # is: it is consulted once per materialisation and never in a hot loop, so the single dynamic
    # dispatch it costs buys nothing back for multiplying the concrete spec types.
    combinestage::AbstractCombineStage
    # `axis` is a required keyword (as on every other spec); a derived regime layer (e.g. a
    # temperature anomaly) declares what it measures, while a mask - which is never paired with a
    # tolerance - says `NicheAxis` to state that it is claiming nothing.
    #
    # **Whether the result holds class codes comes from `axis`, and is not declared separately.**
    # It matters only alongside `combinestage = CombineOnSourceGrid()`, where the combine's own
    # result is what gets sampled - on the default path the layers are sampled first, so nothing ever
    # interpolates a class code. The two remain independent: `gsp / gsl` must collapse early (a
    # ratio does not commute with regridding) and produces perfectly ordinary continuous values.
    function ConstructedRasterSpec(combine, layerargs...;
                                   axis::Type{A},
                                   combinestage::AbstractCombineStage = CombineOnTargetGrid()) where {A <:
                                                                                                      NicheAxis}
        return new{A, typeof(combine)}(axis, combine,
                                       _parselayers(layerargs...),
                                       combinestage)
    end
end

# _sharedunit(source, code) / _sharedaxis(source, code)
#
# The unit and niche axis a `SourceSpec` takes when its caller does not state them - read from the
# shipped catalogue, and for a vector of codes only where they agree.
#
# In this file, beside the constructor that asks: the catalogue moved into the parent with it, so
# there is no submodule boundary left to cross.
# The unit/axis a spec can honestly claim: a single layer's own, the one its layers agree on, or the
# neutral value when they do not. `NoUnits`/`NicheAxis` for a disagreeing multi-layer spec is a
# placeholder, not a claim - such a spec is only materialisable through `_parselayers`, which expands
# it into per-layer specs that each carry the right one, and `_read` refuses it otherwise.
# **Option C**: `layerunit` answers what the shipped table declares; a `SourceSpec`'s `unit`
# answers what *materialising it* yields - which, for a layer with an accumulation period whose
# canonical reading is a rate, is the declared amount per day. The two are different questions and
# now have different answers, instead of one field quietly meaning both.
function _sharedunit(source, code::CODE_TYPE)
    rec = layerinfo(source, code)
    return layerrate(rec.unit, rec.period, rec.axis)
end
function _sharedunit(source, codes::AbstractVector)
    us = unique(layerunit(source, c) for c in codes)
    return length(us) == 1 ? only(us) : NoUnits
end
# The axis a spec's layers share, or `NicheAxis` where the catalogue names none. A multi-layer spec
# must resolve to one axis, since a collection built from it is named by axis - so the vector method
# is where a mixed-axis request is caught.
function _sharedaxis(source, code::CODE_TYPE)
    return something(layeraxis(source, code), NicheAxis)
end
function _sharedaxis(source, codes::AbstractVector)
    as = unique(layeraxis(source, c) for c in codes)
    return length(as) == 1 ? something(only(as), NicheAxis) :
           NicheAxis
end

# --- Display ------------------------------------------------------------------
# As in `Spec.jl`: the one-liner is the expression that builds it, with optional arguments shown
# only where they are not at their default.
#
# `ConstructedRasterSpec` is the one that cannot follow the rule, and says so: its `combine` is an
# arbitrary function with no readable spelling, so the line reports what it is built *from* instead.
function Base.show(io::IO, spec::SourceSpec{A}) where {A}
    kw = isempty(spec.readkw) ? "" :
         ", " * join(("$(k) = $(v)" for (k, v) in pairs(spec.readkw)), ", ")
    return print(io,
                 "SourceSpec($(spec.source), $(repr(spec.code))$(kw), axis = $(nameof(A)))")
end

function Base.show(io::IO, spec::ShapeSpec)
    layer = iszero(spec.layer) ? "" : ", layer = $(spec.layer)"
    return print(io, "ShapeSpec($(repr(spec.path))$(layer))")
end

function Base.show(io::IO, spec::ConstructedRasterSpec{A}) where {A}
    n = length(spec.layers)
    return print(io,
                 "ConstructedRasterSpec($(n) layer$(n == 1 ? "" : "s"), axis = $(nameof(A)))")
end

function Base.show(io::IO, ::MIME"text/plain", spec::ConstructedRasterSpec)
    println(io, sprint(show, spec))
    for l in spec.layers
        println(io, "  ", sprint(show, l))
    end
    return nothing
end

# --- Desugaring and labelling a spec -----------------------------------------
# A per-cell supply is rewritten into a combine here, before any grid is decided.

# How a layer is named in that message: the dataset it comes from and the code asked for.
_speclabel(spec::SourceSpec) = "`$(spec.source)` layer `$(spec.code)`"

# A multi-variable `regime`/`supply` is a *tuple* of specs, each of which shapes the grid in its own
# right. A tuple therefore always means "several layers", at every level - which is why the bare
# `(source, code)` pair form had to go: being itself a tuple, it could be told from a multi-layer
# regime by nothing but nesting depth. `_sourcepairnotaspec` refuses one and names `SourceSpec`.
# Each element is also stripped of any `Varying` wrapper: a declared change has no bearing on the
# grid, so every consumer here - `_probecrs`, `_shapesgrid`, `_asraster` - must see the naked spec.
# Unwrapping in this one place covers all of them, and both roles, because they all iterate this
# output. Without it the wrapper would be *silently accepted* rather than rejected: `_shapesgrid`
# would misclassify a wrapped synthetic spec as data-shaping, and `_probecrs`'s `::Any` fallback
# would decline for the whole area, quietly disabling read windowing.
# **A per-cell accumulation period, desugared.** A layer whose `AccumulationPeriod` is
# `percell=<code>` (today only `gsp`, growing-season precipitation over `gsl`, growing-season length)
# holds an *amount*, and the interval it accumulated over varies by cell. As a **regime** that is
# exactly what is wanted - "how much water over the season" - and nothing needs doing. As a **supply**
# it must become a rate, which means reading that other layer and dividing.
#
# **Driven by the catalogue, never by the code name.** The rewrite fires on
# `PerCellAccumulationPeriod` wherever it is declared, so a second such layer needs no change here -
# a hard-coded `:gsp` would be exactly the second source of truth this project spent Step 5b removing.
#
# **`CombineOnSourceGrid()` is required, not a preference.** Division is cell-wise but *nonlinear*,
# so it does not commute with regridding: native cells (100 mm, 50 d) and (100 mm, 100 d) with a
# target straddling both give 1.5 mm/d if divided early and 1.33 mm/d if divided late. Step 1 built
# this stage for exactly this consumer.
#
# `gsl == 0` needs no policy: the division yields `NaN`, and `_coverage` then marks
# the cell inactive. A cell with no growing season has no growing-season water - the right answer, free.
function _desugarsupply(spec::SourceSpec)
    rec = _percellrecord(spec)
    isnothing(rec) && return spec
    divisor = SourceSpec(spec.source, rec.period.code)
    # `Precipitation`, not `GrowingSeasonPrecipitation` (user, 2026-08-05): a Resource-role axis
    # says *which resource*, not which layer it came from, and water is water. The provenance
    # therefore lives in the spec, not in the built layer - and this now matters, because the built
    # layer is `Supply{Precipitation}` *because the axis said so*, not because its values happened
    # to come out as `L/day`.
    # A combine is handed `ClimateRaster`s and must return one - the wrapper carries the source and
    # code that `_sampledeclared` needs to sample the *result* on the early-collapse path, so dividing
    # the bare arrays and returning that would strip the grid provenance and be refused there.
    source = spec.source
    return ConstructedRasterSpec(spec, divisor, axis = _percellaxis(spec, rec),
                                 combinestage = CombineOnSourceGrid()
                                 ) do amount,
                                      period
        return ClimateRaster(source,
                             _perperiod.(amount.array,
                                         period.array))
    end
end

# A tuple/named tuple of supplies desugars member-wise, keeping its names.
_desugarsupply(spec::Union{Tuple, NamedTuple}) = map(_desugarsupply, spec)

# Anything else - a synthetic spec, an already-built supply, a `Varying` wrapper - passes through. A
# `Varying` is left wrapped deliberately: `_expandspecs` unwraps it later for grid decisions, and
# rewriting inside it here would drop the declared change.
_desugarsupply(spec) = spec

# The catalogue record for `spec` when it declares a per-cell period, else `nothing`. Multi-code specs
# are declined rather than guessed at: a per-cell period is a property of one layer, and a stacked
# read has no single divisor.
function _percellrecord(spec::SourceSpec)
    spec.code isa CODE_TYPE || return nothing
    rec = try
        layerinfo(spec.source, spec.code)
    catch
        return nothing
    end
    return rec.period isa PerCellAccumulationPeriod ? rec : nothing
end

# The niche axis a per-cell layer's *rate* reading belongs to. Only water exists today; anything else
# is refused by name rather than silently given the wrong axis, since guessing here would build a
# supply of the wrong resource.
function _percellaxis(spec::SourceSpec, rec::LayerRecord)
    rec.axis <: WaterAxis && return Precipitation
    return error("`$(spec.code)` accumulates over the `$(rec.period.code)` layer, so as a supply it " *
                 "is a rate - but its axis `$(nameof(rec.axis))` has no rate reading defined here. " *
                 "Add one to `_percellaxis` naming the resource it supplies.")
end
