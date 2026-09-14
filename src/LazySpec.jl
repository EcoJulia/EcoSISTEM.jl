# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The data-backed layer recipes: name a catalogued source, a vector file, or a function combining
# other specs. None holds any data - each is resolved against the target grid at build time. Also
# what can be said about a recipe's files without reading them: their provenance records and
# checksums.

using DimensionalData

import Extents

"""
    AbstractLazySpec <: AbstractSpec

Abstract supertype of the lazy, data-backed / derived specs. Resolved against the target grid at
build time and usable in *either* role - a regime/supply layer or an active mask: [`RasterSpec`](@ref)
(read raster data - a catalogued layer, or a file that belongs to no dataset), `ShapeSpec` (a
vector file), `ConstructedRasterSpec` (combine child specs by a function).
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
    RasterSpec{A <: NicheAxis, U} <: AbstractLazySpec

Name raster data to be read, without reading it - a layer of a catalogued data source, written
`SourceSpec(source, code; ...)`, or a raster file that belongs to no dataset, written
`RasterFileSpec(path; axis, ...)`. It holds **no** grid array: the read, the cut and the resample
happen only when it is materialised onto a decided grid, so a global file costs the window the
study area needs rather than the globe, and refining an area re-reads nothing.

Four things say what is read. **What** the data is: the `source` type and the `code` of the
layer, which is what the shipped catalogue is keyed on for a unit and an axis, and which the caller
may override; a file that belongs to no dataset has no code, so it must be told both. **Where** the
files are: an explicit list, or `nothing` for a source that resolves its own, as a catalogued
dataset does. **How** it is read: `cut` windows the read to a box, `scale` coarsens by a whole
number of cells a side, and `fn` says how a block is reduced. And **how the files become one array
in time**: `times` labels the slices where the files cannot, `atend` says what the series does past
its last slice, and `calendar` what its coordinates mean.

# Fields

  - `source`: the data source, a type marked [`EcoSISTEM.IsRasterData`](@ref) - a dataset such
    as `WorldClim{BioClim}`, or [`SyntheticData`](@ref) for a file that belongs to no dataset.
  - `code`: which layer of the source, as one [`CODE_TYPE`](@ref) or a vector of them; `nothing`
    for a bare file, and never for a catalogued source, since a whole-dataset spec resolves the
    dataset's own code list at construction so that every layer's identity is known before
    anything is read and each can keep its own unit.
  - `files`: the files to read, in time order, each a local path, the
    [`EcoSISTEM.CachedAsset`](@ref) a URL becomes, or a [`EcoSISTEM.CDSRequest`](@ref) fetched on
    first use; `nothing` where the source resolves them itself. Each file is opened by the
    backend its source's catalogue row names, so a netCDF file without an extension reads.
  - `unit`: the physical unit attached on read. `NoUnits` where a multi-layer spec's layers
    disagree, as a neutral placeholder, and for a file given none.
  - `cut`: an `Extents.Extent` of `°` intervals to window the read to, or `nothing` for all of it.
    Applied before the pixels are fetched, so cutting a global layer to one country costs the
    country.
  - `scale`: an integer factor to coarsen by on read, each block of `scale × scale` cells becoming
    one, or `nothing` where the study area may choose one to suit its grid. A stated scale, `1`
    included, stands. A coarsened read of a whole file is memoised on disk, so the cost is paid
    once per file rather than once per study area; a `cut` alongside skips the memo and the cost.
  - `fn`: how a block is reduced to one cell, or `nothing` for the axis to decide - the most
    frequent class (ties to the smallest code) for a `TypologyAxis`, whose codes must not be
    averaged, and the mean for any other. On a grid much coarser than the layer the reduction
    runs in two stages, on the layer's own lattice and then onto the grid, so it must compose: a
    mean, a maximum or a minimum do exactly, the most frequent class approximately.
  - `times`: the coordinates of the third axis, one per slice, replacing what the files say -
    for files carrying no time coordinate, or to label them in elapsed time rather than calendar
    dates; `nothing` to take the files' own, or monthly ordinals where they have none.
  - `atend`: what the series the stack becomes does past its last slice, an
    [`AbstractSeriesEnd`](@ref); [`RepeatAtEnd`](@ref) by default, which a climatology of even
    months cycles by. A series of real dates cannot cycle evenly and says [`HoldAtEnd`](@ref) or
    [`ErrorAtEnd`](@ref).
  - `calendar`: what the coordinates mean, an [`AbstractSeriesCalendar`](@ref), or `nothing` to
    infer it: real dates give a [`DatedSeries`](@ref), anything else an [`UndatedSeries`](@ref); a
    monthly climatology says [`MonthOfYearSeries`](@ref) itself.
  - `readkw`: keywords the source needs to resolve its files, kept for the eventual read -
    `month = 1:12` for a monthly climatology, say. A spec nested inside a
    [`ConstructedRasterSpec`](@ref) can therefore carry its own read options.

# Type parameters

  - `A`: the niche axis, which is what matches the layer to species niches. A type parameter
    rather than a field because it is dispatched on, which is what lets a layer's meaning be
    checked at compile time rather than looked up.
  - `U`: the type of `unit`.
"""
struct RasterSpec{A <: NicheAxis, U} <: EcoSISTEM.AbstractLazySpec
    source::Type
    code::Union{Nothing, CODE_TYPE, Vector{CODE_TYPE}}
    files::Union{Nothing,
                 Vector{Union{String, EcoSISTEM.CachedAsset,
                              EcoSISTEM.CDSRequest}}}
    unit::U
    cut::Union{Nothing, Extents.Extent}
    scale::Union{Nothing, Int}
    # `nothing` for *decide from the axis*. Untyped beyond that because it is consulted once per
    # materialisation and never in a hot loop.
    fn::Union{Nothing, Function}
    times::Union{Nothing, AbstractVector}
    atend::AbstractSeriesEnd
    calendar::Union{Nothing, AbstractSeriesCalendar}
    readkw::NamedTuple

    # The one route to `new`. The two spellings, `SourceSpec` and `RasterFileSpec`, decide the
    # fields and come here; the checks that hold whichever spelling was used live here with it.
    function RasterSpec{A, U}(source::Type, code, files, unit::U, cut,
                              scale::Union{Nothing, Integer}, fn, times,
                              atend::AbstractSeriesEnd,
                              calendar::Union{Nothing, AbstractSeriesCalendar},
                              readkw::NamedTuple) where {A <: NicheAxis, U}
        (isnothing(scale) || scale >= 1) ||
            error("`scale` coarsens by a whole number of cells per side, so it must be at " *
                  "least 1; got $scale.")
        isnothing(files) || !isempty(files) ||
            error("`files` names no file; give at least one path or URL.")
        isnothing(times) || !isempty(times) ||
            error("`times` names no time; give one per slice, or leave it out.")
        return new{A, U}(source, code, files, unit, cut,
                         isnothing(scale) ? nothing : Int(scale), fn, times,
                         atend, calendar, readkw)
    end
end

"""
    SourceSpec(source, code = nothing, unit = nothing; axis, files = nothing, file = nothing,
               directory = nothing, cut = nothing, scale = nothing, fn = nothing,
               times = nothing, atend = RepeatAtEnd(), calendar = nothing, readkw...)

Name a layer of a catalogued data source, without reading it, as a [`RasterSpec`](@ref) - one whose
`files` the source resolves for itself, or one reading files you name, which the catalogue still
describes: `SourceSpec(ERA, "t2m", file = path)` reads variable `t2m` out of an ERA5 netCDF file
with the unit, axis and accumulation period the `ERA` table gives it, and
`SourceSpec(CRUTS, "tavg", directory = dir)` reads a directory of monthly GeoTIFFs as one series.

Passing **no** `code` describes the *whole* dataset - every layer read into one multi-band raster,
which is the form [`ConstructedRasterSpec`](@ref) uses for a bare dataset, such as all the
land-cover class bands for `compress_landcover`.

# Arguments

  - `source`: the dataset type to read from, `WorldClim{BioClim}` and the like.
  - `code`: which layer of it, as one [`CODE_TYPE`](@ref) or a vector of them. Omit it for the whole
    dataset.
  - `unit`: the physical unit to attach on read. Defaults to the layer's own, from the shipped table.
  - `axis`: the niche axis, which is what matches the layer to species niches. Defaults from the
    shipped table, and to [`NicheAxis`](@ref) where the table names none.
  - `files`, `file`, `directory`: where the data is, when the source does not fetch it itself -
    a vector of files in time order, one file, or a directory whose raster files are taken in name
    order; each a local path, a URL, or a [`EcoSISTEM.CDSRequest`](@ref). At most one of the three.
    Required for a source whose catalogue row says its files are not fetched for it (`ERA`, `CERA`,
    `CRUTS`), and refused for none. A source whose row says its files come by plain download
    (`TwentyCR`) fetches the layer's own file into the asset cache on first use unless one of the
    three names a copy already held.
  - `cut`, `scale`, `fn`: the read options, as the [`RasterSpec`](@ref) fields of those names. Give
    a `scale` a `cut` as well where the whole world is not needed: on its own a `scale` coarsens the
    *whole* source file however small a result is wanted, because the aggregated form is memoised
    per file, and for a global dataset that first read can need many gigabytes.
  - `times`, `atend`, `calendar`: how the files become one series, as the [`RasterSpec`](@ref)
    fields of those names. A netCDF file's own dates are kept unless `times` replaces them.
  - any other keyword: kept as a pass-through argument for the eventual read, so
    `SourceSpec(WorldClim{Climate}, :wind, month = 1:12)` reads the twelve monthly layers and
    `month = 1` just the one.
"""
const SourceSpec = RasterSpec

# Defined with the type's alias name and taking `axis` as a runtime *type* argument (see
# `GradientSpec`'s inner constructor comment). Omitting `code` gives the whole dataset. `unit`/`axis`
# are resolved in the *body* rather than as signature defaults because their defaults are
# shipped-table lookups keyed on `code`. The three read options are named so that they land in
# their own fields; any keyword after them is captured as a pass-through read keyword.
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
                             files = nothing, file = nothing,
                             directory = nothing,
                             cut = nothing,
                             scale::Union{Nothing, Integer} = nothing,
                             fn::Union{Nothing, Function} = nothing,
                             times = nothing,
                             atend::AbstractSeriesEnd = RepeatAtEnd(),
                             calendar::Union{Nothing, AbstractSeriesCalendar} = nothing,
                             readkw...) where {S; IsRasterData{S}}
    c = isnothing(code) ? EcoSISTEM._alllayercodes(S) :
        code isa AbstractVector ? collect(CODE_TYPE, code) : code
    u = !isnothing(unit) ? unit : _sharedunit(S, c)
    a = !isnothing(axis) ? axis : _sharedaxis(S, c)
    fs = _specfiles(S, c, files, file, directory)
    _checkfetchable(S, fs)
    return RasterSpec{a, typeof(u)}(S, c, fs, u, cut, scale, fn, times, atend,
                                    calendar, NamedTuple(readkw))
end

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

# One-liner is the spelling that rebuilds the spec, which the files decide: `SourceSpec(...)` for a
# spec whose source resolves them, `RasterFileSpec(...)` for one that names them. Dispatched on the
# `files` field rather than branched, so each spelling is its own method.
Base.show(io::IO, spec::RasterSpec) = _showspec(io, spec, spec.files)

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
    ShapeSpec(path::AbstractString; layer = 0, coverage = AllTerritories(), outline = true)

Name an active-area mask taken from the polygons of a vector file, without reading it. It holds
**no** geometry: the read, any download, the dissolve of its features into connected pieces of
ground, the reprojection into the target grid's own CRS and the cell-membership test all happen
when it is materialised onto a decided grid, as for [`SourceSpec`](@ref); `read(spec)` gives the
pieces themselves.

# Arguments

  - `path`: a shapefile, GeoJSON, GeoPackage, or any format GDAL reads. A path ending in `.zip` is
    read directly, with no need to unzip first, and a URL is downloaded into
    `EcoSISTEM.assetdir(owner = ShapeSpec)` as an [`EcoSISTEM.CachedAsset`](@ref) the first time it
    is needed. A URL must name a **self-contained** file - a `.zip`, a GeoJSON or a GeoPackage.
    Only the one named file is fetched, so a bare remote `.shp` cannot be read: its `.shx`, `.dbf`
    and `.prj` companions never arrive and GDAL refuses the result. Point at the zip the shapefile
    is published in instead, or download the set by hand and give the local `.shp` path.
  - `layer`: which layer of the file, 0-indexed. Every polygon feature in it is used.
  - `coverage`: how much of the ground the file covers to take, once its features are dissolved
    into connected pieces - [`AllTerritories`](@ref), the default and everything the file holds,
    [`LargestLandmass`](@ref) for the principal piece, or [`LandmassesAbove`](@ref) for every
    piece clearing a threshold.
  - `outline`: `true`, the default, activates only the cells whose centres fall inside the
    file's polygons. `false` activates every cell in their bounding box instead, as for
    [`NaturalEarthSpec`](@ref).
"""
struct ShapeSpec{C <: EcoSISTEM.AbstractCoverage} <: AbstractShapeSpec
    path::Union{String, EcoSISTEM.CachedAsset}
    layer::Int
    coverage::C
    outline::Bool
    # A leading URL scheme (`scheme://...`) marks `path` as a download, deferred to a `CachedAsset`;
    # anything else is taken to be an already-local path, used as-is.
    function ShapeSpec(path::AbstractString; layer::Integer = 0,
                       coverage::EcoSISTEM.AbstractCoverage = AllTerritories(),
                       outline::Bool = true)
        p = occursin(r"^[a-zA-Z][a-zA-Z0-9+.-]*://", path) ?
            EcoSISTEM.CachedAsset(ShapeSpec, path) : String(path)
        return new{typeof(coverage)}(p, Int(layer), coverage, outline)
    end
end

function Base.show(io::IO, spec::ShapeSpec)
    print(io, "ShapeSpec(", repr(spec.path))
    iszero(spec.layer) || print(io, ", layer = ", spec.layer)
    EcoSISTEM._isdefaultcoverage(spec.coverage) ||
        print(io, ", coverage = ", spec.coverage)
    spec.outline || print(io, ", outline = false")
    return print(io, ")")
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

Combine several shapes into one mask - the union of the United Kingdom, Ireland and the Isle of
Man, a country with an island group cut out of it, or a study area of your own buffered by a
distance - or transform one.

Shapes combine as **geometry**, so the result is exact and carries no resolution of its own: the
grid is still decided afterwards, and nothing is rasterised twice; `read(spec)` gives the pieces
of ground the result is.

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
    [`ShapeDifference`](@ref), the last taking every later member away from the first, each
    wanting two or more members; or how one is transformed - [`ShapeBuffer`](@ref),
    [`ShapeSimplify`](@ref) or [`ShapeConvexHull`](@ref), each wanting exactly one; or a function
    handed one geometry per member and returning one.
  - `members`: any shape specs - a [`ShapeSpec`](@ref) of your own, a [`NaturalEarthSpec`](@ref),
    or a nested `ConstructedShapeSpec` - as many as the operation wants.
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
    ConstructedRasterSpec(combine, layers...; axis, combinestage = CombineOnTargetGrid(),
                          atend = RepeatAtEnd(), calendar = nothing)

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
  - `atend`, `calendar`: for a combine returning a raster with a time axis, what the series it
    becomes does past its last slice and what its coordinates mean, as on [`SourceSpec`](@ref).

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
    # How a three-dimensional result becomes a series, as on `RasterSpec`: a derived stack reaches
    # `_setseries!` exactly as a read one does.
    atend::AbstractSeriesEnd
    calendar::Union{Nothing, AbstractSeriesCalendar}
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
                                   combinestage::AbstractCombineStage = CombineOnTargetGrid(),
                                   atend::AbstractSeriesEnd = RepeatAtEnd(),
                                   calendar::Union{Nothing,
                                                   AbstractSeriesCalendar} = nothing) where {A <:
                                                                                             NicheAxis}
        return new{A, typeof(combine)}(axis, combine,
                                       _parselayers(layerargs...),
                                       combinestage, atend, calendar)
    end
end

# `ConstructedRasterSpec` is the one that cannot follow the rule, and says so: its `combine` is an
# arbitrary function with no readable spelling, so the line reports what it is built *from* instead.
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

"""
    RasterFileSpec(path::AbstractString; axis, unit = NoUnits, source = SyntheticData,
                   cut = nothing, scale = nothing, fn = nothing, times = nothing,
                   atend = RepeatAtEnd(), calendar = nothing)

Name a raster **file** as a layer, without reading it: a GeoTIFF or anything else GDAL reads, given
by path or URL, that is not a layer of a catalogued dataset. Builds a [`RasterSpec`](@ref) whose
`files` is that one file and whose `code` is `nothing`, so it is read, cut and resampled exactly as
a [`SourceSpec`](@ref) is.

This is the lazy counterpart of [`in_memory_raster`](@ref): that wraps a raster already read into
memory, in full and cached nowhere; this holds only the path. Prefer a [`SourceSpec`](@ref) where the
file *is* a layer of a catalogued dataset, since the catalogue then supplies the unit and axis. A file
has no catalogue, so this spec must be told both.

Its values must be **intensive** - a rate, a density, a state, a fraction - because a layer reaches
the grid by averaging the cells that cover each grid cell, and a count per cell would be averaged
too and read as a density. Divide a count by the cell's area before naming the file.

It is a layer spec. To use a file as a `within` mask, say what in it marks a cell active with a
[`ConstructedRasterSpec`](@ref) over it - `ConstructedRasterSpec(r -> .!isnan.(r), spec, axis = NicheAxis)`.

# Arguments

  - `path`: the file, or a URL naming a **self-contained** file (a `.zip` is read directly). A URL
    is downloaded into `EcoSISTEM.assetdir(owner = RasterSpec)` as an
    [`EcoSISTEM.CachedAsset`](@ref) the first time it is needed.
  - `axis`: the [`NicheAxis`](@ref) the values are on. Required: pass `NicheAxis` itself for data
    whose meaning is not being claimed.
  - `unit`: the physical unit the file's values are in, attached on read. Defaults to `NoUnits`, so
    a file of temperatures in kelvin needs `unit = K` to become a temperature regime.
  - `source`: the data source recorded on the raster - [`SyntheticData`](@ref) by default, for a
    file that belongs to no catalogued dataset.
  - `cut`, `scale`, `fn`: the read options, as the [`RasterSpec`](@ref) fields of those names.
  - `times`, `atend`, `calendar`: for a file with a third axis, how its slices become a series, as
    the [`RasterSpec`](@ref) fields of those names.
"""
function RasterFileSpec(path::AbstractString; axis::Type{A}, unit = NoUnits,
                        source::Type = SyntheticData, cut = nothing,
                        scale::Union{Nothing, Integer} = nothing,
                        fn::Union{Nothing, Function} = nothing,
                        times = nothing,
                        atend::AbstractSeriesEnd = RepeatAtEnd(),
                        calendar::Union{Nothing, AbstractSeriesCalendar} = nothing) where {A <:
                                                                                           NicheAxis}
    return RasterSpec{A, typeof(unit)}(source, nothing, [_fileentry(path)],
                                       unit, cut, scale, fn, times, atend,
                                       calendar, NamedTuple())
end

# == Functions ======================================================================================

function provenance(spec::RasterSpec)
    entries = isnothing(spec.files) ? _presentfiles(spec) : spec.files
    return Union{Nothing, InputRecord}[provenance(_localpath(e))
                                       for e in entries]
end

function provenance(spec::ShapeSpec)
    return Union{Nothing, InputRecord}[provenance(_localpath(spec.path))]
end

function provenance(spec::NaturalEarthSpec)
    path = _localpath(_nesource(_checklevel(spec.level)))
    return Union{Nothing, InputRecord}[_regionrecord(path)]
end

function provenance(spec::ConstructedShapeSpec)
    return reduce(vcat, (provenance(m) for m in spec.members),
                  init = Union{Nothing, InputRecord}[])
end

function provenance(spec::ConstructedRasterSpec)
    return reduce(vcat, (_layerprovenance(l) for l in spec.layers),
                  init = Union{Nothing, InputRecord}[])
end

"""
    verifyassets(spec::RasterSpec)

Check every file a spec reads that is present against the checksum its provenance record holds,
erroring on the first that differs - a truncated or replaced copy - and return how many were
checked. A file with no record, or whose record carries no checksum, has nothing to check
against and is not counted; nothing is fetched.

# Arguments

  - `spec`: the spec whose files to check.
"""
function verifyassets(spec::RasterSpec)
    entries = isnothing(spec.files) ? _presentfiles(spec) : spec.files
    checked = 0
    for path in _localpath.(entries)
        isfile(path) || continue
        checked += _verifyfile(path)
    end
    return checked
end

# The files a source has on disk for a spec that resolves its own, fetching nothing.
function _presentfiles(spec::RasterSpec)
    readkw = (; _getrasterkw(spec.source)..., spec.readkw...)
    return _localfiles(spec.source, spec.code; readkw...)
end

# The record of every file a read of `spec` came from, once the read has made them present: the
# record beside each where one of ours is, else the file's name, all with the catalogue row's
# facts where the spec names a layer of a catalogued dataset. Nothing is hashed.
function _readinputs(spec::RasterSpec; role::Symbol = :habitat)
    entries = isnothing(spec.files) ? _presentfiles(spec) : spec.files
    return InputRecord[_inputrecord(spec, path, role)
                       for path in _localpath.(entries) if isfile(path)]
end

# One file's record for a read of `spec`: its sidecar's, or its name under `role`, with the row's
# facts.
function _inputrecord(spec::RasterSpec, path::AbstractString, role::Symbol)
    rec = _specrecord(spec)
    record = something(_ourrecordat(path),
                       InputRecord(role = role,
                                   dataset = isnothing(rec) ? "file" :
                                             rec.dataset,
                                   code = spec.code, path = basename(path)))
    return isnothing(rec) ? record : _withcatalogue(record, rec, path)
end

# A combination member's entries: a data layer's own, and none for a synthetic one, which reads no
# file. Deliberately untyped, as the fallback for every member a combination accepts.
function _layerprovenance(layer::Union{RasterSpec, ConstructedRasterSpec})
    return provenance(layer)
end

_layerprovenance(::Any) = Union{Nothing, InputRecord}[]

# A Natural Earth zip's record, or `nothing` where it has none, with the source's row filled in.
function _regionrecord(path::AbstractString)
    record = provenance(path)
    isnothing(record) && return nothing
    return _withcatalogue(record, datasetinfo(NaturalEarthLevel), path)
end

# `record` with what the catalogue row knows filled in where it has nothing: the dataset's name,
# its DOI, licence and citation, and the version the file states, or else the row's.
function _withcatalogue(record::InputRecord, rec::DatasetRecord,
                        path::AbstractString)
    orrow(have, row) = isempty(have) ? row : have
    version = isempty(record.version) ?
              something(_fileversion(rec, path), rec.version) : record.version
    return InputRecord(role = record.role,
                       dataset = record.dataset == "file" ? rec.dataset :
                                 record.dataset,
                       code = record.code, path = record.path, url = record.url,
                       request = record.request, job = record.job,
                       fetched = record.fetched, bytes = record.bytes,
                       sha256 = record.sha256,
                       doi = orrow(record.doi, rec.doi),
                       licence = orrow(record.licence, rec.licence),
                       version = version,
                       citation = orrow(record.citation, rec.citation))
end

# The `files` field from whichever of the three location keywords was given: `files` as they are,
# `file` as a one-entry list, `directory` as its raster files in name order. A path with a URL
# scheme becomes a `CachedAsset`, as on `ShapeSpec`. Given none, a source whose row says `https`
# names the layer's own download, owned by the source so its cache directory is its own.
function _specfiles(S, code, files, file, directory)
    given = count(!isnothing, (files, file, directory))
    given <= 1 ||
        error("give one of `files`, `file` or `directory`, not $given of them.")
    isnothing(files) || return _fileentry.(collect(files))
    isnothing(file) || return [_fileentry(file)]
    isnothing(directory) || return _directoryfiles(S, code, directory)
    return _httpsfile(S, code)
end

# The layer's download as a one-entry `files` list, for a source whose row says `https` and a
# layer whose row names a `File`; `nothing` for any other source, or a layer naming none.
function _httpsfile(S, code)
    rec = EcoSISTEM._datasetrecord(S)
    (isnothing(rec) || rec.fetch !== :https || !(code isa CODE_TYPE)) &&
        return nothing
    url = layerinfo(S, code).file
    isnothing(url) && return nothing
    return [EcoSISTEM.CachedAsset(S, url)]
end

_fileentry(x::Union{EcoSISTEM.CachedAsset, EcoSISTEM.CDSRequest}) = x

function _fileentry(path::AbstractString)
    return occursin(r"^[a-zA-Z][a-zA-Z0-9+.-]*://", path) ?
           EcoSISTEM.CachedAsset(RasterSpec, path) : String(path)
end

# The raster files in `dir`, in name order, which is the time order every archive read this way
# names them in. Filtered by the extension the source's catalogue row implies - a netCDF file may
# also carry none, as a Climate Data Store download does - so a directory holding the files and
# their provenance sidecars reads cleanly; a source with no row takes every regular file. A netCDF
# directory is filtered further to the files whose header holds the code's variable, since one
# directory commonly holds every variable of an archive.
function _directoryfiles(S, code, dir::AbstractString)
    isdir(dir) || error("`directory = $(repr(dir))` is not a directory.")
    rec = EcoSISTEM._datasetrecord(S)
    exts = isnothing(rec) ? nothing :
           rec.format === :netCDF ? (".nc", "") : (".tif", ".tiff")
    names = filter(sort(readdir(dir))) do f
        startswith(f, ".") && return false
        isfile(joinpath(dir, f)) || return false
        return isnothing(exts) || lowercase(last(splitext(f))) in exts
    end
    if !isnothing(rec) && rec.format === :netCDF && code isa CODE_TYPE
        names = filter(f -> _holdsvariable(joinpath(dir, f), Symbol(code)),
                       names)
    end
    isempty(names) &&
        error("`directory = $(repr(dir))` holds no " *
              (isnothing(exts) ? "files" : join(exts, "/") * " files") *
              (code isa CODE_TYPE ? " holding `$code`" : "") * " for `$S`.")
    return _fileentry.(joinpath.(dir, names))
end

# Whether a netCDF file holds a variable of that name: a header open, no pixels, and a file that
# cannot be opened as one does not.
function _holdsvariable(path::AbstractString, name::Symbol)
    return try
        EcoSISTEM._lazyopen(path, source = Rasters.NCDsource(), name = name)
        true
    catch
        false
    end
end

# A source that does not fetch its own files must be told where they are, and is refused here,
# where the keyword was omitted, rather than at read time.
_checkfetchable(S, ::AbstractVector) = nothing

function _checkfetchable(S, ::Nothing)
    rec = EcoSISTEM._datasetrecord(S)
    (isnothing(rec) || rec.fetch === :getraster) && return nothing
    how = rec.fetch === :cds ?
          " - or a `CDSRequest` entry in `files`, fetched on first use" :
          rec.fetch === :https ?
          " - its layer table names no `File` to download for this layer" :
          ""
    return error("`$S` does not fetch its own files (its catalogue row says `$(rec.fetch)`), so " *
                 "say where they are: `file = path`, `files = [...]` or `directory = dir`$how.")
end

# A spec's path as text, for a label or a cache key: the string itself, or a download's URL.
_pathtext(path::AbstractString) = String(path)
_pathtext(asset::EcoSISTEM.CachedAsset) = asset.url
_pathtext(request::EcoSISTEM.CDSRequest) = request.path

# What a spec's stack does in time, for `_setseries!`: its own `atend` and `calendar`, or the
# defaults for a spec that has none.
_seriespolicy(spec::RasterSpec) = (atend = spec.atend, calendar = spec.calendar)

function _seriespolicy(spec::ConstructedRasterSpec)
    return (atend = spec.atend, calendar = spec.calendar)
end

_seriespolicy(::Any) = (atend = RepeatAtEnd(), calendar = nothing)

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
    return EcoSISTEM._foldedunit(layerrate(rec.unit, rec.period, rec.axis),
                                 EcoSISTEM._thickness(rec), rec.axis)
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
# The read and series options a spec states, as `name = value` for either spelling.
function _readoptions(spec::RasterSpec)
    opts = String[]
    isnothing(spec.cut) || push!(opts, "cut = $(spec.cut)")
    isnothing(spec.scale) || push!(opts, "scale = $(spec.scale)")
    isnothing(spec.fn) || push!(opts, "fn = $(nameof(spec.fn))")
    isnothing(spec.times) ||
        push!(opts, "times = <$(length(spec.times)) times>")
    spec.atend isa RepeatAtEnd || push!(opts, "atend = $(spec.atend)")
    isnothing(spec.calendar) || push!(opts, "calendar = $(spec.calendar)")
    return opts
end

function _showspec(io::IO, spec::RasterSpec{A}, ::Nothing) where {A}
    kw = ["$(k) = $(v)" for (k, v) in pairs(spec.readkw)]
    return print(io, "SourceSpec($(spec.source), $(repr(spec.code))",
                 join(", " .* vcat(kw, _readoptions(spec))),
                 ", axis = $(nameof(A)))")
end

# A spec naming files is spelled `RasterFileSpec` where it names no layer, and `SourceSpec` with a
# location keyword where the catalogue describes what it reads.
function _showspec(io::IO, spec::RasterSpec{A},
                   files::AbstractVector) where {A}
    where_ = length(files) == 1 ? "file = " * repr(_pathtext(only(files))) :
             "files = " * repr(_pathtext.(files))
    if isnothing(spec.code)
        unit = spec.unit === NoUnits ? "" : ", unit = $(spec.unit)"
        source = spec.source === SyntheticData ? "" :
                 ", source = $(spec.source)"
        path = length(files) == 1 ? repr(_pathtext(only(files))) : where_
        return print(io, "RasterFileSpec(", path, unit, source,
                     join(", " .* _readoptions(spec)), ", axis = $(nameof(A)))")
    end
    kw = ["$(k) = $(v)" for (k, v) in pairs(spec.readkw)]
    return print(io, "SourceSpec($(spec.source), $(repr(spec.code)), ", where_,
                 join(", " .* vcat(kw, _readoptions(spec))),
                 ", axis = $(nameof(A)))")
end

# --- Desugaring and labelling a spec -----------------------------------------
# A per-cell supply is rewritten into a combine here, before any grid is decided.

# How a layer is named in that message: the dataset and the code asked for, or the file named.
_speclabel(spec::RasterSpec) = _speclabel(spec, spec.files)

function _speclabel(spec::RasterSpec, ::Nothing)
    return "`$(spec.source)` layer `$(spec.code)`"
end

function _speclabel(::RasterSpec, files::AbstractVector)
    return "file `$(_pathtext(first(files)))`"
end

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
function _desugarsupply(spec::RasterSpec)
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
function _percellrecord(spec::RasterSpec)
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
function _percellaxis(spec::RasterSpec, rec::LayerRecord)
    rec.axis <: WaterAxis && return Precipitation
    return error("`$(spec.code)` accumulates over the `$(rec.period.code)` layer, so as a supply it " *
                 "is a rate - but its axis `$(nameof(rec.axis))` has no rate reading defined here. " *
                 "Add one to `_percellaxis` naming the resource it supplies.")
end
