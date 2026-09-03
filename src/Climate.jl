# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a raster of climate data IS - the values, the layer code identifying them, and the
# broadcasting that lets a raster be combined like its values and stay a raster.

using DimensionalData

using SimpleTraits

using Unitful

using Dates: Dates

using RecipesBase

using DimensionalData: At, Ti

using EcoSISTEM.Units

using EcoSISTEM.Units: _monthindex

import Base: size, length, eltype

"""
    AbstractClimate

Abstract supertype of all climate data.
"""
abstract type AbstractClimate end

"""
    ClimateRaster{Source, Code, Array} <: AbstractClimate

Type for climate data derived from `RasterDataSource`s.

# Fields

`array` is the grid itself, a `DimensionalData` array over `(Y, X)` - plus a third dimension for a
multi-layer or monthly read. **Prefer to operate on the raster rather than reaching for this**: a
raster broadcasts and yields a raster, so `lc .!= code` and `sum(bands)` work directly, which is what
lets a [`ConstructedRasterSpec`](@ref) combine avoid naming an array type at all.

`code` is the layer this holds, when it is one identifiable layer of its source - or `nothing` for a
whole-dataset read or a raster derived by arbitrary arithmetic. It is what lets a *materialised*
layer still be looked up in the shipped table, which `iscategorical` needs to know whether
resampling it may interpolate: the answer lives in the `ValueType` column, per layer, and a source
type alone cannot supply it (BioClimPlus holds all three value types at once).

**Write it however reads best - `4`, `:bio4` or `"bio4"` - and it is stored as the one spelling
that dataset prefers**, via `_preferredcode`. So two rasters of one layer always
compare equal, and a code naming no layer is refused here rather than surfacing later.

**There is deliberately no `valuetype` field.** Whether the values are class codes or measurements
is a property of the layer's **niche axis**, and the axis lives on the spec that declares what the
layer means - see `iscategorical`. A copy stored here could contradict it, which is the one thing a
declaration must not be able to do.
"""
struct ClimateRaster{S, C, A <: DimensionalData.AbstractDimArray} <:
       AbstractClimate
    array::A
    code::C
    # The sole constructor, and the only door: defining it suppresses Julia's generated
    # parameterised constructors, so `ClimateRaster{SomeUnmarkedType, ...}(...)` has no method. That
    # matters because a trait cannot bound a struct's parameter (`struct T{S; Tr{S}}` is a syntax
    # error), so `S` is free *in the type* and the guarantee lives here instead.
    @traitfn function ClimateRaster(::Type{S}, a::A,
                                    code::C = nothing) where {A <:
                                                              DimensionalData.AbstractDimArray,
                                                              C, S;
                                                              RasterDataAcceptableCode{S,
                                                                                       C}}
        preferred = _preferredcode(S, code)
        return new{S, typeof(preferred), A}(a, preferred)
    end
end

# --- Broadcasting: a raster behaves like its values and stays a raster -------
#
# **This is what lets a `ConstructedRasterSpec` combine be written without naming the array type at
# all** - `compress_landcover(lc) .!= landcoverclass(:open_water)` rather than
# `ClimateRaster(EarthEnv{LandCover}, compress_landcover(lc).array .!= ...)`. A combine is *user* code,
# and the array type is an implementation detail we intend to stay free to change; requiring it to be
# named on the way in and constructed on the way out made every such change everyone's problem.
#
# The wrapper adds **no arithmetic of its own**: every raster in the expression is replaced by its
# array and the work is done by DimensionalData exactly as before, so a broadcast means precisely
# what it meant when it was written out by hand.
struct ClimateRasterStyle <: Broadcast.BroadcastStyle end

# == Functions ==================================================================================

# --- What a layer IS, and how one is declared ----------------------------------------------------
#
# **Included very early**, so what it declares constrains what it may use: `NicheAxis`
# (`Ecology.jl`), `SimpleTraits` and `DimensionalData` are all it has. Everything that *reads* a
# raster - `SourceSpec` (`LazySpec.jl`), the readers (`datasetread.jl`, `erareaders.jl`), the shipped
# tables (`LayerCatalogue.jl`) - comes later and constructs the types declared here.
#
# **Nothing here names `RasterDataSources`**, which is a weak dependency, and widening
# `ClimateRaster`'s source parameter onto the `IsRasterData` trait is what makes that possible. The
# dataset-shaped hooks below - `_codetype`, `_preferredcode`, `_isdatasettype` - are declared here
# and answered by `EcoSISTEMRasterDataSourcesExt`.

using DimensionalData

using Unitful

using Dates: Dates

using SimpleTraits

using RecipesBase

# --- Display ------------------------------------------------------------------
# A raster carries a whole array, so without this the default prints it: measured at 1 757
# characters for a 5 x 5 synthetic one, and unbounded in the grid.
#
# The source and the layer code are what identify a raster - they are the two things that decide
# what its values *mean* - so both are on the line, with the code omitted where there is none.
# Nothing here touches the values.
# The source is rendered with its parameters intact and its module qualifier dropped. `nameof` will
# not do: it strips parameters, which collapses `DerivedData{EarthEnv{LandCover}}` to `DerivedData`
# and so throws away exactly the lineage that type exists to record.
_sourcename(::Type{S}) where {S} = replace(string(S), "EcoSISTEM." => "")

function Base.show(io::IO, r::ClimateRaster{S}) where {S}
    dims = join(size(r.array), " × ")
    code = isnothing(r.code) ? "" : ", $(repr(r.code))"
    return print(io, "ClimateRaster{$(_sourcename(S))}($(dims)$(code))")
end

function Base.show(io::IO, ::MIME"text/plain", r::ClimateRaster{S}) where {S}
    println(io, sprint(show, r))
    println(io, "  source    ", S)
    println(io, "  code      ", isnothing(r.code) ? "none" : repr(r.code))
    print(io, "  array     ", join(size(r.array), " × "), " ", eltype(r.array))
    return nothing
end

# --- The netCDF climate sources ----------------------------------------------
# `ERA`, `CERA` and `CRUTS` are **data sources**, not containers: a reader returns a
# `ClimateRaster{ERA}`, exactly as it returns a `ClimateRaster{WorldClim{BioClim}}` for a
# RasterDataSources dataset. They are `EcoSISTEMSource` subtypes because the archives are ours to
# describe rather than `RasterDataSources`', and that supertype is what grants `IsRasterData`
# (`@traitimpl IsRasterData{EcoSISTEMSource}` below), so each can name a source in a `SourceSpec`.
#
# They stay here rather than in `DataSource.jl` because that file holds the *generic* markers -
# "synthetic", "derived from" - while these name particular archives, which is climate data's
# business.

"""
    ERA <: EcoSISTEMSource

The ERA5 reanalysis archive, as a data source. Read one with
`read(ERA, file, param)`, which returns a [`ClimateRaster`](@ref)`{ERA}`.

Fieldless: a source names *where data came from*, and the data itself lives in the
[`ClimateRaster`](@ref) that carries it.
"""
struct ERA <: EcoSISTEMSource end

"""
    CERA <: EcoSISTEMSource

The CERA-20C reanalysis archive, as a data source. Read one with
`read(CERA, dir, file, param)`, which returns a [`ClimateRaster`](@ref)`{CERA}` - the archive is one
file per decade, and the reader concatenates them along time.

Fieldless, as [`ERA`](@ref) is.
"""
struct CERA <: EcoSISTEMSource end

"""
    CRUTS <: EcoSISTEMSource

The CRU TS archive, as a data source. Read one with `read(CRUTS, dir, var_name)`, which returns a
[`ClimateRaster`](@ref)`{CRUTS}`.

CRU TS has no layer table of its own, so its variable codes and units are taken from
`WorldClim{Climate}`'s.

Fieldless, as [`ERA`](@ref) is.
"""
struct CRUTS <: EcoSISTEMSource end

# The time-axis guard the three container types used to carry in their inner constructors, kept as a
# function because `ClimateRaster` has no such check and the readers are the only door. Every reader
# and plot recipe for these sources slices by time, so a third dimension that is not time is a
# failure later and further away.
function _timeseriesraster(::Type{S},
                           array::DimensionalData.AbstractDimArray) where {S}
    _istimeaxis(eltype(DimensionalData.lookup(array, 3))) ||
        error("Third dimension of array must be time")
    return ClimateRaster(S, array)
end

"""
    in_memory_raster(raster::ClimateRaster; axis)

Wrap a raster you already hold as a layer spec, declaring what its values mean.

**Prefer naming the source.** A [`SourceSpec`](@ref) lets the package read only the window the
study area needs, cache the result between layers, and take the layer's unit, axis, accumulation
period and value type from the shipped catalogue. An in-memory raster gives all of that up: it is
read in full by whoever built it, cached nowhere, and describes itself only by the `axis` given
here. Reach for this when the data genuinely did not come from a catalogued source - something
computed elsewhere, or read by hand - not as the ordinary way to build a layer.

**It exists because a raster is refused as a spec, and rightly so**: a raster carries values and
possibly a layer code, but no niche axis, so nothing about it says whether those numbers are a
temperature, a rainfall rate or a cover fraction. The spec is where that is declared, which is
exactly what a raster cannot do - so this is the pathway, and `axis` is the whole of what it adds.

# Arguments

  - `raster`: the [`ClimateRaster`](@ref) to wrap. It is returned by the spec's combine verbatim, so
    it is neither re-read nor re-projected before the study area samples it.
  - `axis`: the [`NicheAxis`](@ref) the values are on - what makes them matchable against a species'
    tolerances. Required: pass `NicheAxis` itself for data whose meaning is not being claimed, but a
    layer meant to pair with a tolerance needs a real one.
"""
function in_memory_raster(raster::ClimateRaster;
                          axis::Type{<:NicheAxis})
    return ConstructedRasterSpec(() -> raster, axis = axis)
end

# --- Which types may name a raster's source, and how that source names its layers ----------------
#
# **Three questions, deliberately kept apart.** `IsRasterData` asks *may this type name a source at
# all*; `RasterDataAcceptableCode` asks *is this a plausible shape for a layer name*; and
# `_preferredcode` asks *which layer is that, and what does this dataset call it*. Conflating
# them is what the old `S <: RDS.RasterDataSource` bound did: one `<:` answered the first two by
# accident and could not ask the third at all, so a raster with **no** dataset - a derived layer, a
# synthetic field, an ECMWF download - had to name a dataset it did not come from, and a code naming
# no layer was accepted without complaint.
#
# **A struct parameter cannot carry a trait bound** (`struct Bad{S; IsRasterData{S}}` is a syntax
# error), so `ClimateRaster`'s `S` is unconstrained *in the type*. The guarantee comes from the single
# inner constructor instead: it is the only door, it dispatches on the trait, and defining it
# suppresses Julia's generated parameterised constructors - so `ClimateRaster{SomethingElse, ...}` is
# expressible as a type but **unconstructible**.

# `RasterDataSources`' own hierarchy is marked in `EcoSISTEMRasterDataSourcesExt`, which is the only
# place that package is visible. This file marks only what it defines itself.

@traitimpl IsRasterData{EcoSISTEMSource}

@traitimpl RasterDataAcceptableCode{S, C} < - _acceptablecode(S, C)

# _derivedfrom(source)
#
# Return the source that a layer derived from `source` should carry - `DerivedData{source}`, or
# `source` itself when it is already derived.
#
# **This is how the collapse is actually enforced, and it must be at the *type* level.** A
# `DerivedData` is never instantiated: it appears only as [`ClimateRaster`](@ref)'s first type
# parameter, so a normalising *constructor* would never run and
# `DerivedData{DerivedData{DerivedData{...}}}` would accumulate one wrapper per operation - silently,
# since each nesting is a perfectly good type. Measured before this existed: three chained
# broadcasts gave three levels.
#
# Idempotent by construction, so it needs no recursion: the first method catches an
# already-derived source, and nothing else can produce one.
# The sources a lineage records, flattened out of any `DerivedData`/`Tuple` wrapper.
_origins(::Type{DerivedData{S}}) where {S} = _origins(S)

_origins(::Type{S}) where {S <: Tuple} = Tuple(S.parameters)

_origins(::Type{S}) where {S} = (S,)

function _derivedfrom(sources...)
    isempty(sources) && return DerivedData{SyntheticData}
    for s in sources
        _checkderivable(s)
    end
    # **Sorted, and that is not cosmetic**: without it `lc .+ cl` and `cl .+ lc` would be
    # different types for the same computation - the argument-order dependence this whole change
    # exists to remove, reappearing one level up. Deduplicated so re-deriving from a source already
    # in the lineage does not grow it.
    origins = sort(collect(Set(Iterators.flatten(_origins(s) for s in sources))),
                   by = string)
    return length(origins) == 1 ? DerivedData{first(origins)} :
           DerivedData{Tuple{origins...}}
end

@traitfn _checkderivable(::Type{S}) where {S; IsRasterData{S}} = nothing

# Deriving from something that could not have been a source in the first place is a mistake worth
# naming, rather than a `DerivedData` of a type nothing else will accept.
@traitfn function _checkderivable(::Type{S}) where {S; !IsRasterData{S}}
    return error("`$S` is not raster data, so nothing can be derived from it. If it names a " *
                 "raster source, mark it with `@traitimpl IsRasterData{$S}`.")
end

# Worth the extra method: without it the failure is a bare `MethodError` mentioning
# `Not{RasterDataAcceptableCode{...}}`, which leaks SimpleTraits' internals and names no remedy.
@traitfn function ClimateRaster(::Type{S}, a,
                                code::C) where {C, S;
                                                !RasterDataAcceptableCode{S,
                                                                          C}}
    return error("a layer of `$S` is named by a code, a vector of codes, or `nothing`; " *
                 "got a `$C`.")
end

# _codetype(source)
#
# Return the type this source's layer codes are stored as, or `Nothing` if it has none.
#
# **Derived from the source package rather than decreed here**, so it cannot drift from it:
# `eltype(RasterDataSources.layers(T))`, which is an `Int` for BioClim and land cover and a `Symbol`
# for the rest. **It is therefore per-dataset, not one type for all of them** - the alternative,
# canonicalising every dataset to a `Symbol`, would disagree with `RasterDataSources` (and with what
# this package already stores) for three of its eight datasets.
#
# `Nothing` for every source with no layers of its own to name - which is **both** a synthetic field
# and anything derived. A derived raster is not the layer it came from: the combine is free to change
# what the values *are*, so inheriting the parent's code would attach that layer's whole catalogue row
# - its unit, axis and accumulation period - to a quantity none of them describe. `iscategorical` says
# the same thing: *"a raster with no code at all is different and legitimate - every synthetic or
# derived layer is one"*.
#
# **`_codetype(S) === Nothing` is also the "has no catalogue" test**, since a source with no layer
# codes has no row to look up - one fact, not two.
# **Codeless unless a source says otherwise.** A source this file has never heard of - a foreign
# type marked with `IsRasterData` - has no layers we can name, which is the honest default; it says
# otherwise with a `_codetype` method of its own. The `RasterDataSources` method lives in
# `EcoSISTEMRasterDataSourcesExt`, which is the only place that package is visible.
_codetype(::Type) = Nothing

# _alllayercodes(source)
#
# Every layer code `source` has, as a `Vector{CODE_TYPE}` - what a whole-dataset `SourceSpec(dataset)`
# expands to, so each layer's identity (and therefore its own unit) is known before anything is read.
#
# The exact sibling of `_codetype` above, and for the same reason: the answer is the source
# package's to give, not ours to decree. The `RasterDataSources` method is
# `collect(CODE_TYPE, RasterDataSources.layers(T))` and lives in the extension beside it.
# A source with no layers of its own to name - synthetic, derived, or a foreign type marked with
# `IsRasterData` but not taught about - cannot answer, and says so rather than returning an empty
# list, which would silently describe a dataset with nothing in it.
function _alllayercodes(::Type{S}) where {S}
    return error("`$S` does not say what layers it has, so a whole-dataset spec cannot be " *
                 "expanded; name the layer you want, or give `$S` an `_alllayercodes` method.")
end

# Is this code *shape* plausible for a raster of source `S` - one code, a list of them, or none?
# `Nothing` is always allowed: even a catalogued source produces uncoded rasters, e.g. a
# whole-dataset read whose layers are stacked rather than identified singly.
#
# **Shape only, deliberately.** Whether a particular value names a real layer is a question about
# the catalogue, not about types, and belongs to `_rasterdatapreferredcode` below. Trying to make a
# type-level test carry it is what forces a hand-written list of admissible types that still cannot
# reject `:not_a_layer`. A bare `AbstractVector` is required rather than
# `Vector{Symbol}`/`Vector{Int}`: mixing spellings (`[:bio4, 3]`) gives a `Vector{Any}`, which is
# legitimate and is normalised element by element.
function _acceptablecode(S, C)
    return C === Nothing ||
           C <: Union{Symbol, Integer, AbstractString} ||
           C <: AbstractVector
end

# _preferredcode(source, code)
#
# Return `code` as `source`'s own preferred spelling, or throw if it names no layer of that source.
#
# **`RasterDataSources` is the authority, and the catalogue is the index into it.** The catalogue
# records every spelling a layer answers to (`bio4` is `["4", "bio4"]`; the EarthEnv texture measures
# are `["Contrast", "contrast"]`), and `RasterDataSources.layers` says which of those the package
# itself uses - `4` for BioClim, `:Contrast` for habitat heterogeneity. So the preferred code is the
# alias that appears in `layers`, and it is neither uniformly an `Int` nor uniformly a `Symbol`: it is
# whatever that dataset natively calls its layers.
#
# **Every construction pays the lookup, on purpose.** Skipping it for an already-canonical code
# would also skip the validation, and reading data is not a hot path - so a code that names nothing is
# refused where it is written, rather than surfacing later as an empty read.
function _preferredcode end

_preferredcode(::Type{S}, ::Nothing) where {S} = nothing

function _preferredcode(::Type{S}, codes::AbstractVector) where {S}
    return [_preferredcode(S, c) for c in codes]
end

# The `::_codetype(S)` annotation is what keeps the constructor inferable: without it the stored
# parameter depends on a value the compiler cannot see, and `ClimateRaster`'s return type is `Any`.
function _preferredcode(::Type{S}, code) where {S}
    _codetype(S) === Nothing &&
        error("a `$S` raster holds no layer codes, so it cannot be given `$(repr(code))`. " *
              "A derived or synthetic layer is identified by its spec's `axis`, not by a code.")
    # Reachable only if a source declares a `_codetype` without a way to resolve a spelling to it.
    return error("`$S` declares its layer codes as `$(_codetype(S))` but supplies no way to " *
                 "resolve one; `_preferredcode` needs a method.")
end

# Replace every raster in a broadcast tree by its array, leaving the tree's shape untouched.
_unwrapped(x) = x

_unwrapped(raster::ClimateRaster) = raster.array

function _unwrapped(bc::Broadcast.Broadcasted)
    return Broadcast.broadcasted(bc.f, map(_unwrapped, bc.args)...)
end

# Every raster in a broadcast tree, in argument order, so the result can be rebuilt from them.
_rastersof(x) = ()

_rastersof(raster::ClimateRaster) = (raster,)

function _rastersof(bc::Broadcast.Broadcasted)
    return reduce((a, b) -> (a..., b...), map(_rastersof, bc.args), init = ())
end

# Rebuild a raster around a broadcast result.
#
# **Dropping a disagreeing `code` is the point, not a shortcoming.** Adding three land-cover bands
# gives a quantity that is none of the three, so inheriting one input's code would attach a claim the
# values cannot support; a `ConstructedRasterSpec`'s own `axis` is what declares what a derived layer means.
# Agreement is the interesting case anyway - masking a single layer keeps its code, which is what
# leaves a derived **mask** identifiable.
_sourceof(::ClimateRaster{S}) where {S} = S

function _rewrap(array, rasters::Tuple{ClimateRaster, Vararg{ClimateRaster}})
    # **Nothing is inherited - not the source, not the code, not the value type.** A combine is
    # free to change all three: `argmax` over twelve continuous land-cover bands gives a *class code*,
    # so even a unanimous `:continuous` among its inputs would be the wrong answer for its output.
    # An earlier version kept `code` and `valuetype` "where every input agrees", which is exactly
    # the reasoning that put a false source on every derived raster - agreement among the inputs says
    # nothing about the output.
    # What a derived layer *is* comes from its spec's declared `axis`; what it came *from* is in
    # `DerivedData`'s parameter.
    # **`DerivedData{S}`, not `S`** - and agreement is the wrong test for the *source*, which is
    # why this differs from the two fields above. Every input to `sum(bands)` agrees it is
    # `EarthEnv{LandCover}`, and the sum is still not land cover; multiply those same bands by an
    # incident flux and the result is solar radiation. What a combine produces is a new quantity whose
    # meaning comes from its spec's declared `axis`, so **every** broadcast result is derived,
    # agreement or not. `DerivedData` keeps the lineage without claiming to be the thing.
    # Resampling is the opposite case and must *keep* the source: it moves the same data onto
    # another grid rather than computing anything - see `_sampledeclared`/`_declare`.
    # **Every input's source, not the first one's.** Taking `first(rasters)` made the result depend
    # on argument order - `lc .+ cl` and `cl .+ lc` disagreed - and silently dropped the other
    # lineages.
    return ClimateRaster(_derivedfrom(map(_sourceof, rasters)...), array)
end

# Is this argument a *dataset type* - one that names a whole catalogued source, and so wants the
# following argument read as its layer code(s)?
#
# **A hook, because the answer needs a weak dependency.** This file cannot name `RasterDataSources`,
# so `EcoSISTEMRasterDataSourcesExt` supplies the sole `RasterDataSource` methods for both of these -
# the same method-less-hook pattern the package uses for `retrieve_era5`.
_isdatasettype(_) = false

# Build the spec for one layer of a dataset. Only ever called when `_isdatasettype` said yes, so the
# fallback exists to give a comprehensible error if something declares the first without the second.
function _datasetspec(dataset, code)
    return error("`$dataset` says it is a dataset type but supplies no way to build a layer spec " *
                 "from it; `_datasetspec` needs a method.")
end

# Normalise `ConstructedRasterSpec`'s trailing layer arguments to a `Vector{AbstractSpec}`.
#
# **Three shapes are accepted**, and the middle one is why this needs a parser rather than a `map`:
# - any `AbstractSpec` - a `SourceSpec`, a nested `ConstructedRasterSpec`, or a **synthetic** spec;
# - a dataset type, optionally followed by its code(s) - a scalar, or a vector/tuple giving one
# single-layer spec per code; a bare dataset means the whole dataset as one multi-band read;
# - nothing else.
#
# The lookahead is what makes `ConstructedRasterSpec(f, EarthEnv{LandCover}, [:a, :b])` legible, and it
# is also why a following argument that is itself a spec or a type must *not* be eaten as a code.
function _parselayers(args...)
    layers = AbstractSpec[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a isa AbstractSpec
            push!(layers, a)
            i += 1
        elseif _isdatasettype(a)
            if i < length(args) && !(args[i + 1] isa Type) &&
               !(args[i + 1] isa AbstractSpec)
                codes = args[i + 1]
                for c in (codes isa Union{AbstractVector, Tuple} ? codes :
                          (codes,))
                    push!(layers, _datasetspec(a, c))
                end
                i += 2
            else
                push!(layers, _datasetspec(a, nothing))
                i += 1
            end
        else
            error("ConstructedRasterSpec layer argument $i: expected a layer spec or a dataset type, " *
                  "got $(typeof(a)).")
        end
    end
    return layers
end

# --- A raster on its own grid ------------------------------------------------

# The combined raster put on `target`. Rewrapped as a `ClimateRaster` because `_sampledata` hands
# back a bare `DimArray`, while everything downstream - `_materialisefield`, and the enclosing
# combine itself - takes the wrapper; its source and code travel with it.
function _sampledeclared(combined::ClimateRaster{S}, target,
                         axis::Type{<:NicheAxis}) where {S}
    return ClimateRaster(S,
                         _sampledata(combined, target, name = "layer",
                                     categorical = iscategorical(combined,
                                                                 axis)),
                         combined.code)
end

# A combine whose result is not a `ClimateRaster` - a bare mask or array - has no grid provenance to
# sample *from*, so the early path cannot place it. The late path can, because there each layer was
# already on the target before the combine ran.
function _sampledeclared(combined, target, axis)
    return error("a `CombineOnSourceGrid` combine must return a `ClimateRaster`, so that its " *
                 "result carries the grid it was computed on and can be sampled onto the study " *
                 "area; this one returned a `$(typeof(combined))`. Wrap it " *
                 "(`ClimateRaster(dataset, array)`), or leave the spec on the default " *
                 "`CombineOnTargetGrid`, where the layers are sampled before they are combined.")
end

# A raster's own `(Y, X)` grid as a `Rasters.Raster` template, carrying its coordinates exactly as
# the raster states them - an exactly no-op resample for the raster itself.
#
# The coordinates keep their unit, like every other grid's, and `_reproject` undresses them for GDAL
# at the one boundary that needs it. Stripping them here would make the template the odd one out.
function _owngrid(raster::ClimateRaster)
    yx = dims(raster.array, (Y, X))
    return Rasters.Raster(zeros(length.(yx)), yx)
end

# --- What a raster or a layer stack IS ----------------------------------------
# One question, asked of whatever the caller happens to be holding: an axis, a dataset and a layer
# code, a stack of codes, or a raster - answered by `iscategorical` throughout. Only the stack method
# can throw, and only because layers that disagree have no answer rather than a `false` one.

# Degrees north (latitude) and east (longitude) of a raster's cell centres.
function _latvals(raster::ClimateRaster)
    return parent(DimensionalData.lookup(raster.array, Y))
end

# A raster's longitude coordinates as a plain vector, unwrapped from its lookup. Used where the
# values are compared or rewrapped rather than indexed, which a `Dimension` cannot do directly.
function _longvals(raster::ClimateRaster)
    return parent(DimensionalData.lookup(raster.array, X))
end

# Does this raster hold class codes rather than measurements? It decides two things, and they are the
# same question: whether the regime built from it is a `CategoricalRegime`, and whether resampling it
# must take the nearest class (`:mode`) instead of interpolating.
#
# Named for *categorical*, not "discrete", because the shipped `ValueType` column distinguishes
# three things and only one of them forbids averaging. Its `discrete` layers - a day-of-year, a count
# of growing-degree days - are ordinary numbers that average perfectly well; only its `categorical`
# ones (`kg0`-`kg5`, the Köppen-Geiger/Wissmann/Thornthwaite/Troll-Pfaffen typologies) are class
# codes whose mean is meaningless. The old name `_isdiscrete` therefore pointed at the wrong half of
# the split. EarthEnv land cover looks like a counterexample and is not: its bands are each a
# per-class continuous % cover, categorical only after a reducer such as `compress_landcover`
# collapses them to one winning class.
#
# **Whether a layer holds class codes is a property of its AXIS**, so this asks `iscategorical`
# and nothing declares it separately. Measured across the shipped catalogue: none of its 33 axes
# carries more than one value type, so a per-layer `valuetype` was always a copy of what the axis
# already said - and a copy that could contradict it.
#
# **Two arities, and the difference between them is the whole point.** Asking about a raster alone
# means "work it out from my own code"; asking with an axis means "the axis is this". A single method
# with an `axis = NicheAxis` default could not tell those apart, because "I was not told" and "the
# axis is `NicheAxis`" would be the same value. Dispatch removes the collision rather than moving it:
# there is no magic value left to test for.
#
# **`NicheAxis` gets its own method**, which is not the same thing as the sentinel above. It is the
# *absence* of a named axis rather than a claim that the values are continuous - `iscategorical`'s
# root fallback is a majority answer, not a declaration - so where real information exists it wins.
# It has to: a `SourceSpec` derives its axis from its codes and falls back to `NicheAxis` when
# they disagree, so an ordinary multi-layer stack such as `SourceSpec(WorldClim{BioClim})` arrives
# here carrying it. Answering `false` would bilinearly interpolate the class codes of any
# all-categorical stack spanning two typologies. A *named* axis stays authoritative.
#
# A raster with no code at all cannot answer from the catalogue, and `false` is right rather than an
# error: nothing about it claims to be class-coded, and a caller who does mean class codes says so -
# `in_memory_raster(raster, axis = LandCoverTypology)`.
function iscategorical(raster::ClimateRaster{S}) where {S}
    isnothing(raster.code) && return false
    return iscategorical(S, raster.code)
end

iscategorical(::ClimateRaster, axis::Type{<:NicheAxis}) = iscategorical(axis)

iscategorical(raster::ClimateRaster, ::Type{NicheAxis}) = iscategorical(raster)

# ---------------------------------------------------------------------------
# The concrete data-source types
# ---------------------------------------------------------------------------
# `ERA`, `CERA` and `CRUTS` are siblings of `ClimateRaster` under `AbstractClimate`, so they live
# beside it; their readers are in `erareaders.jl`. The plot recipes travel with the types they plot,
# as they do for `Ecosystem` and `Layer`.

# Whether a `Ti` axis eltype represents time: either a real calendar coordinate (an anchored source's
# own dates) or a `Unitful.Time` offset (a climatology's ticks, or a caller-supplied override).
#
# **`Dates.TimeType`, never `Dates.AbstractDateTime`** - the narrower name looks right and is not.
# Julia branches `AbstractTime -> TimeType -> {AbstractDateTime -> DateTime, Date, Time}`, so
# `AbstractDateTime` excludes a plain `Date` (a daily series could not be constructed at all) and
# every `CFTime` calendar (`AbstractCFDateTime <: Dates.TimeType`), which is how NCDatasets decodes a
# netCDF whose `calendar` attribute is not Gregorian. `TimeType` is also the level the simulation side
# already uses - `Ecosystem.epoch`, `DatedSeries.start`, `_giventimes` - so this brings the two into
# line rather than widening past them.
function _istimeaxis(::Type{T}) where {T}
    return T <: Unitful.Time || T <: Dates.TimeType
end

# --- Plot recipes -----------------------------------------------------------------------------------
#
# **Here, not in an extension.** A `@recipe` needs only `RecipesBase`, which this package depends
# on hard, so these are defined unconditionally and are simply inert until a plotting backend is
# loaded - the same arrangement `src/Ecosystem.jl` and `src/Layer.jl` already use.
# **An extension whose content names nothing from its own trigger package is a trap**: its trigger
# and its dependencies drift apart, and nothing surfaces it until some unrelated change does. A
# recipe names only `RecipesBase`, so it belongs here.
#
# The monthly-climate recipes (`ClimateRaster{WorldClim{Climate}}`/`{CHELSA{Climate}}`) are in
# `EcoSISTEMRasterDataSourcesExt`, with everything else keyed on a dataset.

# Recipe for plotting ERA and CERA data from a particular time period.
@recipe function f(era::ClimateRaster{<:Union{ERA, CERA}}, time::Unitful.Time)
    tm = ustrip.(uconvert(year, time))
    yr = floor(Int64, tm)
    ind = round(Int64, (tm - yr) / (1 / 12))
    mnth = Dates.monthabbr(ind + 1)
    A = transpose(ustrip.(era.array[Ti(At(time))]))
    x = ustrip.(parent(DimensionalData.lookup(era.array, 1)))
    y = ustrip.(parent(DimensionalData.lookup(era.array, 2)))
    seriestype := :heatmap
    grid --> false
    title --> "$yr $mnth"
    return x, y, A
end

@recipe function f(era::ClimateRaster{<:Union{ERA, CERA}}, time::Unitful.Time,
                   xrange, yrange)
    tm = ustrip.(uconvert(year, time))
    yr = floor(Int64, tm)
    ind = round(Int64, (tm - yr) / (1 / 12))
    mnth = Dates.monthabbr(ind + 1)
    A = transpose(ustrip.(era.array[Y(xrange), X(yrange), Ti(At(time))]))
    step1 = ustrip(parent(DimensionalData.lookup(era.array, 1))[2] -
                   parent(DimensionalData.lookup(era.array, 1))[1])
    step2 = ustrip(parent(DimensionalData.lookup(era.array, 2))[2] -
                   parent(DimensionalData.lookup(era.array, 2))[1])
    x = ustrip(xrange.left):step1:ustrip(xrange.right)
    y = ustrip(yrange.left):step2:ustrip(yrange.right)
    seriestype := :heatmap
    grid --> false
    title --> "$yr $mnth"
    return x, y, A
end

# **Parked, and needing more than this package has.** Uncommenting `getprofile` would need `Plots`
# itself, for `histogram` and `px`, **and** `IndexedTables` - neither of which this package depends
# on, and a `@recipe` cannot supply either. Kept rather than deleted, as `SizeDemand` is, because it
# records an intent; but it cannot be revived without deciding how those two arrive.
#=
"""
    getprofile(spp_names::Vector{String}, data::IndexedTable, variable_name::String, dims::Tuple{Int64, Int64} = (1,1))

Function to plot climate profiles for a vector of species names, `spp_names`, using a JuliaDB table of GBIF records, `data`, column containing climate variable of interest, `var`, and dimensions over which it should be plotted, `dims`.
"""
function getprofile(spp_names::Vector{String}, data::IndexedTable, var::Symbol,
                    label::String, dims::Tuple{Int64, Int64} = (1, 1))
    # Check for dimensions to be greater or the same as the length of species, or for all to be plotted in one plot pane.
    (dims[1] * dims[2] >= length(spp_names) || dims == (1, 1)) ||
        error("Dimensions not big enough for number of species")
    first_plot = true
    # Loop through species names, adding profiles to plot
    hist = histogram(layout = dims)
    for i in eachindex(spp_names)
        # Filter data for species in question
        spp = filter(p -> p[:species] == spp_names[i], data)
        # Select climate variable of interest
        vals = select(spp, var)
        # Remove NaNs
        res = vcat(vals...)
        res = res[.!isnan.(res)]
        # If dimensions are (1,1)...
        # ... all plotted to the same subplot window, else individual for each species
        sp = ifelse(dims == (1, 1), 1, i)
        # ... put legend on the left, otherwise none
        lg = ifelse(dims == (1, 1), :left, :none)
        # ... title as species name, else blank
        title = ifelse(dims == (1, 1), "", spp_names[i])
        # For the first species plot histogram, otherwise add to previous plot
        histogram!(res, bins = -20:2:30, grid = false, xlabel = label,
                   label = spp_names[i], subplot = sp, legend = lg,
                   top_margin = 20px,
                   bottom_margin = 20px, title = title)
    end
    return hist
end
=#
