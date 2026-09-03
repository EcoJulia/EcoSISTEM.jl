# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Where a layer's data comes from - the marker hierarchy a spec reads from, and the two
# SimpleTraits markers that say whether a source is raster data and what codes it accepts.

using SimpleTraits

"""
    EcoSISTEMSource

Abstract supertype of the sources this package defines itself, for rasters that come from no
dataset - see [`SyntheticData`](@ref) and [`DerivedData`](@ref).
"""
abstract type EcoSISTEMSource end

"""
    SyntheticData

Source of a raster that was **generated**, not read: a synthetic supply or regime field built from a
spec. It has no dataset, no layer codes and no catalogue row, and says so rather than borrowing a
real dataset's name.
"""
struct SyntheticData <: EcoSISTEMSource end

"""
    DerivedData{S}

Source of a raster **computed from** `S` but no longer `S` - the result of a [`ConstructedRasterSpec`](@ref)
combine, or of an operation like `compress_landcover` that turns a dataset's layers into a quantity
that is none of them. Summing eight land-cover bands does not give land cover, and what the result
*means* comes from the spec's declared `axis` rather than from its inputs, so
`DerivedData{EarthEnv{LandCover}}` records the lineage without claiming to be it.

**Never write this type out - call `_derivedfrom(S)`**, which is what collapses nesting. Deriving
from derived data is still just derived data, and no consumer asks how many steps there were; but
`DerivedData{DerivedData{S}}` stays constructible, so the collapse is that function's job and not
the type's.

`S` is a **phantom** parameter, never instantiated: it appears only as [`ClimateRaster`](@ref)'s
first type argument.
"""
struct DerivedData{S} <: EcoSISTEMSource end

"""
    IsRasterData{X}

Trait marking a type as one that may name a raster's data source.

Marking an abstract type covers its whole subtree, so `RasterDataSources`' hierarchy and this
package's own each take one line. **A third party's raster type needs one `@traitimpl` line and no
change here**, which is the point of a trait rather than a supertype bound.
"""
@traitdef IsRasterData{X}

"""
    RasterDataAcceptableCode{S, C}

Trait holding when `C` is a shape in which a layer of source `S` may be *named* - a scalar code, a
vector of them, or `Nothing`.

**This admits a spelling; it does not confirm a layer exists.** Rejecting `:not_a_layer` needs the
catalogue and happens there. `C` records the *canonical* type a code resolves to, never the shape the
caller wrote - see [`CODE_TYPE`](@ref).
"""
@traitdef RasterDataAcceptableCode{S, C}

"""
    CODE_TYPE

The type of a single `RasterDataSources` layer code **as a caller may write it**: an `Int` layer
number, a `Symbol` key or a `String`. All three name the same layer - `layerinfo(4)`,
`layerinfo(:bio4)` and `layerinfo("bio4")` return one record.

**This is the input vocabulary, not what gets stored.** A [`ClimateRaster`](@ref) resolves whatever
it is given to that dataset's preferred spelling (`_preferredcode`), so two rasters of one layer
cannot disagree about how it is named; [`SourceSpec`](@ref)'s `code` field takes it for the same
reason.
"""
const CODE_TYPE = Union{Symbol, Int, AbstractString}
