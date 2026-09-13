# SPDX-License-Identifier: LGPL-3.0-or-later

# --- What `RasterDataSources` contributes to the general layer machinery --------------------------
#
# `src/Spec.jl` and `src/LazySpec.jl` define what a raster and a spec *are*, and deliberately name no dataset
# package; `src/LayerCatalogue.jl` holds the shipped catalogue, keyed on a bare `Type`. These
# are the methods that teach both of them about `RasterDataSources`. Every one is the **sole** method
# for its hook - the parent declares a fallback and nothing more, so this extension *adds* rather
# than overwrites.

# Marking the abstract type covers the whole `RasterDataSources` hierarchy in one line.
SimpleTraits.@traitimpl EcoSISTEM.IsRasterData{RDS.RasterDataSource}

# Derived from the source package rather than declared here, so it cannot drift from it: an `Int`
# for BioClim and land cover, a `Symbol` for the rest.
function EcoSISTEM._codetype(::Type{T}) where {T <: RDS.RasterDataSource}
    return eltype(RDS.layers(T))
end

# Every layer this dataset has, in the source package's own order - what `SourceSpec(dataset)` (no
# code) expands to. Read from `RasterDataSources` rather than from the shipped table for the same
# reason as `_codetype`: the source package is the authority on what it can fetch, and the catalogue
# is an index into it.
function EcoSISTEM._alllayercodes(::Type{T}) where {T <: RDS.RasterDataSource}
    return collect(EcoSISTEM.CODE_TYPE, RDS.layers(T))
end

# Resolve any spelling a caller may write - `4`, `:bio4`, `"bio4"` - to the one this dataset uses.
# The `::Nothing` method is repeated here rather than inherited: `(::Type{<:RDS...}, ::Any)` and the
# parent's `(::Type{S}, ::Nothing)` are equally specific, so without it an uncoded raster of a
# catalogued source is an ambiguity error.
EcoSISTEM._preferredcode(::Type{<:RDS.RasterDataSource}, ::Nothing) = nothing
function EcoSISTEM._preferredcode(::Type{S},
                                  codes::AbstractVector) where {S <:
                                                                RDS.RasterDataSource}
    return [EcoSISTEM._preferredcode(S, c) for c in codes]
end
function EcoSISTEM._preferredcode(::Type{S},
                                  code) where {S <: RDS.RasterDataSource}
    rec = layerinfo(S, code)     # throws, naming the dataset, when the code is unknown
    for l in RDS.layers(S)
        string(l) in rec.aliases && return l::EcoSISTEM._codetype(S)
    end
    # Reachable only if the shipped table and `RasterDataSources` disagree about a layer's names,
    # which is a packaging error rather than a user one - so it says so.
    return error("`$(repr(code))` is in the shipped table for `$S` but matches none of its " *
                 "`RasterDataSources.layers`; the table and the source package disagree.")
end

# What `_parselayers` dispatches through to turn a bare dataset (optionally with codes) into specs.
EcoSISTEM._isdatasettype(::Type{<:RDS.RasterDataSource}) = true
function EcoSISTEM._datasetspec(dataset::Type{<:RDS.RasterDataSource},
                                ::Nothing)
    return SourceSpec(dataset)
end
function EcoSISTEM._datasetspec(dataset::Type{<:RDS.RasterDataSource}, code)
    return SourceSpec(dataset, code)
end

# --- The one per-dataset fact that is about `getraster` rather than about the data --------------
#
# Everything the readers need to know about a dataset - which axis its files stack on, the
# ceiling behind a published-scale check, its CRS and extent - is in the shipped catalogue, keyed
# on the type alone, so the parent answers it. What remains here is a fact about this package's
# `getraster`: WorldClim monthly climate must name its months, since `getraster` has no default for
# `month` and a read without one is an `UndefKeywordError`.
EcoSISTEM._getrasterkw(::Type{<:RDS.WorldClim{RDS.Climate}}) = (month = 1:12,)
