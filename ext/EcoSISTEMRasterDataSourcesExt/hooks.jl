# SPDX-License-Identifier: LGPL-3.0-or-later

# --- What `RasterDataSources` contributes to the general layer machinery --------------------------
#
# `src/Spec.jl` and `src/LazySpec.jl` define what a raster and a spec *are*, and deliberately name no dataset
# package; `src/LayerCatalogue.jl` holds the shipped catalogue, keyed on a bare `Type`. These
# are the methods that teach both of them about `RasterDataSources`. Every one is the **sole** method
# for its hook — the parent declares a fallback and nothing more, so this extension *adds* rather
# than overwrites.

# Marking the abstract type covers the whole `RasterDataSources` hierarchy in one line.
SimpleTraits.@traitimpl EcoSISTEM.IsRasterData{RDS.RasterDataSource}

# Derived from the source package rather than declared here, so it cannot drift from it: an `Int`
# for BioClim and land cover, a `Symbol` for the rest.
function EcoSISTEM._codetype(::Type{T}) where {T <: RDS.RasterDataSource}
    return eltype(RDS.layers(T))
end

# Every layer this dataset has, in the source package's own order — what `SourceSpec(dataset)` (no
# code) expands to. Read from `RasterDataSources` rather than from the shipped table for the same
# reason as `_codetype`: the source package is the authority on what it can fetch, and the catalogue
# is an index into it.
function EcoSISTEM._alllayercodes(::Type{T}) where {T <: RDS.RasterDataSource}
    return collect(EcoSISTEM.CODE_TYPE, RDS.layers(T))
end

# Resolve any spelling a caller may write — `4`, `:bio4`, `"bio4"` — to the one this dataset uses.
# The `::Nothing` method is repeated here rather than inherited: `(::Type{<:RDS…}, ::Any)` and the
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
    # which is a packaging error rather than a user one — so it says so.
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

# --- Catalogue: the one entry keyed on a concrete dataset ------------------------------------------
#
# The rest of `LayerCatalogue.jl` stays in the parent, because it uses its `Type` argument only to find
# a shipped CSV. This is the single exception — a per-dataset correction, and therefore a fact about
# `RasterDataSources` rather than about the catalogue's shape. The parent's `::Type` fallback of 1.0
# stays where it is; this adds to it.
function EcoSISTEM._documentedceiling(::Type{<:RDS.WorldClim{<:RDS.BioClim}},
                                      code)
    return 100.0
end

# --- Per-dataset facts the readers dispatch on ------------------------------------------------------
#
# **These are the specialisations of `datasetread.jl`'s fallbacks, and two of the three fail SILENTLY
# if they go missing** — which is exactly what happened while this extension was being built: they
# were deleted from the parent with only a comment saying they had moved here, and `extras_canonical`
# (the one gate that reads real data) is what caught it.

# Kind of axis stacked when combining multiple files of one layer: a monthly series (`Ti`) rather
# than the `Dim{:layer}` band index the parent defaults to.
# CHELSA's `Climate` layers are the same shape as WorldClim's — one file per month — so they stack
# on `Ti` too. Without this method they fall through to the `Dim{:layer}` fallback and are read as
# twelve unrelated *bands*, which means a `SourceSpec(CHELSA{Climate}, …)` is **not recognised as a
# series at all** and cannot drive a layer through time. It fails silently: the data loads, it
# simply is not a time series.
EcoSISTEM._stackaxis(::Type{<:RDS.WorldClim{RDS.Climate}}) = Ti
EcoSISTEM._stackaxis(::Type{<:RDS.CHELSA{RDS.Climate}}) = Ti

# Default read-time block-aggregation factor: land cover is coarsened 10×.
# Also silent if lost — the read succeeds at full resolution and every downstream number moves.
EcoSISTEM._defaultscale(::Type{<:RDS.EarthEnv{<:RDS.LandCover}}) = 10

# Default keywords forwarded to `getraster`: WorldClim monthly climate must name its months.
# The one of the three that fails loudly — `getraster` has no default for `month`, so a read
# without this is an `UndefKeywordError` rather than a wrong answer.
EcoSISTEM._getrasterkw(::Type{<:RDS.WorldClim{RDS.Climate}}) = (month = 1:12,)
