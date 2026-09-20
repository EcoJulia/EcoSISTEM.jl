# SPDX-License-Identifier: LGPL-3.0-or-later

# --- Deprecations whose signatures name a dataset --------------------------------------------------
#
# A method that dispatches on a `RasterDataSources` type cannot be defined without that package
# loaded, so these are here rather than in `src/deprecations.jl`, labelled and ordered as that file
# is.

using Statistics: mean

# ---------------------------------------------------------------------------
# `read` on a dataset type -> `read(::RasterSpec)`
#
# `Base.read` on RasterDataSources' types is a method on a generic and types this package owns
# neither of; the owned spelling is `read(SourceSpec(T, layers; ...))`. The two methods keep their
# old behaviour for one release - bare magnitudes, a whole dataset of differing units in one array
# - and warn, since the spec form attaches the layer's unit and refuses a mixed-unit stack. They
# remain piracies until they go, which `test/clean_Aqua.jl` names.
#
# Deprecated in v0.8.0.
# ---------------------------------------------------------------------------
function Base.read(T::Type{<:RDS.RasterDataSource}, layers = RDS.layers(T);
                   cut = nothing, scale = 1, fn = _defaultfn(T),
                   axis = _readaxis(T, layers), kw...)
    Base.depwarn("`read($T, layers; ...)` is deprecated; read " *
                 "`SourceSpec($T, layers; ...)` instead, which attaches the layer's unit.",
                 :read)
    rasterkw = (; _getrasterkw(T)..., kw...)
    raw = getraster(T, layers; rasterkw...)
    out = _readraw(T, raw; cut = cut, scale = scale, fn = fn,
                   slices = get(rasterkw, :month, nothing), axis = axis)
    return _rescalepublished(T, layers, out)
end

function Base.read(T::Type{RDS.CHELSA{RDS.Climate}}, dir::AbstractString,
                   var_name::AbstractString; scale = 1, fn = mean,
                   cut = nothing)
    Base.depwarn("`read(CHELSA{Climate}, dir, var_name; ...)` is deprecated; read " *
                 "`SourceSpec(CHELSA{Climate}, var_name, directory = dir; ...)` instead.",
                 :read)
    u = layerunit(T, var_name)
    return ClimateRaster(T,
                         _readmonthlydir(dir, u; scale = scale, fn = fn,
                                         cut = cut))
end
