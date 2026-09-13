# SPDX-License-Identifier: LGPL-3.0-or-later

# --- What RasterDataSources answers about a dataset's files ---------------------------------------
#
# The three hooks whose sole methods need `getraster`: where a dataset's files are (`_fetchfiles`,
# which downloads them on first use), the first of them opened lazily (`_lazysource`), and its CRS
# (`sourcecrs`). Reading is the parent's - `read(::RasterSpec)` in `src/rasters.jl` - and never
# names this package; a spec's files reach it through `_fetchfiles`.

# **Documented on the stub in `src/datasetread.jl`, not here** - and it was briefly
# documented in *both*, which rendered two `sourcecrs` entries with the same HTML anchor. A name that
# can take a parent stub must be documented there and only there.
function EcoSISTEM.sourcecrs(T::Type{<:RDS.RasterDataSource},
                             layers = RDS.layers(T);
                             cut = nothing, scale = nothing, fn = nothing,
                             kw...)
    c = Rasters.crs(EcoSISTEM._lazysource(T, layers; kw...))
    return _isblankcrs(c) ? nothing : c
end

# The first file the source resolves to, opened lazily. `nothing` for `layers` means the whole
# dataset, as it does for a spec without a code. The read options are accepted and ignored, so a
# spec's stored keywords can be splatted in unchanged.
function EcoSISTEM._lazysource(T::Type{<:RDS.RasterDataSource},
                               layers = RDS.layers(T);
                               cut = nothing, scale = nothing, fn = nothing,
                               kw...)
    raw = getraster(T, something(layers, RDS.layers(T));
                    _getrasterkw(T)..., kw...)
    return Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        return Raster(_firstfile(raw); lazy = true)
    end
end

# Where a dataset's layers are, through `getraster`, which downloads them on first use.
function EcoSISTEM._fetchfiles(T::Type{<:RDS.RasterDataSource}, code; kw...)
    return getraster(T, something(code, RDS.layers(T)); kw...)
end
