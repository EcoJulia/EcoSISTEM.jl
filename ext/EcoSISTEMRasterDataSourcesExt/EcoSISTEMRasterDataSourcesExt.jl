# SPDX-License-Identifier: LGPL-3.0-or-later

# EcoSISTEMRasterDataSourcesExt
#
# Everything in EcoSISTEM that genuinely needs `RasterDataSources`: the hooks that teach the general
# layer machinery what a dataset is, the two read entry points that download, the land-cover
# operations, the monthly-climate plot recipes, and a home for deprecations whose signatures name a
# dataset (`deprecations.jl`, empty at present).
#
# **This module is a manifest, as `EcoSISTEMMPIExt` is** - the parts are ordinary files it
# includes, laid out the way `src/EcoSISTEM.jl` is, because one file would be unreadable.
#
# **No docstring belongs in here.** `docs/src/api.md`'s `@autodocs` cannot see inside an extension,
# so every public name whose implementation moved keeps its docstring on a method-less stub in the
# parent.
module EcoSISTEMRasterDataSourcesExt

using RasterDataSources
const RDS = RasterDataSources

using EcoSISTEM
using EcoSISTEM.Units
using EcoSISTEM.Units: _monthindex
using Dates: Dates
# Named explicitly because they are `public` or private rather than exported - `using` brings in
# only the exported names, and a moved body calls its neighbours by bare name.
using EcoSISTEM: ClimateRaster, CODE_TYPE, NicheAxis, Supply,
                 _cellareas, _cancel, landcoverclass
using EcoSISTEM: SourceSpec
using EcoSISTEM: compress_landcover, layerinfo, sourcecrs
# The dataset-reading pipeline lives in `src/datasetread.jl`, so these are the parent's.
using EcoSISTEM: _filelist, _firstfile, _getrasterkw, _isblankcrs

using DimensionalData
using DimensionalData: dims, At
using Rasters
using Rasters: Raster, X, Y, Ti
using RecipesBase
using SimpleTraits
using Unitful
using Unitful.DefaultSymbols

# Point `RasterDataSources` at its own subdirectory of EcoSISTEM's scratch space, keeping downloads
# under that scratch lifecycle - unless the user has already set `RASTERDATASOURCES_PATH`.
#
# **Here rather than in the parent's `__init__`**: the variable configures `RasterDataSources`
# itself, so it is meaningless until that package is loaded, and setting it from the parent would be
# the one remaining reason for the parent to name that package at all.
function __init__()
    if !haskey(ENV, "RASTERDATASOURCES_PATH")
        ENV["RASTERDATASOURCES_PATH"] = _SETPATH[] = EcoSISTEM.assetdir(RasterDataSources)
    end
    return nothing
end

# The value `__init__` gave `RASTERDATASOURCES_PATH`, or `nothing` where it was set already, so
# that moving the asset cache moves only a path this extension chose.
const _SETPATH = Ref{Union{Nothing, String}}(nothing)

function EcoSISTEM._repointrasterdatasources()
    (isnothing(_SETPATH[]) ||
     get(ENV, "RASTERDATASOURCES_PATH", nothing) != _SETPATH[]) &&
        return nothing
    ENV["RASTERDATASOURCES_PATH"] = _SETPATH[] = EcoSISTEM.assetdir(RasterDataSources)
    return nothing
end

include("hooks.jl")
include("read.jl")
include("landcover.jl")
include("supply.jl")
include("recipes.jl")
include("deprecations.jl")

end
