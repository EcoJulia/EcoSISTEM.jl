# SPDX-License-Identifier: LGPL-3.0-or-later

# EcoSISTEMRasterDataSourcesExt
#
# Everything in EcoSISTEM that genuinely needs `RasterDataSources`: the hooks that teach the general
# layer machinery what a dataset is, the two read entry points that download, the land-cover
# operations, the monthly-climate plot recipes, and the deprecations whose signatures name a dataset.
#
# **This module is a manifest, as `EcoSISTEMMPIExt` is** - the parts are ordinary files it
# includes, laid out the way `src/EcoSISTEM.jl` is, because one file would be unreadable.
#
# **A `#` comment, not a docstring, deliberately.** `docs/src/api.md` names this module in its
# `@autodocs` list (see below), so a docstring here would be rendered as a public API entry - and
# this is maintainer prose about the layout, not something a user needs.
#
# **Almost no docstring belongs in here either.** Every public name whose implementation moved
# keeps its docstring on a method-less stub in the parent, because `@autodocs` cannot otherwise see
# it. The **one** exception is `Base.read`'s three dataset methods: a stub for them would mean the
# parent defining `read(::Type, ...)`, pirating `Base.read` for every type in Julia. Naming this module
# in `api.md` is what keeps those three in the manual, and it is the only reason it is named there.
module EcoSISTEMRasterDataSourcesExt

using RasterDataSources
const RDS = RasterDataSources

using EcoSISTEM
using EcoSISTEM.Units
using EcoSISTEM.Units: _monthindex
using Dates: Dates
# Named explicitly because they are `public` or private rather than exported - `using` brings in
# only the exported names, and a moved body calls its neighbours by bare name.
using EcoSISTEM: AbstractSupply, ClimateRaster, CODE_TYPE, ConstructedSpec,
                 NicheAxis, Supply, NicheAxis, _cellareas, _derivedfrom,
                 cancel, landcoverclass
using EcoSISTEM: SourceSpec
using EcoSISTEM: compress_landcover, layerinfo, layerunit, sourcecrs
# `CRUTS` only: this extension supplies `read(::Type{CRUTS}, ...)`, and `ERA`'s readers live in the
# parent's `erareaders.jl` because they need netCDF rather than a catalogued raster source.
using EcoSISTEM: CRUTS
# The ONLY remaining reference to the submodule, and it is confined to `deprecations.jl`: five
# deprecated per-source wrappers (`Worldclim_bioclim` and friends) plus `readworldclim` and
# `readCHELSA_monthly` are declared in `ClimatePref/deprecations.jl`, which is where the deprecations
# stay, so their methods have to be attached there. Everything else this extension touches is the
# parent's.
using EcoSISTEM: ClimatePref
# The dataset-reading pipeline moved from `ClimatePref` to `src/datasetread.jl` in v0.5.0, so these
# are imported from the parent now. `_documentedceiling` and `_isblankcrs` were always the parent's.
using EcoSISTEM: _defaultfn, _defaultscale, _documentedceiling,
                 _filelist, _firstfile, _getrasterkw, _isblankcrs,
                 _readraw, _readmonthlydir, _readsource,
                 _rescalepublished, _stackaxis

using DimensionalData
using DimensionalData: dims, At
using Rasters
using Rasters: Raster, X, Y, Ti
using RecipesBase
using SimpleTraits
using Statistics: mean
using Unitful
using Unitful.DefaultSymbols

import Base: read

# Point `RasterDataSources` at its own subdirectory of EcoSISTEM's scratch space, keeping downloads
# under that scratch lifecycle - unless the user has already set `RASTERDATASOURCES_PATH`.
#
# **Here rather than in the parent's `__init__`**: the variable configures `RasterDataSources`
# itself, so it is meaningless until that package is loaded, and setting it from the parent would be
# the one remaining reason for the parent to name that package at all.
function __init__()
    get!(ENV, "RASTERDATASOURCES_PATH") do
        return EcoSISTEM.assetdir(RasterDataSources)
    end
    return nothing
end

include("hooks.jl")
include("read.jl")
include("landcover.jl")
include("supply.jl")
include("recipes.jl")
include("deprecations.jl")

end
