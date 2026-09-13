# SPDX-License-Identifier: LGPL-3.0-or-later

module ClimatePref

export ClimateRaster, ERA, CERA, CRUTS
# Deprecated climate type names, retained for backward compatibility.
export Worldclim_bioclim, CHELSA_bioclim, Landcover,
       Worldclim_monthly, CHELSA_monthly

# The parent module owns most of what this one exports. The catalogue, the raster geometry, the
# dataset-reading pipeline, the geographic vocabulary and the phylogenetic trait models are all
# dataset-agnostic machinery that types such as `SourceSpec` need, so none of it can sit in a
# submodule included after them. Each is imported back below, because an `export` or `public` line
# here can only name something this module can see.
export SourceSpec, ShapeSpec, ConstructedRasterSpec
public in_memory_raster
export CombineOnTargetGrid, CombineOnSourceGrid
public AbstractCombineStage, AbstractClimate, EcoSISTEMSource
public CODE_TYPE

# A `using` is module-scoped, so every name any file in this submodule needs is stated here rather
# than inherited from a neighbour's imports. Only the deprecations live here now, so what is
# imported is what they forward to, plus every name this module re-exports or declares public -
# an `export` or `public` line can only name something the module can see.
using DimensionalData
using Unitful
# The bare module name, so `deprecations.jl` can write `EcoSISTEM.convert_coords`; a geographic
# region is an `Extents.Extent`, which the deprecated positional-extent readers build.
using EcoSISTEM
import Extents

# Re-exported: the spec types, the raster types and their sources, the combine stages.
using EcoSISTEM: SourceSpec, ShapeSpec, ConstructedRasterSpec, in_memory_raster
using EcoSISTEM: AbstractClimate, ClimateRaster, CODE_TYPE, ERA, CERA, CRUTS,
                 EcoSISTEMSource, AbstractCombineStage, CombineOnTargetGrid,
                 CombineOnSourceGrid
# Re-exported: the region vocabulary and the phylogenetic trait models (declared in
# `src/extensions.jl`, implemented in `ext/EcoSISTEMPhyloExt`).
using EcoSISTEM: boundingbox, extract_values
using EcoSISTEM: Brownian, varcovar, fitbrownian
# `import`, not `using`: `deprecations.jl` below **extends** `readfile` with the deprecated
# positional-extent method, and a `using`-imported name cannot be extended.
import EcoSISTEM: readfile
# Re-exported: the two hooks declared in `src/extensions.jl` with their
# `EcoSISTEMRasterDataSourcesExt` methods, and the catalogue's public surface.
using EcoSISTEM: compress_landcover, sourcecrs
using EcoSISTEM: LayerRecord, AxisNode, DatasetRecord,
                 AbstractAccumulationPeriod,
                 ConstantAccumulationPeriod, PerSliceAccumulationPeriod,
                 PerCellAccumulationPeriod
using EcoSISTEM: datasetinfo, layeraxes, layeraxis, layerinfo, layerrate,
                 layersbyaxis, layerunit
export layerunit, layeraxis
public layerinfo, layersbyaxis, layeraxes, LayerRecord, AxisNode
public datasetinfo, DatasetRecord
public layerrate
public AbstractAccumulationPeriod, ConstantAccumulationPeriod,
       PerSliceAccumulationPeriod, PerCellAccumulationPeriod

export readworldclim, readERA,
       readCERA, readfile, readCHELSA_monthly, readCRUTS
public boundingbox, sourcecrs

public extract_values

export compress_landcover

export Brownian, fitbrownian, varcovar

# The deprecated climate API - constructors and readers - collected last so that everything it
# forwards to is already defined. The dataset-typed subset (the five per-source wrapper constructors,
# `readworldclim` and `readCHELSA_monthly`) is in `EcoSISTEMRasterDataSourcesExt` instead, since none
# of those can be defined without `RasterDataSources`. The names are exported here either way.
include("deprecations.jl")

end
