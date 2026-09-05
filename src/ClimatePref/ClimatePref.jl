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
# than inherited from a neighbour's imports. A geographic region is an `Extents.Extent`:
# `boundingbox` returns one and `cut` takes one.
using DimensionalData
using Rasters
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Dates: Dates
# The bare module name, so `deprecations.jl` can write `EcoSISTEM.convert_coords`; and the three
# spec types this module re-exports.
using EcoSISTEM
using EcoSISTEM: SourceSpec, ShapeSpec, ConstructedRasterSpec
import Extents

using EcoSISTEM: _crsunit, _isangle, _stacklayers, _isblankcrs
# The shared "which index contains this coordinate" helper, from `src/cellgeometry.jl`.
using EcoSISTEM: _axisindexat
# Three pieces of the dataset-reading pipeline (`src/datasetread.jl`), called by the `ERA`/`CERA`
# readers and by `boundingbox` below.
using EcoSISTEM: _applycut, _locus, _rastertodimarray
using EcoSISTEM: boundingbox
using EcoSISTEM: ERA, CERA, CRUTS, _istimeaxis
using EcoSISTEM: extract_values
# The phylogenetic trait models: declared in `src/extensions.jl`, implemented in
# `ext/EcoSISTEMPhyloExt`.
using EcoSISTEM: Brownian, varcovar, fitbrownian
using EcoSISTEM: NicheAxis, TemperatureAxis
using EcoSISTEM: AbstractClimate, ClimateRaster, CODE_TYPE,
                 ConstructedRasterSpec,
                 ShapeSpec, AbstractCombineStage, CombineOnTargetGrid,
                 CombineOnSourceGrid, EcoSISTEMSource, SyntheticData,
                 DerivedData, IsRasterData, in_memory_raster
using EcoSISTEM.Units: _monthindex
using RecipesBase
using SimpleTraits
import Base: size, length, eltype
# `import`, not `using`: `deprecations.jl` below **extends** `readfile` with the deprecated
# positional-extent method, and a `using`-imported name cannot be extended.
import EcoSISTEM: readfile
# Declared in `src/extensions.jl` with their `EcoSISTEMRasterDataSourcesExt` methods, so that the
# whole extension's surface sits in one place. Re-exported below.
using EcoSISTEM: compress_landcover, sourcecrs
using EcoSISTEM: LayerRecord, AxisNode, AbstractAccumulationPeriod,
                 ConstantAccumulationPeriod, PerSliceAccumulationPeriod,
                 PerCellAccumulationPeriod
using EcoSISTEM: _CATEGORIES, _OPTIONAL_COLUMNS, _PERIOD_CATEGORIES,
                 _PERIOD_INHERITED_CATEGORIES,
                 _PERSLICE_PERIODS, _REQUIRED_COLUMNS, _SUPERSCRIPTS,
                 _VALUETYPES,
                 _axischain, _axisnames, _axisnode, _catalogue,
                 _checkaxishomogeneity, _checkaxisunit, _checkperiod,
                 _checkrangerows,
                 _checkschema, _checkunitdimension, _checkupstreamscale,
                 _datasettype,
                 _documentedceiling, _layerfile, _layerrow, _layertable,
                 _leafaxes, _normdim, _parsecategory, _parseperiod,
                 _parsepublishedscale, _parsetemporal, _parsevaluetype,
                 _periodcode,
                 _perioddivisors, _periodphrase, _periodsagree, _periodwording,
                 _persliceperiod, _publishedscale, _readdivisors,
                 _refreshaxisnames!,
                 _resolveaxis, _showtree, _stridedsample, _unitstr,
                 layeraxes, layeraxis, layerinfo, layerrate,
                 layersbyaxis, layerunit, partition_eq, _sharedunit, _sharedaxis
export layerunit, layeraxis
public layerinfo, layersbyaxis, layeraxes, LayerRecord, AxisNode
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
