# SPDX-License-Identifier: LGPL-3.0-or-later

module ClimatePref

# ERA extension
function retrieve_era5 end
export retrieve_era5

include("ClimateTypes.jl")
export ClimateRaster, ERA, CERA, CRUTS, Reference
# Deprecated climate type names, retained for backward compatibility.
export Worldclim_bioclim, CHELSA_bioclim, Landcover,
       Worldclim_monthly, CHELSA_monthly

include("LayerUnits.jl")
export layerunit, layeraxis
public layerinfo, layersbyaxis, layeraxes, LayerRecord, AxisNode

include("ReadData.jl")
export searchdir, readworldclim, readERA,
       readCERA, readfile, readCHELSA_monthly, readCRUTS
public boundingbox

include("ExtractClimate.jl")
export extractvalues

include("DataCleaning.jl")
export create_reference, upresolution, downresolution, downresolution!,
       compressLandCover

include("PhyloModels.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

# Deprecated climate API (constructors + readers), collected last so the types/readers it forwards to are
# defined. The `@deprecate`d names are (re-)exported by the macro; `readworldclim`/`readfile` are exported
# above with the live readers.
include("deprecations.jl")

end
