# SPDX-License-Identifier: LGPL-3.0-or-later

module ClimatePref

# ERA extension
function retrieve_era5 end
export retrieve_era5

include("ClimateTypes.jl")
export ClimateRaster, Worldclim_monthly, ERA, CERA, CRUTS,
       CHELSA_monthly, Reference

include("ReadData.jl")
export searchdir, readworldclim, readERA,
       readCERA, readfile, readCHELSA_monthly, readCRUTS

include("ExtractClimate.jl")
export extractvalues

include("DataCleaning.jl")
export create_reference, upresolution, downresolution, downresolution!,
       compressLC

include("PhyloModels.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

end
