module ClimatePref

# ERA extension
function retrieve_era5 end
export retrieve_era5

include("ClimateTypes.jl")
export Worldclim_monthly, Worldclim_bioclim, ERA, CERA, CRUTS, CHELSA_bioclim,
       CHELSA_monthly, Reference, Landcover

include("ReadData.jl")
export read, searchdir, readworldclim, readbioclim, readERA,
       readCERA, readfile, readCHELSA_monthly, readCHELSA_bioclim, readCRUTS,
       readlc

include("ExtractClimate.jl")
export extractvalues

include("DataCleaning.jl")
export create_reference, upresolution, downresolution, downresolution!,
       compressLC

include("PhyloModels.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

end
