@warn "This functionality remains under development!"

include("ClimateTypes.jl")
export Worldclim_monthly, Worldclim_bioclim, ERA, CERA, CRUTS, CHELSA_bioclim, CHELSA_monthly, Reference, Landcover

include("ReadData.jl")
export read, searchdir, readworldclim, readbioclim, readERA,
readCERA, readfile, readCHELSA_monthly, readCHELSA_bioclim, readCRUTS,
readMet, readlc

include("ExtractClimate.jl")
export extractvalues

include("DataCleaning.jl")
export create_reference, upresolution, downresolution, downresolution!, compressLC

include("Plotting.jl")
export getprofile

include("PhyloModels.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

