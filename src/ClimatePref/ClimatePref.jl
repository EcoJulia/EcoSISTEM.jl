using Requires
function __init__()
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
        println("Creating ECMWF interface ...")
        include("ERA_interim_tools.jl")
        export retrieve_era_interim
        include("ECMWF_tools.jl")
        export retrieve_ECMWF
    end
end

include("ClimateTypes.jl")
export Worldclim, Bioclim, ERA, CERA, CRUTS, CHELSA_bioclim, CHELSA_monthly, Reference

include("ReadData.jl")
export read, searchdir, readworldclim, readbioclim, readERA, 
readCERA, readfile, readCHELSA_monthly, readCHELSA_bioclim

include("ReadGBIF.jl")
export ReadGBIF

include("ReadTPL.jl")
export ReadTPL

include("ExtractClimate.jl")
export extractvalues

include("DataCleaning.jl")
export create_reference, gardenmask, genus_worldclim_average,
    genus_worldclim_monthly, upresolution, downresolution

include("Conversion.jl")
export worldclim_to_DB, era_to_DB, CHELSA_to_DB

include("Plotting.jl")
export getprofile

include("PhyloModels.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

include("SimpleSDMInterface.jl")
export Worldclim, Bioclim, CHELSA_bioclim, Landcover